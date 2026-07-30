//! Top-level ingress, parsing, and output orchestration.

use std::{
    any::Any,
    fs::File,
    marker::PhantomData,
    panic::{AssertUnwindSafe, catch_unwind},
    path::PathBuf,
    time::Instant,
};

use needletail::{
    errors::{ParseError, ParseErrorKind},
    parse_fastx_reader,
    parser::SequenceRecord,
};

use crate::{
    adapter::{AdapterPreset, TrimAdaptersTransform},
    cli::{Cli, Ingress, InvalidFastqPolicy, UiPolicy},
    ena::{Accession, EnaClient, EnaInput},
    error::{
        EnaContentProblem, IndeterminateInputError, InternalError, IoError, MalformedInputError,
        Result, RunError, UnavailableInputError, UsageError,
    },
    filter::{MaxNsFilter, MinEntropyFilter, MinLengthFilter, MinMeanQualityFilter},
    output::{OutputArgs, PairedOutputHandle, SingleOutputHandle, UnitOutput},
    pair_merge::MergePairsTransform,
    plan::{
        BuildPlan, Execute, Execution, Logical, OrphanPolicy, Plan, RecordPair, TransformArena,
    },
    progress::ProgressReporter,
    quality::{QualityBinCount, QualityBinTransform, QualityTrimTransform},
    record::{InputSource, InvalidFastqReport, MateSide, ReadStats, RecordProvenance, RecordView},
    report::{self, RunContext as RunSummaryContext, RunLayout},
};

struct SingleEnd;

struct PairedEnd;

trait FastqRunLayout {
    type Source: RunSource;
    type Readers;
    type Output;
}

impl FastqRunLayout for SingleEnd {
    type Source = SingleSource;
    type Readers = Box<dyn std::io::Read + Send>;
    type Output = SingleOutputHandle;
}

impl FastqRunLayout for PairedEnd {
    type Source = PairedSource;
    type Readers = PairedReaders;
    type Output = PairedOutputHandle;
}

trait RunSource {
    fn input_label(&self) -> String;
    fn input_origin(&self, mate: &'static str) -> InputOrigin<'_>;
    fn summary_context(&self) -> RunSummaryContext;
}

#[derive(Clone, Copy)]
enum InputOrigin<'source> {
    Ena(&'source Accession),
    Local(&'source std::path::Path),
}

enum SingleSource {
    Ena { accession: Accession },
    Local { input: PathBuf },
}

enum PairedSource {
    Ena { accession: Accession },
    LocalInterleaved { input: PathBuf },
    LocalSplit { input1: PathBuf, input2: PathBuf },
}

enum PairedReaders {
    Interleaved(Box<dyn std::io::Read + Send>),
    Split {
        left: Box<dyn std::io::Read + Send>,
        right: Box<dyn std::io::Read + Send>,
    },
}

impl RunSource for SingleSource {
    fn input_label(&self) -> String {
        match self {
            Self::Ena { accession } => format!("ena:{accession}"),
            Self::Local { input } => format!("local:{}", input.display()),
        }
    }

    fn input_origin(&self, _mate: &'static str) -> InputOrigin<'_> {
        match self {
            Self::Ena { accession } => InputOrigin::Ena(accession),
            Self::Local { input } => InputOrigin::Local(input),
        }
    }

    fn summary_context(&self) -> RunSummaryContext {
        match self {
            Self::Ena { accession } => RunSummaryContext {
                ingress_mode: report::IngressMode::Ena,
                layout: RunLayout::Single,
                accession: Some(accession.to_string()),
                input1: None,
                input2: None,
            },
            Self::Local { input } => RunSummaryContext {
                ingress_mode: report::IngressMode::Local,
                layout: RunLayout::Single,
                accession: None,
                input1: Some(input.display().to_string()),
                input2: None,
            },
        }
    }
}

impl SingleSource {
    fn provenance(&self) -> RecordProvenance<'_> {
        match self {
            Self::Ena { accession } => RecordProvenance {
                source: InputSource::Ena {
                    accession: accession.as_str(),
                },
                mate: None,
            },
            Self::Local { input } => RecordProvenance {
                source: InputSource::LocalSingle { input },
                mate: None,
            },
        }
    }
}

impl RunSource for PairedSource {
    fn input_label(&self) -> String {
        match self {
            Self::Ena { accession } => format!("ena:{accession}"),
            Self::LocalInterleaved { input } => format!("local-interleaved:{}", input.display()),
            Self::LocalSplit { input1, input2 } => {
                format!("local-paired:{}|{}", input1.display(), input2.display())
            }
        }
    }

    fn input_origin(&self, mate: &'static str) -> InputOrigin<'_> {
        match self {
            Self::Ena { accession } => InputOrigin::Ena(accession),
            Self::LocalInterleaved { input } => InputOrigin::Local(input),
            Self::LocalSplit { input2, .. } if mate == "right" => InputOrigin::Local(input2),
            Self::LocalSplit { input1, .. } => InputOrigin::Local(input1),
        }
    }

    fn summary_context(&self) -> RunSummaryContext {
        match self {
            Self::Ena { accession } => RunSummaryContext {
                ingress_mode: report::IngressMode::Ena,
                layout: RunLayout::Paired,
                accession: Some(accession.to_string()),
                input1: None,
                input2: None,
            },
            Self::LocalInterleaved { input } => RunSummaryContext {
                ingress_mode: report::IngressMode::Local,
                layout: RunLayout::Paired,
                accession: None,
                input1: Some(input.display().to_string()),
                input2: None,
            },
            Self::LocalSplit { input1, input2 } => RunSummaryContext {
                ingress_mode: report::IngressMode::Local,
                layout: RunLayout::Paired,
                accession: None,
                input1: Some(input1.display().to_string()),
                input2: Some(input2.display().to_string()),
            },
        }
    }
}

impl PairedSource {
    fn provenance(&self, mate: MateSide) -> RecordProvenance<'_> {
        match self {
            Self::Ena { accession } => RecordProvenance {
                source: InputSource::Ena {
                    accession: accession.as_str(),
                },
                mate: Some(mate),
            },
            Self::LocalInterleaved { input } => RecordProvenance {
                source: InputSource::LocalInterleavedPaired { input },
                mate: Some(mate),
            },
            Self::LocalSplit { input1, input2 } => RecordProvenance {
                source: InputSource::LocalPaired { input1, input2 },
                mate: Some(mate),
            },
        }
    }
}

struct RunContext<L: FastqRunLayout> {
    source: L::Source,
    readers: L::Readers,
    output: L::Output,
    _layout: PhantomData<L>,
}

type SingleEndContext = RunContext<SingleEnd>;
type PairedEndContext = RunContext<PairedEnd>;

struct RunConfig {
    min_length: usize,
    max_ns: usize,
    min_mean_q: f64,
    min_entropy: f64,
    trim_min_q: u8,
    bin_qualities: Option<QualityBinCount>,
    adapter_preset: AdapterPreset,
    merge_pairs: bool,
    passthrough: bool,
    merge_min_overlap: usize,
    merge_max_mismatch_rate: f32,
    merge_min_correction_delta_q: u8,
    invalid_fastq_policy: InvalidFastqPolicy,
    progress_every: u64,
    summary: Option<PathBuf>,
    invalid_fastq_report: Option<PathBuf>,
}

impl From<&Cli> for RunConfig {
    fn from(cli: &Cli) -> Self {
        Self {
            min_length: cli.min_length,
            max_ns: cli.max_ns,
            min_mean_q: cli.min_mean_q,
            min_entropy: cli.min_entropy,
            trim_min_q: cli.trim_min_q,
            bin_qualities: cli.bin_qualities,
            adapter_preset: cli.adapter_preset,
            merge_pairs: cli.merge_pairs,
            passthrough: cli.passthrough,
            merge_min_overlap: cli.merge_min_overlap,
            merge_max_mismatch_rate: cli.merge_max_mismatch_rate,
            merge_min_correction_delta_q: cli.merge_min_correction_delta_q,
            invalid_fastq_policy: cli.invalid_fastq_policy,
            progress_every: cli.progress_every,
            summary: cli.summary.clone(),
            invalid_fastq_report: cli.invalid_fastq_report.clone(),
        }
    }
}

impl RunConfig {
    fn validate_layout(&self, layout: RunLayout) -> Result<()> {
        if self.merge_pairs && layout == RunLayout::Single {
            return Err(UsageError::MergeRequiresPairedInput.into());
        }

        Ok(())
    }

    fn build_plan(&self, layout: RunLayout) -> Result<Plan<Execution>> {
        self.validate_layout(layout)?;

        let plan = Plan::<Logical>::new();
        if self.passthrough {
            return Ok(plan.orphan_policy(OrphanPolicy::DropPair).compile());
        }

        let plan = if self.merge_pairs && layout == RunLayout::Paired {
            plan.merge_pairs(crate::pair_merge::MergePairsConfig {
                min_overlap: self.merge_min_overlap,
                max_mismatch_rate: self.merge_max_mismatch_rate,
                min_correction_delta_q: self.merge_min_correction_delta_q,
            })?
        } else {
            plan
        };
        let plan = plan.max_ns(self.max_ns);
        let plan = match self.adapter_preset.catalog() {
            Some(catalog) => plan.trim_adapters(catalog),
            None => plan,
        };

        let plan = plan
            .quality_trim(self.trim_min_q)
            .min_length(self.min_length)
            .min_mean_q(self.min_mean_q)
            .min_entropy(self.min_entropy);

        let plan = match self.bin_qualities {
            Some(count) => plan.quality_bin(count),
            None => plan,
        };

        Ok(plan.orphan_policy(OrphanPolicy::DropPair).compile())
    }
}

/// Run the CLI-selected ingress path to completion.
///
/// # Errors
///
/// Returns an error when ingress resolution, reader construction, parsing, writing, or output
/// finalization fails.
pub(crate) fn run(cli: &Cli) -> Result<()> {
    cli.init_tracing()?;
    tracing::info!(version = env!("CARGO_PKG_VERSION"), "starting nuclease");
    let config = RunConfig::from(cli);
    let ui = cli.ui_policy();

    match cli.ingress()? {
        Ingress::LocalSingle { fastq } => {
            config.validate_layout(RunLayout::Single)?;
            SingleEndContext::open_local(fastq, cli.output_args())?.run(&config, &ui)
        }
        Ingress::LocalInterleavedPaired { fastq } => {
            config.validate_layout(RunLayout::Paired)?;
            PairedEndContext::open_local_interleaved(fastq, cli.output_args())?.run(&config, &ui)
        }
        Ingress::LocalSplitPaired { r1, r2 } => {
            config.validate_layout(RunLayout::Paired)?;
            PairedEndContext::open_local_split(r1, r2, cli.output_args())?.run(&config, &ui)
        }
        Ingress::Ena { accession } => run_ena(&accession, cli.output_args(), &config, &ui),
    }
}

fn run_ena(
    accession: &Accession,
    output_args: OutputArgs,
    config: &RunConfig,
    ui: &UiPolicy,
) -> Result<()> {
    let client = EnaClient::new()?;
    let input = client.resolve(accession)?;

    match input {
        EnaInput::Single(fastq) => {
            config.validate_layout(RunLayout::Single)?;
            let stream = client.stream(fastq)?;
            SingleEndContext::open_ena(accession, stream, output_args)?.run(config, ui)
        }
        EnaInput::Paired { left, right } => {
            config.validate_layout(RunLayout::Paired)?;
            let left = client.stream(left)?;
            let right = client.stream(right)?;
            PairedEndContext::open_ena(accession, left, right, output_args)?.run(config, ui)
        }
    }
}

impl RunContext<SingleEnd> {
    fn open_local(fastq: PathBuf, output_args: OutputArgs) -> Result<Self> {
        let reader =
            File::open(&fastq).map_err(|source| UnavailableInputError::OpenLocalFastq {
                path: fastq.clone(),
                source,
            })?;
        let output = output_args.resolve_single()?;
        Ok(Self {
            source: SingleSource::Local { input: fastq },
            readers: Box::new(reader),
            output,
            _layout: PhantomData,
        })
    }

    fn open_ena(
        accession: &Accession,
        reader: impl std::io::Read + Send + 'static,
        output_args: OutputArgs,
    ) -> Result<Self> {
        let output = output_args.resolve_single()?;
        Ok(Self {
            source: SingleSource::Ena {
                accession: accession.clone(),
            },
            readers: Box::new(reader),
            output,
            _layout: PhantomData,
        })
    }

    fn run(self, config: &RunConfig, ui: &UiPolicy) -> Result<()> {
        let Self {
            source,
            readers: reader,
            mut output,
            _layout,
        } = self;
        let mut plan = config.build_plan(RunLayout::Single)?;

        let mut arena = TransformArena::new();
        let mut stats = read_stats(config)?;
        // build out the mutable state needed to run the application loop
        let mut parser = parse_fastx_reader(reader).map_err(|source_error| {
            parser_error(
                source.input_origin("single"),
                source.input_label(),
                "single",
                config.invalid_fastq_policy,
                &stats,
                source_error,
            )
        })?;
        let mut progress = ProgressReporter::new(ui.progress_mode, config.progress_every);
        let started_at = Instant::now();
        let admission = FastqAdmission::<SingleEnd>::new(&source, config.invalid_fastq_policy);

        while let Some(next_record) =
            catch_parser_panic(&source, "single", &stats, || parser.next())?
        {
            let parsed_record = admission.parse("single", &mut stats, next_record)?;
            let Some(record) = admission.single(&parsed_record, &mut stats)? else {
                continue;
            };
            arena.reset();

            let outcome = plan.execute(record, &mut arena, &mut stats)?;
            output.write_outcome(&outcome, &mut stats)?;

            progress.maybe_report(&stats);
        }

        progress.finish();
        stats.finish_invalid_fastq_report()?;
        output.finish()?;
        let summary =
            report::RunSummary::from_stats(source.summary_context(), &stats, started_at.elapsed());
        if ui.show_summary {
            report::print_summary(&summary);
        }
        if let Some(path) = &config.summary {
            report::write_summary_json(path, &summary)?;
        }
        Ok(())
    }
}

impl RunContext<PairedEnd> {
    fn open_local_interleaved(fastq: PathBuf, output_args: OutputArgs) -> Result<Self> {
        let reader =
            File::open(&fastq).map_err(|source| UnavailableInputError::OpenLocalFastq {
                path: fastq.clone(),
                source,
            })?;
        let output = output_args.resolve_paired()?;
        Ok(Self {
            source: PairedSource::LocalInterleaved { input: fastq },
            readers: PairedReaders::Interleaved(Box::new(reader)),
            output,
            _layout: PhantomData,
        })
    }

    fn open_local_split(r1: PathBuf, r2: PathBuf, output_args: OutputArgs) -> Result<Self> {
        let reader1 = File::open(&r1).map_err(|source| UnavailableInputError::OpenLocalFastq {
            path: r1.clone(),
            source,
        })?;
        let reader2 = File::open(&r2).map_err(|source| UnavailableInputError::OpenLocalFastq {
            path: r2.clone(),
            source,
        })?;
        let output = output_args.resolve_paired()?;
        Ok(Self {
            source: PairedSource::LocalSplit {
                input1: r1,
                input2: r2,
            },
            readers: PairedReaders::Split {
                left: Box::new(reader1),
                right: Box::new(reader2),
            },
            output,
            _layout: PhantomData,
        })
    }

    fn open_ena(
        accession: &Accession,
        r1: impl std::io::Read + Send + 'static,
        r2: impl std::io::Read + Send + 'static,
        output_args: OutputArgs,
    ) -> Result<Self> {
        let output = output_args.resolve_paired()?;
        Ok(Self {
            source: PairedSource::Ena {
                accession: accession.clone(),
            },
            readers: PairedReaders::Split {
                left: Box::new(r1),
                right: Box::new(r2),
            },
            output,
            _layout: PhantomData,
        })
    }

    #[allow(
        clippy::too_many_lines,
        reason = "paired run orchestration intentionally keeps split and interleaved parser loops local to RunContext"
    )]
    fn run(self, config: &RunConfig, ui: &UiPolicy) -> Result<()> {
        let Self {
            source,
            readers,
            mut output,
            _layout,
        } = self;
        let mut plan = config.build_plan(RunLayout::Paired)?;
        let mut arena = TransformArena::new();
        let mut stats = read_stats(config)?;
        let mut progress = ProgressReporter::new(ui.progress_mode, config.progress_every);
        let started_at = Instant::now();
        let admission = FastqAdmission::<PairedEnd>::new(&source, config.invalid_fastq_policy);

        match readers {
            PairedReaders::Split { left, right } => {
                let mut parser_r1 = parse_fastx_reader(left).map_err(|source_error| {
                    parser_error(
                        source.input_origin("left"),
                        source.input_label(),
                        "left",
                        config.invalid_fastq_policy,
                        &stats,
                        source_error,
                    )
                })?;
                let mut parser_r2 = parse_fastx_reader(right).map_err(|source_error| {
                    parser_error(
                        source.input_origin("right"),
                        source.input_label(),
                        "right",
                        config.invalid_fastq_policy,
                        &stats,
                        source_error,
                    )
                })?;

                loop {
                    let next_r1 = catch_parser_panic(&source, "left", &stats, || parser_r1.next())?;
                    let next_r2 =
                        catch_parser_panic(&source, "right", &stats, || parser_r2.next())?;

                    match (next_r1, next_r2) {
                        (Some(record_r1), Some(record_r2)) => {
                            let parsed_r1 = admission.parse("left", &mut stats, record_r1)?;
                            let parsed_r2 = admission.parse("right", &mut stats, record_r2)?;
                            let Some(pair) = admission.pair(&parsed_r1, &parsed_r2, &mut stats)?
                            else {
                                continue;
                            };

                            arena.reset();
                            let outcome = plan.execute(pair, &mut arena, &mut stats)?;
                            output.write_outcome(&outcome, &mut stats)?;
                            progress.maybe_report(&stats);
                        }
                        (None, None) => break,
                        _ => return Err(record_count_error(&source, &stats)),
                    }
                }
            }
            PairedReaders::Interleaved(reader) => {
                let mut parser = parse_fastx_reader(reader).map_err(|source_error| {
                    parser_error(
                        source.input_origin("left"),
                        source.input_label(),
                        "left",
                        config.invalid_fastq_policy,
                        &stats,
                        source_error,
                    )
                })?;
                let mut left_buffer = InterleavedLeftBuffer::default();

                loop {
                    let next_left = catch_parser_panic(&source, "left", &stats, || parser.next())?;
                    let Some(left_record) = next_left else {
                        break;
                    };

                    let parsed_left = admission.parse("left", &mut stats, left_record)?;
                    let left =
                        admission.buffered_left_record(&parsed_left, &stats, &mut left_buffer)?;

                    let next_right =
                        catch_parser_panic(&source, "right", &stats, || parser.next())?;
                    let Some(right_record) = next_right else {
                        return Err(record_count_error(&source, &stats));
                    };
                    let parsed_right = admission.parse("right", &mut stats, right_record)?;
                    let right =
                        admission.paired_record(&parsed_right, MateSide::Right, "right", &stats)?;

                    let Some(pair) = admission.admit_pair(left, right, &mut stats)? else {
                        continue;
                    };

                    arena.reset();
                    let outcome = plan.execute(pair, &mut arena, &mut stats)?;
                    output.write_outcome(&outcome, &mut stats)?;
                    progress.maybe_report(&stats);
                }
            }
        }

        progress.finish();
        stats.finish_invalid_fastq_report()?;
        output.finish()?;
        let summary =
            report::RunSummary::from_stats(source.summary_context(), &stats, started_at.elapsed());
        if ui.show_summary {
            report::print_summary(&summary);
        }
        if let Some(path) = &config.summary {
            report::write_summary_json(path, &summary)?;
        }
        Ok(())
    }
}

#[derive(Default)]
struct InterleavedLeftBuffer {
    header: Vec<u8>,
    sequence: Vec<u8>,
    quality: Vec<u8>,
}

impl InterleavedLeftBuffer {
    fn copy_from<'buffer>(
        &'buffer mut self,
        header: &[u8],
        sequence: &[u8],
        quality: &[u8],
    ) -> RecordView<'buffer> {
        self.header.clear();
        self.sequence.clear();
        self.quality.clear();

        self.header.extend_from_slice(header);
        self.sequence.extend_from_slice(sequence);
        self.quality.extend_from_slice(quality);

        RecordView::new(&self.header, &self.sequence, &self.quality)
    }
}

struct FastqAdmission<'source, L: FastqRunLayout> {
    source: &'source L::Source,
    policy: InvalidFastqPolicy,
    _layout: PhantomData<L>,
}

impl<'source, L: FastqRunLayout> FastqAdmission<'source, L> {
    fn new(source: &'source L::Source, policy: InvalidFastqPolicy) -> Self {
        Self {
            source,
            policy,
            _layout: PhantomData,
        }
    }

    fn parse<'record>(
        &self,
        mate: &'static str,
        stats: &mut ReadStats,
        next_record: std::result::Result<SequenceRecord<'record>, ParseError>,
    ) -> Result<SequenceRecord<'record>> {
        match next_record {
            Ok(record) => Ok(record),
            Err(error) => self.parser_error(mate, stats, error),
        }
    }

    fn parser_error<T>(
        &self,
        mate: &'static str,
        stats: &mut ReadStats,
        error: ParseError,
    ) -> Result<T> {
        let source = self.source.input_label();
        let parser_error_kind = format!("{:?}", error.kind);
        let parser_error_message = error.to_string();
        let parser_error_line = (error.position.line > 0).then_some(error.position.line);

        stats.record_invalid_parse_error(self.policy, |context| {
            context.parse_error(
                &source,
                mate,
                parser_error_kind.clone(),
                parser_error_message.clone(),
                parser_error_line,
            )
        })?;

        if self.policy == InvalidFastqPolicy::WarnDrop {
            tracing::warn!(
                source,
                mate,
                parser_error_kind,
                parser_error = parser_error_message,
                "invalid FASTQ parser error is unrecoverable; stopping instead of dropping and continuing"
            );
        }

        Err(parser_error(
            self.source.input_origin(mate),
            source,
            mate,
            self.policy,
            stats,
            error,
        ))
    }
}

impl FastqAdmission<'_, SingleEnd> {
    fn single<'record>(
        &'record self,
        parsed_record: &'record SequenceRecord<'_>,
        stats: &mut ReadStats,
    ) -> Result<Option<RecordView<'record>>> {
        let sequence = parsed_record.raw_seq();
        let quality = parsed_record
            .qual()
            .ok_or_else(|| missing_quality_error(self.source, "single", stats))?;
        let record = RecordView::new(parsed_record.id(), sequence, quality)
            .with_provenance(self.source.provenance());

        stats.record_seen(sequence.len());

        record.validate(self.policy, stats)
    }
}

impl<'source> FastqAdmission<'source, PairedEnd> {
    fn pair<'record>(
        &'record self,
        parsed_r1: &'record SequenceRecord<'_>,
        parsed_r2: &'record SequenceRecord<'_>,
        stats: &mut ReadStats,
    ) -> Result<Option<RecordPair<'record>>> {
        let left = self.paired_record(parsed_r1, MateSide::Left, "left", stats)?;
        let right = self.paired_record(parsed_r2, MateSide::Right, "right", stats)?;

        self.admit_pair(left, right, stats)
    }

    fn paired_record<'record>(
        &'record self,
        parsed_record: &'record SequenceRecord<'_>,
        mate: MateSide,
        mate_label: &'static str,
        stats: &ReadStats,
    ) -> Result<RecordView<'record>> {
        let sequence = parsed_record.raw_seq();
        let quality = parsed_record
            .qual()
            .ok_or_else(|| missing_quality_error(self.source, mate_label, stats))?;

        Ok(RecordView::new(parsed_record.id(), sequence, quality)
            .with_provenance(self.source.provenance(mate)))
    }

    fn buffered_left_record<'record>(
        &'record self,
        parsed_record: &SequenceRecord<'_>,
        stats: &ReadStats,
        buffer: &'record mut InterleavedLeftBuffer,
    ) -> Result<RecordView<'record>>
    where
        'source: 'record,
    {
        let sequence = parsed_record.raw_seq();
        let quality = parsed_record
            .qual()
            .ok_or_else(|| missing_quality_error(self.source, "left", stats))?;

        Ok(buffer
            .copy_from(parsed_record.id(), sequence, quality)
            .with_provenance(self.source.provenance(MateSide::Left)))
    }

    fn admit_pair<'record>(
        &self,
        left: RecordView<'record>,
        right: RecordView<'record>,
        stats: &mut ReadStats,
    ) -> Result<Option<RecordPair<'record>>> {
        stats.record_seen(left.sequence().len());
        stats.record_seen(right.sequence().len());
        stats.pairs_seen += 1;

        left.validate_pair(right, self.policy, stats)
    }
}

fn catch_parser_panic<T, S: RunSource>(
    source: &S,
    mate: &'static str,
    stats: &ReadStats,
    operation: impl FnOnce() -> T,
) -> Result<T> {
    catch_unwind(AssertUnwindSafe(operation)).map_err(|panic| {
        let panic = panic_message(&panic);
        match source.input_origin(mate) {
            InputOrigin::Ena(accession) => IndeterminateInputError::ParserPanic {
                accession: accession.clone(),
                mate,
                reads_seen: stats.reads_seen,
                pairs_seen: stats.pairs_seen,
                panic,
            }
            .into(),
            InputOrigin::Local(_) => InternalError::LocalParserPanic {
                source_label: source.input_label(),
                mate,
                reads_seen: stats.reads_seen,
                pairs_seen: stats.pairs_seen,
                panic,
            }
            .into(),
        }
    })
}

fn missing_quality_error<S: RunSource>(
    source: &S,
    mate: &'static str,
    stats: &ReadStats,
) -> RunError {
    match source.input_origin(mate) {
        InputOrigin::Ena(accession) => IndeterminateInputError::Content {
            accession: accession.to_string(),
            problem: EnaContentProblem::MissingQuality { mate },
        }
        .into(),
        InputOrigin::Local(_) => MalformedInputError::MissingQuality {
            source_label: source.input_label(),
            mate,
            reads_seen: stats.reads_seen,
            pairs_seen: stats.pairs_seen,
        }
        .into(),
    }
}

fn panic_message(panic: &Box<dyn Any + Send>) -> String {
    if let Some(message) = panic.downcast_ref::<&str>() {
        (*message).to_owned()
    } else if let Some(message) = panic.downcast_ref::<String>() {
        message.clone()
    } else {
        "<non-string panic>".to_owned()
    }
}

fn read_stats(config: &RunConfig) -> Result<ReadStats> {
    let mut stats = ReadStats::default();
    if let Some(path) = &config.invalid_fastq_report {
        stats.set_invalid_fastq_report(InvalidFastqReport::create(path)?);
    }
    Ok(stats)
}

fn parser_error(
    origin: InputOrigin<'_>,
    source_label: String,
    mate: &'static str,
    policy: InvalidFastqPolicy,
    stats: &ReadStats,
    source: ParseError,
) -> RunError {
    match origin {
        InputOrigin::Ena(accession) => IndeterminateInputError::Parser {
            accession: accession.clone(),
            mate: match mate {
                "left" => Some(MateSide::Left),
                "right" => Some(MateSide::Right),
                _ => None,
            },
            reads_seen: stats.reads_seen,
            pairs_seen: stats.pairs_seen,
            source,
        }
        .into(),
        InputOrigin::Local(path) if source.kind == ParseErrorKind::Io => IoError::LocalFastqRead {
            path: path.to_path_buf(),
            source,
        }
        .into(),
        InputOrigin::Local(_) => MalformedInputError::LocalParser {
            source_label,
            mate,
            policy: policy.to_string(),
            reads_seen: stats.reads_seen,
            pairs_seen: stats.pairs_seen,
            invalid_reads: stats.invalid_reads,
            invalid_pairs: stats.invalid_pairs,
            parser_error_kind: format!("{:?}", source.kind),
            source,
        }
        .into(),
    }
}

fn record_count_error(source: &PairedSource, stats: &ReadStats) -> RunError {
    match source {
        PairedSource::Ena { accession } => IndeterminateInputError::Content {
            accession: accession.to_string(),
            problem: EnaContentProblem::RecordCount {
                complete_pairs_seen: stats.pairs_seen,
            },
        }
        .into(),
        PairedSource::LocalInterleaved { .. } => MalformedInputError::InterleavedRecordCount {
            source_label: source.input_label(),
            complete_pairs_seen: stats.pairs_seen,
            reads_seen: stats.reads_seen,
        }
        .into(),
        PairedSource::LocalSplit { .. } => MalformedInputError::PairedRecordCount {
            source_label: source.input_label(),
            complete_pairs_seen: stats.pairs_seen,
            reads_seen: stats.reads_seen,
        }
        .into(),
    }
}

#[cfg(test)]
mod tests {
    use std::{
        error::Error as _,
        fs::File,
        io::{self, Cursor, Write as _},
        marker::PhantomData,
        path::Path,
    };

    use color_eyre::{Result, eyre::bail};
    use needletail::{
        errors::{ErrorPosition, ParseError, ParseErrorKind},
        parse_fastx_reader,
        parser::Format,
    };
    use tempfile::tempdir;

    use crate::{
        adapter::AdapterPreset,
        cli::{Cli, InvalidFastqPolicy},
        ena::Accession,
        error::{
            EnaContentProblem, IndeterminateInputError, IoError, MalformedInputError, RunError,
        },
        output::{
            InterleavedOutput, OutputArgs, OutputEncoding, OutputFormat, PairedRecordOutput,
            SingleOutput, SingleRecordOutput, StreamSink,
        },
        record::{ReadStats, RecordView},
    };

    use super::{
        InputOrigin, PairedEndContext, PairedSource, RunConfig, SingleSource,
        missing_quality_error, parser_error, record_count_error,
    };

    fn single_output_for_vec(format: OutputFormat) -> SingleOutput<StreamSink<Vec<u8>>> {
        SingleOutput::new(StreamSink::new(Vec::new(), format))
    }

    fn interleaved_output_for_vec(format: OutputFormat) -> InterleavedOutput<StreamSink<Vec<u8>>> {
        InterleavedOutput::new(StreamSink::new(Vec::new(), format))
    }

    fn test_cli() -> Cli {
        Cli {
            ena: None,
            input: None,
            paired: false,
            in1: None,
            in2: None,
            min_length: 50,
            max_ns: 4,
            min_mean_q: 20.0,
            trim_min_q: 20,
            bin_qualities: None,
            adapter_preset: AdapterPreset::None,
            merge_pairs: false,
            passthrough: false,
            merge_min_overlap: 10,
            merge_max_mismatch_rate: 0.2,
            merge_min_correction_delta_q: 0,
            min_entropy: 0.0,
            output_format: OutputFormat::Fastq,
            output_encoding: None,
            invalid_fastq_policy: InvalidFastqPolicy::Error,
            out: None,
            out1: None,
            out2: None,
            progress_every: 100_000,
            summary: None,
            invalid_fastq_report: None,
            verbose: 0,
            quiet: 0,
        }
    }

    #[test]
    fn single_end_passthrough_preserves_fastq_bytes() -> Result<()> {
        let temp = tempdir()?;
        let input = temp.path().join("reads.fastq");
        let expected = b"@read1\nACGT\n+\nIIII\n@read2\nTGCA\n+\nJJJJ\n";
        write_fixture(&input, expected)?;

        let reader = File::open(&input)?;
        let output = single_output_for_vec(OutputFormat::Fastq);
        let mut output = output;
        let mut parser = parse_fastx_reader(reader)?;
        while let Some(parsed_record) = parser.next() {
            let parsed_record = parsed_record?;
            let record = RecordView::new(
                parsed_record.id(),
                parsed_record.raw_seq(),
                parsed_record
                    .qual()
                    .expect("FASTQ parser must provide quality scores"),
            );
            output.write_record(record)?;
        }

        assert_eq!(output.into_inner().into_inner(), expected);
        Ok(())
    }

    #[test]
    fn single_end_passthrough_preserves_rich_header_content() -> Result<()> {
        let temp = tempdir()?;
        let input = temp.path().join("reads.fastq");
        let expected = concat!(
            "@read1 sample=alpha lane=3 umi:ACGT-TGCA extra text\n",
            "ACGTN\n",
            "+\n",
            "IIIII\n",
            "@instrument:1:FCID:2:2104:15343:197393 1:N:0:NTTGTA\n",
            "TGCA\n",
            "+\n",
            "!~AB\n"
        )
        .as_bytes();
        write_fixture(&input, expected)?;

        let reader = File::open(&input)?;
        let output = single_output_for_vec(OutputFormat::Fastq);
        let mut output = output;
        let mut parser = parse_fastx_reader(reader)?;
        while let Some(parsed_record) = parser.next() {
            let parsed_record = parsed_record?;
            let record = RecordView::new(
                parsed_record.id(),
                parsed_record.raw_seq(),
                parsed_record
                    .qual()
                    .expect("FASTQ parser must provide quality scores"),
            );
            output.write_record(record)?;
        }

        assert_eq!(output.into_inner().into_inner(), expected);
        Ok(())
    }

    #[test]
    fn single_end_passthrough_can_emit_fasta() -> Result<()> {
        let temp = tempdir()?;
        let input = temp.path().join("reads.fastq");
        write_fixture(
            &input,
            b"@read1 sample=alpha\nACGT\n+\nIIII\n@read2 sample=beta\nTGCA\n+\nJJJJ\n",
        )?;

        let reader = File::open(&input)?;
        let output = single_output_for_vec(OutputFormat::Fasta);
        let mut output = output;
        let mut parser = parse_fastx_reader(reader)?;
        while let Some(parsed_record) = parser.next() {
            let parsed_record = parsed_record?;
            let record = RecordView::new(
                parsed_record.id(),
                parsed_record.raw_seq(),
                parsed_record
                    .qual()
                    .expect("FASTQ parser must provide quality scores"),
            );
            output.write_record(record)?;
        }

        assert_eq!(
            output.into_inner().into_inner(),
            b">read1 sample=alpha\nACGT\n>read2 sample=beta\nTGCA\n"
        );
        Ok(())
    }

    #[test]
    fn paired_passthrough_emits_interleaved_fastq() -> Result<()> {
        let r1 = Cursor::new(b"@r1/1\nAAAA\n+\nIIII\n@r2/1\nCCCC\n+\nJJJJ\n".as_slice());
        let r2 = Cursor::new(b"@r1/2\nTTTT\n+\nKKKK\n@r2/2\nGGGG\n+\nLLLL\n".as_slice());
        let output = interleaved_output_for_vec(OutputFormat::Fastq);
        let mut output = output;
        let mut parser_r1 = parse_fastx_reader(r1)?;
        let mut parser_r2 = parse_fastx_reader(r2)?;
        loop {
            match (parser_r1.next(), parser_r2.next()) {
                (Some(parsed_r1), Some(parsed_r2)) => {
                    let parsed_r1 = parsed_r1?;
                    let parsed_r2 = parsed_r2?;
                    output.write_pair(
                        RecordView::new(
                            parsed_r1.id(),
                            parsed_r1.raw_seq(),
                            parsed_r1
                                .qual()
                                .expect("FASTQ parser must provide quality scores"),
                        ),
                        RecordView::new(
                            parsed_r2.id(),
                            parsed_r2.raw_seq(),
                            parsed_r2
                                .qual()
                                .expect("FASTQ parser must provide quality scores"),
                        ),
                    )?;
                }
                (None, None) => break,
                _ => bail!("paired FASTQ inputs have different record counts"),
            }
        }

        assert_eq!(
            output.into_inner().into_inner(),
            b"@r1/1\nAAAA\n+\nIIII\n@r1/2\nTTTT\n+\nKKKK\n@r2/1\nCCCC\n+\nJJJJ\n@r2/2\nGGGG\n+\nLLLL\n"
        );
        Ok(())
    }

    #[test]
    fn paired_passthrough_fails_when_record_counts_differ() -> Result<()> {
        let temp = tempdir()?;
        let out = temp.path().join("interleaved.fastq");
        let r1 = Cursor::new(b"@r1/1\nAAAA\n+\nIIII\n@r2/1\nCCCC\n+\nJJJJ\n".as_slice());
        let r2 = Cursor::new(b"@r1/2\nTTTT\n+\nKKKK\n".as_slice());

        let output_args = OutputArgs::new(
            OutputFormat::Fastq,
            Some(OutputEncoding::Plain),
            Some(out),
            None,
            None,
        );
        let output = output_args.resolve_paired()?;
        let cli = test_cli();
        let ui = cli.ui_policy();
        let config = RunConfig::from(&cli);
        let error = PairedEndContext {
            source: PairedSource::LocalSplit {
                input1: "reads_1.fastq.gz".into(),
                input2: "reads_2.fastq.gz".into(),
            },
            readers: super::PairedReaders::Split {
                left: Box::new(r1),
                right: Box::new(r2),
            },
            output,
            _layout: PhantomData,
        }
        .run(&config, &ui)
        .expect_err("mismatched paired inputs should fail");

        let RunError::MalformedInput(error) = error else {
            panic!("local count mismatch should be malformed input");
        };
        assert!(matches!(
            *error,
            MalformedInputError::PairedRecordCount {
                complete_pairs_seen: 1,
                reads_seen: 2,
                ..
            }
        ));
        Ok(())
    }

    #[test]
    fn missing_local_input_retains_not_found_source() -> Result<()> {
        let temp = tempdir()?;
        let missing = temp.path().join("missing.fastq");
        let output = OutputArgs::new(
            OutputFormat::Fastq,
            Some(OutputEncoding::Plain),
            None,
            None,
            None,
        );
        let Err(error) = PairedEndContext::open_local_interleaved(missing, output) else {
            panic!("missing local input should fail before output opens");
        };
        assert!(
            error.source().is_some(),
            "typed root should preserve the filesystem source"
        );

        let RunError::UnavailableInput(crate::error::UnavailableInputError::OpenLocalFastq {
            source,
            ..
        }) = error
        else {
            panic!("missing path should retain its typed unavailable-input cause");
        };
        assert_eq!(source.kind(), io::ErrorKind::NotFound);
        Ok(())
    }

    #[test]
    fn parser_errors_are_classified_by_origin_and_typed_kind() -> Result<()> {
        let accession = Accession::new("SRR35939766")?;
        let stats = ReadStats {
            reads_seen: 7,
            pairs_seen: 3,
            ..ReadStats::default()
        };

        let parser_errors = [
            ParseError::from(io::Error::other("remote read failed")),
            ParseError::new_unknown_format(b'!'),
            ParseError::new_invalid_start(b'!', ErrorPosition::default(), Format::Fastq),
            ParseError::new_invalid_separator(b'!', ErrorPosition::default()),
            ParseError::new_unequal_length(4, 1, ErrorPosition::default()),
            ParseError::new_unexpected_end(ErrorPosition::default(), Format::Fastq),
            ParseError::new_empty_file(),
        ];

        for source in parser_errors.iter().cloned() {
            let expected_kind = source.kind.clone();
            let error = parser_error(
                InputOrigin::Ena(&accession),
                format!("ena:{accession}"),
                "single",
                InvalidFastqPolicy::Error,
                &stats,
                source,
            );
            let RunError::IndeterminateInput(IndeterminateInputError::Parser {
                accession: observed_accession,
                mate,
                reads_seen,
                pairs_seen,
                source,
            }) = error
            else {
                panic!("all ENA parser errors should remain indeterminate");
            };
            assert_eq!(observed_accession, accession);
            assert_eq!(mate, None);
            assert_eq!(reads_seen, 7);
            assert_eq!(pairs_seen, 3);
            assert_eq!(source.kind, expected_kind);
        }

        let local_path = Path::new("reads.fastq");
        let error = parser_error(
            InputOrigin::Local(local_path),
            "local:reads.fastq".to_owned(),
            "single",
            InvalidFastqPolicy::Error,
            &stats,
            ParseError::from(io::Error::other("local read failed")),
        );
        let RunError::Io(IoError::LocalFastqRead { source, .. }) = error else {
            panic!("local parser I/O should retain the I/O category");
        };
        assert_eq!(source.kind, ParseErrorKind::Io);

        for source in parser_errors
            .into_iter()
            .filter(|source| source.kind != ParseErrorKind::Io)
        {
            let expected_kind = source.kind.clone();
            let error = parser_error(
                InputOrigin::Local(local_path),
                "local:reads.fastq".to_owned(),
                "single",
                InvalidFastqPolicy::Error,
                &stats,
                source,
            );
            assert!(
                error.source().is_some(),
                "typed root should preserve the parser source"
            );
            let RunError::MalformedInput(error) = error else {
                panic!("local structural parser error should be malformed input");
            };
            assert!(matches!(
                *error,
                MalformedInputError::LocalParser { source, .. }
                    if source.kind == expected_kind
            ));
        }
        Ok(())
    }

    #[test]
    fn record_count_errors_preserve_origin_and_local_layout() -> Result<()> {
        let stats = ReadStats {
            reads_seen: 3,
            pairs_seen: 1,
            ..ReadStats::default()
        };
        let accession = Accession::new("SRR35939766")?;

        assert!(matches!(
            record_count_error(&PairedSource::Ena { accession }, &stats),
            RunError::IndeterminateInput(IndeterminateInputError::Content {
                problem: EnaContentProblem::RecordCount {
                    complete_pairs_seen: 1
                },
                ..
            })
        ));
        assert!(matches!(
            record_count_error(
                &PairedSource::LocalInterleaved {
                    input: "reads.fastq".into(),
                },
                &stats,
            ),
            RunError::MalformedInput(error)
                if matches!(*error, MalformedInputError::InterleavedRecordCount { .. })
        ));
        assert!(matches!(
            record_count_error(
                &PairedSource::LocalSplit {
                    input1: "reads_1.fastq".into(),
                    input2: "reads_2.fastq".into(),
                },
                &stats,
            ),
            RunError::MalformedInput(error)
                if matches!(*error, MalformedInputError::PairedRecordCount { .. })
        ));
        Ok(())
    }

    #[test]
    fn missing_quality_errors_preserve_input_origin() -> Result<()> {
        let stats = ReadStats {
            reads_seen: 7,
            pairs_seen: 3,
            ..ReadStats::default()
        };
        let accession = Accession::new("SRR35939766")?;

        assert!(matches!(
            missing_quality_error(&SingleSource::Ena { accession }, "single", &stats),
            RunError::IndeterminateInput(IndeterminateInputError::Content {
                problem: EnaContentProblem::MissingQuality { mate: "single" },
                ..
            })
        ));
        assert!(matches!(
            missing_quality_error(
                &SingleSource::Local {
                    input: "reads.fastq".into(),
                },
                "single",
                &stats,
            ),
            RunError::MalformedInput(error)
                if matches!(*error, MalformedInputError::MissingQuality { .. })
        ));
        Ok(())
    }

    fn write_fixture(path: &Path, bytes: &[u8]) -> Result<()> {
        let mut file = File::create(path)?;
        file.write_all(bytes)?;
        Ok(())
    }
}
