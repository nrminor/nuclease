//! Top-level ingress, parsing, and output orchestration.

use needletail::{
    errors::{ParseError, ParseErrorKind},
    parse_fastx_reader,
    parser::{Format, SequenceRecord},
};

use std::{fs::File, marker::PhantomData, path::PathBuf, time::Instant};

use crate::{
    adapter::{AdapterPreset, TrimAdaptersTransform},
    cli::{AdmissionPolicy, Cli, Ingress, UiPolicy},
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
    record::{
        AdmissionEvent, AdmissionReport, InputSource, MateSide, ReadStats, RecordProvenance,
        RecordView,
    },
    report::{self, RunContext as RunSummaryContext, RunLayout},
};

const ADMISSION_WARNING_LIMIT: u64 = 5;

struct SingleEnd;

struct PairedEnd;

trait RecordLayout {
    type Source: RunSource;
    type Readers;
    type Output;
}

impl RecordLayout for SingleEnd {
    type Source = SingleSource;
    type Readers = Box<dyn std::io::Read + Send>;
    type Output = SingleOutputHandle;
}

impl RecordLayout for PairedEnd {
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

struct RunContext<L: RecordLayout> {
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
    admission_policy: AdmissionPolicy,
    progress_every: u64,
    summary: Option<PathBuf>,
    invalid_input_report: Option<PathBuf>,
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
            admission_policy: cli.on_invalid_input,
            progress_every: cli.progress_every,
            summary: cli.summary.clone(),
            invalid_input_report: cli.invalid_input_report.clone(),
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
                config.admission_policy,
                &stats,
                source_error,
            )
        })?;
        let mut progress = ProgressReporter::new(ui.progress_mode, config.progress_every);
        let started_at = Instant::now();
        let admission = RecordAdmission::new(&source, config.admission_policy);

        while let Some(next_record) = parser.next() {
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
        let output_result = output.finish();
        let report_result = stats.finish_admission_report();
        output_result?;
        report_result?;
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
        let admission = RecordAdmission::<PairedEnd>::new(&source, config.admission_policy);

        match readers {
            PairedReaders::Split { left, right } => {
                let mut parser_r1 = parse_fastx_reader(left).map_err(|source_error| {
                    parser_error(
                        source.input_origin("left"),
                        source.input_label(),
                        "left",
                        config.admission_policy,
                        &stats,
                        source_error,
                    )
                })?;
                let mut parser_r2 = parse_fastx_reader(right).map_err(|source_error| {
                    parser_error(
                        source.input_origin("right"),
                        source.input_label(),
                        "right",
                        config.admission_policy,
                        &stats,
                        source_error,
                    )
                })?;

                loop {
                    let next_r1 = parser_r1.next();
                    let next_r2 = parser_r2.next();

                    match (next_r1, next_r2) {
                        (Some(record_r1), Some(record_r2)) => {
                            let parsed_r1 = admission.parse("left", &mut stats, record_r1)?;
                            let left = admission.paired_record(
                                &parsed_r1,
                                MateSide::Left,
                                "left",
                                &stats,
                            )?;
                            stats.record_seen(left.sequence().len());

                            let parsed_r2 = admission.parse("right", &mut stats, record_r2)?;
                            let right = admission.paired_record(
                                &parsed_r2,
                                MateSide::Right,
                                "right",
                                &stats,
                            )?;
                            stats.record_seen(right.sequence().len());

                            let Some(pair) = admission.admit_pair(left, right, &mut stats)? else {
                                continue;
                            };

                            arena.reset();
                            let outcome = plan.execute(pair, &mut arena, &mut stats)?;
                            output.write_outcome(&outcome, &mut stats)?;
                            progress.maybe_report(&stats);
                        }
                        (None, None) => break,
                        (Some(record_r1), None) => {
                            let parsed = admission.parse("left", &mut stats, record_r1)?;
                            admission.missing_mate(&parsed, MateSide::Left, "left", &mut stats)?;
                        }
                        (None, Some(record_r2)) => {
                            let parsed = admission.parse("right", &mut stats, record_r2)?;
                            admission.missing_mate(
                                &parsed,
                                MateSide::Right,
                                "right",
                                &mut stats,
                            )?;
                        }
                    }
                }
            }
            PairedReaders::Interleaved(reader) => {
                let mut parser = parse_fastx_reader(reader).map_err(|source_error| {
                    parser_error(
                        source.input_origin("left"),
                        source.input_label(),
                        "left",
                        config.admission_policy,
                        &stats,
                        source_error,
                    )
                })?;
                let mut left_buffer = InterleavedLeftBuffer::default();

                loop {
                    let next_left = parser.next();
                    let Some(left_record) = next_left else {
                        break;
                    };

                    let parsed_left = admission.parse("left", &mut stats, left_record)?;
                    let left =
                        admission.buffered_left_record(&parsed_left, &stats, &mut left_buffer)?;
                    stats.record_seen(left.sequence().len());

                    let next_right = parser.next();
                    let Some(right_record) = next_right else {
                        admission.missing_mate_record(left, &mut stats)?;
                        break;
                    };
                    let parsed_right = admission.parse("right", &mut stats, right_record)?;
                    let right =
                        admission.paired_record(&parsed_right, MateSide::Right, "right", &stats)?;
                    stats.record_seen(right.sequence().len());

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
        let output_result = output.finish();
        let report_result = stats.finish_admission_report();
        output_result?;
        report_result?;
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

struct RecordAdmission<'source, L: RecordLayout> {
    source: &'source L::Source,
    policy: AdmissionPolicy,
    _layout: PhantomData<L>,
}

impl<'source, L: RecordLayout> RecordAdmission<'source, L> {
    fn new(source: &'source L::Source, policy: AdmissionPolicy) -> Self {
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
            Ok(record) if record.format() == Format::Fastq => Ok(record),
            Ok(record) => Err(unsupported_format_error(
                self.source.input_origin(mate),
                self.source.input_label(),
                mate,
                record.format(),
            )),
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
        let parser_error_kind = parse_error_kind_name(&error.kind);
        let parser_error_message = error.to_string();
        let parser_error_line = (error.position.line > 0).then_some(error.position.line);

        if matches!(
            &error.kind,
            ParseErrorKind::Io | ParseErrorKind::UnknownFormat | ParseErrorKind::EmptyFile
        ) {
            return Err(parser_error(
                self.source.input_origin(mate),
                source,
                mate,
                self.policy,
                stats,
                error,
            ));
        }

        stats.invalid_reads += 1;
        stats.record_admission_event(AdmissionEvent::RecordParseFailure {
            source: source.clone(),
            mate,
            parser_kind: parser_error_kind.to_owned(),
            message: parser_error_message.clone(),
            line: parser_error_line,
            reads_seen: stats.reads_seen,
            pairs_seen: (mate != "single").then_some(stats.pairs_seen),
            continued: false,
        })?;

        if self.policy == AdmissionPolicy::Skip {
            tracing::warn!(
                source,
                mate,
                parser_kind = parser_error_kind,
                parser_error = parser_error_message,
                "record parser error is not recoverable; stopping instead of skipping and continuing"
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

    fn admit_record<'record>(
        &self,
        record: RecordView<'record>,
        stats: &mut ReadStats,
    ) -> Result<Option<RecordView<'record>>> {
        if record.sequence().len() == record.quality().len() {
            return Ok(Some(record));
        }

        self.record_length_mismatch(record, stats)?;

        match self.policy {
            AdmissionPolicy::Error => Err(Self::length_mismatch_error(record)),
            AdmissionPolicy::Skip => {
                Self::warn_length_mismatch(record, stats);
                Ok(None)
            }
        }
    }

    fn record_length_mismatch(&self, record: RecordView<'_>, stats: &mut ReadStats) -> Result<()> {
        stats.invalid_reads += 1;
        stats.record_admission_event(AdmissionEvent::SequenceQualityLengthMismatch {
            source: record.source_display(),
            mate: record.mate_display(),
            header: String::from_utf8_lossy(record.header()).into_owned(),
            sequence_len: record.sequence().len(),
            quality_len: record.quality().len(),
            reads_seen: stats.reads_seen,
            pairs_seen: record
                .provenance()
                .and_then(|provenance| provenance.mate)
                .map(|_| stats.pairs_seen),
            continued: self.policy == AdmissionPolicy::Skip,
        })?;
        Ok(())
    }

    fn length_mismatch_error(record: RecordView<'_>) -> RunError {
        match record.provenance().map(|provenance| provenance.source) {
            Some(InputSource::Ena { accession }) => IndeterminateInputError::Content {
                accession: accession.to_owned(),
                problem: EnaContentProblem::RecordLength {
                    mate: record.mate_display(),
                    header: String::from_utf8_lossy(record.header()).into_owned(),
                    sequence_len: record.sequence().len(),
                    quality_len: record.quality().len(),
                },
            }
            .into(),
            _ => MalformedInputError::RecordLength {
                source_label: record.source_display(),
                mate: record.mate_display(),
                header: String::from_utf8_lossy(record.header()).into_owned(),
                sequence_len: record.sequence().len(),
                quality_len: record.quality().len(),
            }
            .into(),
        }
    }

    fn warn_length_mismatch(record: RecordView<'_>, stats: &mut ReadStats) {
        if stats.should_emit_admission_warning(ADMISSION_WARNING_LIMIT) {
            tracing::warn!(
                source = %record.source_display(),
                mate = %record.mate_display(),
                header = %String::from_utf8_lossy(record.header()),
                sequence_len = record.sequence().len(),
                quality_len = record.quality().len(),
                "skipping input record with mismatched sequence and quality lengths"
            );
        } else if stats.should_emit_admission_suppressed_notice() {
            tracing::warn!("further record-admission warnings suppressed");
        }
    }
}

fn parse_error_kind_name(kind: &ParseErrorKind) -> &'static str {
    match kind {
        ParseErrorKind::Io => "io",
        ParseErrorKind::UnknownFormat => "unknown_format",
        ParseErrorKind::InvalidStart => "invalid_start",
        ParseErrorKind::InvalidSeparator => "invalid_separator",
        ParseErrorKind::UnequalLengths => "unequal_lengths",
        ParseErrorKind::UnexpectedEnd => "unexpected_end",
        ParseErrorKind::EmptyFile => "empty_file",
    }
}

impl RecordAdmission<'_, SingleEnd> {
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

        self.admit_record(record, stats)
    }
}

impl<'source> RecordAdmission<'source, PairedEnd> {
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

    fn missing_mate(
        &self,
        parsed_record: &SequenceRecord<'_>,
        present_mate: MateSide,
        mate_label: &'static str,
        stats: &mut ReadStats,
    ) -> Result<()> {
        let record = self.paired_record(parsed_record, present_mate, mate_label, stats)?;
        stats.record_seen(record.sequence().len());
        self.missing_mate_record(record, stats)
    }

    fn missing_mate_record(&self, record: RecordView<'_>, stats: &mut ReadStats) -> Result<()> {
        let source = record.source_display();
        let mate = record
            .provenance()
            .and_then(|provenance| provenance.mate)
            .ok_or_else(|| InternalError::PlanInvariant {
                detail: "paired record did not carry mate provenance".to_owned(),
            })?;
        let mate_label = record.mate_display();
        let header = String::from_utf8_lossy(record.header()).into_owned();
        let continued = self.policy == AdmissionPolicy::Skip;

        stats.record_admission_event(AdmissionEvent::MissingMate {
            source: source.clone(),
            present_mate: mate,
            header: header.clone(),
            reads_seen: stats.reads_seen,
            pairs_seen: stats.pairs_seen,
            continued,
        })?;

        match self.policy {
            AdmissionPolicy::Error => {
                Err(record_count_error(self.source, mate_label, header, stats))
            }
            AdmissionPolicy::Skip => {
                if stats.should_emit_admission_warning(ADMISSION_WARNING_LIMIT) {
                    tracing::warn!(
                        source,
                        present_mate = mate_label,
                        header,
                        "skipping input record without a mate"
                    );
                } else if stats.should_emit_admission_suppressed_notice() {
                    tracing::warn!("further record-admission warnings suppressed");
                }
                Ok(())
            }
        }
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
        stats.pairs_seen += 1;

        let left_length_mismatch = left.sequence().len() != left.quality().len();
        let right_length_mismatch = right.sequence().len() != right.quality().len();
        if left_length_mismatch || right_length_mismatch {
            stats.invalid_pairs += 1;
            if left_length_mismatch {
                self.record_length_mismatch(left, stats)?;
            }
            if right_length_mismatch {
                self.record_length_mismatch(right, stats)?;
            }

            return match self.policy {
                AdmissionPolicy::Error => {
                    Err(Self::length_mismatch_error(if left_length_mismatch {
                        left
                    } else {
                        right
                    }))
                }
                AdmissionPolicy::Skip => {
                    if left_length_mismatch {
                        Self::warn_length_mismatch(left, stats);
                    }
                    if right_length_mismatch {
                        Self::warn_length_mismatch(right, stats);
                    }
                    Ok(None)
                }
            };
        }

        if left.pair_key() == right.pair_key() {
            return Ok(Some(RecordPair { left, right }));
        }

        stats.invalid_pairs += 1;
        stats.record_admission_event(AdmissionEvent::PairIdentifierMismatch {
            source: left.source_display(),
            left_header: String::from_utf8_lossy(left.header()).into_owned(),
            right_header: String::from_utf8_lossy(right.header()).into_owned(),
            reads_seen: stats.reads_seen,
            pairs_seen: stats.pairs_seen,
            continued: self.policy == AdmissionPolicy::Skip,
        })?;

        match self.policy {
            AdmissionPolicy::Error => Err(Self::mate_identifier_error(left, right)),
            AdmissionPolicy::Skip => {
                Self::warn_identifier_mismatch(left, right, stats);
                Ok(None)
            }
        }
    }

    fn mate_identifier_error(left: RecordView<'_>, right: RecordView<'_>) -> RunError {
        match left.provenance().map(|provenance| provenance.source) {
            Some(InputSource::Ena { accession }) => IndeterminateInputError::Content {
                accession: accession.to_owned(),
                problem: EnaContentProblem::MateIdentifier {
                    left_header: String::from_utf8_lossy(left.header()).into_owned(),
                    right_header: String::from_utf8_lossy(right.header()).into_owned(),
                },
            }
            .into(),
            _ => MalformedInputError::MateIdentifier {
                source_label: left.source_display(),
                left_mate: left.mate_display(),
                right_mate: right.mate_display(),
                left_header: String::from_utf8_lossy(left.header()).into_owned(),
                right_header: String::from_utf8_lossy(right.header()).into_owned(),
            }
            .into(),
        }
    }

    fn warn_identifier_mismatch(
        left: RecordView<'_>,
        right: RecordView<'_>,
        stats: &mut ReadStats,
    ) {
        if stats.should_emit_admission_warning(ADMISSION_WARNING_LIMIT) {
            tracing::warn!(
                source = %left.source_display(),
                left_mate = %left.mate_display(),
                right_mate = %right.mate_display(),
                left_header = %String::from_utf8_lossy(left.header()),
                right_header = %String::from_utf8_lossy(right.header()),
                "skipping paired input with mismatched identifiers"
            );
        } else if stats.should_emit_admission_suppressed_notice() {
            tracing::warn!("further record-admission warnings suppressed");
        }
    }
}

fn missing_quality_error<S: RunSource>(
    source: &S,
    mate: &'static str,
    stats: &ReadStats,
) -> RunError {
    match source.input_origin(mate) {
        InputOrigin::Ena(accession) => IndeterminateInputError::Content {
            accession: accession.to_string(),
            problem: EnaContentProblem::MissingQuality {
                mate,
                reads_seen: stats.reads_seen,
                pairs_seen: stats.pairs_seen,
                invalid_reads: stats.invalid_reads,
                invalid_pairs: stats.invalid_pairs,
            },
        }
        .into(),
        InputOrigin::Local(_) => MalformedInputError::MissingQuality {
            source_label: source.input_label(),
            mate,
            reads_seen: stats.reads_seen,
            pairs_seen: stats.pairs_seen,
            invalid_reads: stats.invalid_reads,
            invalid_pairs: stats.invalid_pairs,
        }
        .into(),
    }
}

fn read_stats(config: &RunConfig) -> Result<ReadStats> {
    let mut stats = ReadStats::default();
    if let Some(path) = &config.invalid_input_report {
        stats.set_admission_report(AdmissionReport::create(path)?);
    }
    Ok(stats)
}

fn parser_error(
    origin: InputOrigin<'_>,
    source_label: String,
    mate: &'static str,
    policy: AdmissionPolicy,
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
        InputOrigin::Local(_)
            if matches!(
                &source.kind,
                ParseErrorKind::UnknownFormat | ParseErrorKind::EmptyFile
            ) =>
        {
            MalformedInputError::UnreadableFastq {
                source_label,
                mate,
                reads_seen: stats.reads_seen,
                pairs_seen: stats.pairs_seen,
                parser_error_kind: parse_error_kind_name(&source.kind),
                source,
            }
            .into()
        }
        InputOrigin::Local(_) => MalformedInputError::LocalParser {
            source_label,
            mate,
            policy: policy.to_string(),
            reads_seen: stats.reads_seen,
            pairs_seen: stats.pairs_seen,
            invalid_reads: stats.invalid_reads,
            invalid_pairs: stats.invalid_pairs,
            parser_error_kind: parse_error_kind_name(&source.kind).to_owned(),
            source,
        }
        .into(),
    }
}

fn unsupported_format_error(
    origin: InputOrigin<'_>,
    source_label: String,
    mate: &'static str,
    format: Format,
) -> RunError {
    match origin {
        InputOrigin::Ena(accession) => IndeterminateInputError::Content {
            accession: accession.to_string(),
            problem: EnaContentProblem::UnsupportedFormat { mate, format },
        }
        .into(),
        InputOrigin::Local(_) => MalformedInputError::UnsupportedFormat {
            source_label,
            mate,
            format,
        }
        .into(),
    }
}

fn record_count_error(
    source: &PairedSource,
    present_mate: &'static str,
    header: String,
    stats: &ReadStats,
) -> RunError {
    match source {
        PairedSource::Ena { accession } => IndeterminateInputError::Content {
            accession: accession.to_string(),
            problem: EnaContentProblem::RecordCount {
                complete_pairs_seen: stats.pairs_seen,
                present_mate,
                header,
            },
        }
        .into(),
        PairedSource::LocalInterleaved { .. } => MalformedInputError::InterleavedRecordCount {
            source_label: source.input_label(),
            present_mate,
            header,
            complete_pairs_seen: stats.pairs_seen,
            reads_seen: stats.reads_seen,
        }
        .into(),
        PairedSource::LocalSplit { .. } => MalformedInputError::PairedRecordCount {
            source_label: source.input_label(),
            present_mate,
            header,
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
        path::{Path, PathBuf},
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
        cli::{AdmissionPolicy, Cli},
        ena::Accession,
        error::{
            EnaContentProblem, IndeterminateInputError, IoError, MalformedInputError, RunError,
        },
        output::{
            InterleavedOutput, OutputArgs, OutputEncoding, OutputFormat, PairedRecordOutput,
            SingleOutput, SingleRecordOutput, StreamSink,
        },
        record::{AdmissionEvent, AdmissionReport, MateSide, ReadStats, RecordView},
    };

    use super::{
        InputOrigin, PairedEnd, PairedEndContext, PairedSource, RecordAdmission, RunConfig,
        SingleEnd, SingleSource, missing_quality_error, parser_error, record_count_error,
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
            on_invalid_input: AdmissionPolicy::Error,
            out: None,
            out1: None,
            out2: None,
            progress_every: 100_000,
            summary: None,
            invalid_input_report: None,
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
                reads_seen: 3,
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
                AdmissionPolicy::Error,
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
            AdmissionPolicy::Error,
            &stats,
            ParseError::from(io::Error::other("local read failed")),
        );
        let RunError::Io(IoError::LocalFastqRead { source, .. }) = error else {
            panic!("local parser I/O should retain the I/O category");
        };
        assert_eq!(source.kind, ParseErrorKind::Io);

        for source in parser_errors.into_iter().filter(|source| {
            !matches!(
                &source.kind,
                ParseErrorKind::Io | ParseErrorKind::UnknownFormat | ParseErrorKind::EmptyFile
            )
        }) {
            let expected_kind = source.kind.clone();
            let error = parser_error(
                InputOrigin::Local(local_path),
                "local:reads.fastq".to_owned(),
                "single",
                AdmissionPolicy::Error,
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
    fn local_parser_source_classification_preserves_typed_source() {
        let stats = ReadStats {
            reads_seen: 7,
            pairs_seen: 3,
            ..ReadStats::default()
        };
        let local_path = Path::new("reads.fastq");

        for source in [
            ParseError::new_unknown_format(b'!'),
            ParseError::new_empty_file(),
        ] {
            let expected_kind = source.kind.clone();
            let error = parser_error(
                InputOrigin::Local(local_path),
                "local:reads.fastq".to_owned(),
                "single",
                AdmissionPolicy::Error,
                &stats,
                source,
            );
            assert!(
                error.source().is_some(),
                "source-classification failures should retain the parser source"
            );
            assert!(matches!(
                error,
                RunError::MalformedInput(error)
                    if matches!(
                        &*error,
                        MalformedInputError::UnreadableFastq { source, .. }
                            if source.kind == expected_kind
                    )
            ));
        }
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
            record_count_error(
                &PairedSource::Ena { accession },
                "left",
                "read2/1".to_owned(),
                &stats,
            ),
            RunError::IndeterminateInput(IndeterminateInputError::Content {
                problem: EnaContentProblem::RecordCount {
                    complete_pairs_seen: 1,
                    present_mate: "left",
                    ref header,
                },
                ..
            }) if header == "read2/1"
        ));
        assert!(matches!(
            record_count_error(
                &PairedSource::LocalInterleaved {
                    input: "reads.fastq".into(),
                },
                "left",
                "read2/1".to_owned(),
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
                "right",
                "read2/2".to_owned(),
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
                problem: EnaContentProblem::MissingQuality { mate: "single", .. },
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

    #[test]
    fn record_length_errors_preserve_input_origin_at_admission() -> Result<()> {
        let local_source = SingleSource::Local {
            input: PathBuf::from("reads.fastq"),
        };
        let local_admission =
            RecordAdmission::<SingleEnd>::new(&local_source, AdmissionPolicy::Error);
        let local_record =
            RecordView::new(b"read1", b"ACGT", b"I").with_provenance(local_source.provenance());
        let mut local_stats = ReadStats::default();
        local_stats.record_seen(local_record.sequence().len());

        let Err(local_error) = local_admission.admit_record(local_record, &mut local_stats) else {
            panic!("local length mismatch should fail admission");
        };
        assert!(matches!(
            local_error,
            RunError::MalformedInput(error)
                if matches!(*error, MalformedInputError::RecordLength { .. })
        ));
        assert_eq!(local_stats.invalid_reads, 1);
        assert_eq!(
            local_stats
                .admission_event_counts
                .get("sequence_quality_length_mismatch"),
            Some(&1)
        );

        let ena_source = SingleSource::Ena {
            accession: Accession::new("SRR35939766")?,
        };
        let ena_admission = RecordAdmission::<SingleEnd>::new(&ena_source, AdmissionPolicy::Error);
        let ena_record =
            RecordView::new(b"read1", b"ACGT", b"I").with_provenance(ena_source.provenance());
        let mut ena_stats = ReadStats::default();
        ena_stats.record_seen(ena_record.sequence().len());

        let Err(ena_error) = ena_admission.admit_record(ena_record, &mut ena_stats) else {
            panic!("ENA length mismatch should fail admission");
        };
        assert!(matches!(
            ena_error,
            RunError::IndeterminateInput(IndeterminateInputError::Content {
                problem: EnaContentProblem::RecordLength { .. },
                ..
            })
        ));
        assert_eq!(ena_stats.invalid_reads, 1);
        Ok(())
    }

    #[test]
    fn mate_identifier_errors_preserve_input_origin_at_admission() -> Result<()> {
        let local_source = PairedSource::LocalInterleaved {
            input: PathBuf::from("reads.fastq"),
        };
        let local_admission =
            RecordAdmission::<PairedEnd>::new(&local_source, AdmissionPolicy::Error);
        let local_left = RecordView::new(b"read1/1", b"A", b"I")
            .with_provenance(local_source.provenance(MateSide::Left));
        let local_right = RecordView::new(b"read2/2", b"T", b"I")
            .with_provenance(local_source.provenance(MateSide::Right));
        let mut local_stats = ReadStats::default();
        local_stats.record_seen(1);
        local_stats.record_seen(1);

        let Err(local_error) =
            local_admission.admit_pair(local_left, local_right, &mut local_stats)
        else {
            panic!("local mate mismatch should fail admission");
        };
        assert!(matches!(
            local_error,
            RunError::MalformedInput(error)
                if matches!(*error, MalformedInputError::MateIdentifier { .. })
        ));
        assert_eq!(local_stats.invalid_pairs, 1);

        let ena_source = PairedSource::Ena {
            accession: Accession::new("SRR35939766")?,
        };
        let ena_admission = RecordAdmission::<PairedEnd>::new(&ena_source, AdmissionPolicy::Error);
        let ena_left = RecordView::new(b"read1/1", b"A", b"I")
            .with_provenance(ena_source.provenance(MateSide::Left));
        let ena_right = RecordView::new(b"read2/2", b"T", b"I")
            .with_provenance(ena_source.provenance(MateSide::Right));
        let mut ena_stats = ReadStats::default();
        ena_stats.record_seen(1);
        ena_stats.record_seen(1);

        let Err(ena_error) = ena_admission.admit_pair(ena_left, ena_right, &mut ena_stats) else {
            panic!("ENA mate mismatch should fail admission");
        };
        assert!(matches!(
            ena_error,
            RunError::IndeterminateInput(IndeterminateInputError::Content {
                problem: EnaContentProblem::MateIdentifier { .. },
                ..
            })
        ));
        assert_eq!(ena_stats.invalid_pairs, 1);
        Ok(())
    }

    #[test]
    fn parser_io_failure_is_not_an_admission_event() -> Result<()> {
        let temp = tempdir()?;
        let report_path = temp.path().join("invalid-input.jsonl");
        let source = SingleSource::Local {
            input: PathBuf::from("reads.fastq"),
        };
        let admission = RecordAdmission::<SingleEnd>::new(&source, AdmissionPolicy::Skip);
        let mut stats = ReadStats::default();
        stats.set_admission_report(AdmissionReport::create(&report_path)?);

        let error = admission
            .parse(
                "single",
                &mut stats,
                Err(ParseError::from(io::Error::other("synthetic read failure"))),
            )
            .expect_err("I/O failure should remain fatal");

        assert!(
            matches!(error, RunError::Io(IoError::LocalFastqRead { .. })),
            "I/O parser failure should be classified as local read I/O"
        );
        assert_eq!(stats.invalid_reads, 0);
        assert!(stats.admission_event_counts.is_empty());
        assert_eq!(std::fs::read_to_string(report_path)?, "");
        Ok(())
    }

    #[test]
    fn parser_source_classification_failures_are_not_admission_events() -> Result<()> {
        let temp = tempdir()?;
        let report_path = temp.path().join("invalid-input.jsonl");
        let source = SingleSource::Local {
            input: PathBuf::from("reads.fastq"),
        };
        let admission = RecordAdmission::<SingleEnd>::new(&source, AdmissionPolicy::Skip);
        let mut stats = ReadStats::default();
        stats.set_admission_report(AdmissionReport::create(&report_path)?);

        for parse_error in [
            ParseError::new_unknown_format(b'!'),
            ParseError::new_empty_file(),
        ] {
            let error = admission
                .parse("single", &mut stats, Err(parse_error))
                .expect_err("source classification failure should remain fatal");
            assert!(
                error
                    .to_string()
                    .contains("input source did not provide a readable FASTQ stream")
            );
        }

        assert_eq!(stats.invalid_reads, 0);
        assert!(stats.admission_event_counts.is_empty());
        assert_eq!(std::fs::read_to_string(report_path)?, "");
        Ok(())
    }

    #[test]
    fn malformed_pair_admission_is_atomic_for_each_invalid_mate() -> Result<()> {
        for (invalid_mate, left_quality, right_quality) in [
            ("left", b"I".as_slice(), b"KKKK".as_slice()),
            ("right", b"IIII".as_slice(), b"K".as_slice()),
        ] {
            for policy in [AdmissionPolicy::Error, AdmissionPolicy::Skip] {
                let source = PairedSource::LocalSplit {
                    input1: PathBuf::from("reads_1.fastq"),
                    input2: PathBuf::from("reads_2.fastq"),
                };
                let admission = RecordAdmission::<PairedEnd>::new(&source, policy);
                let left = RecordView::new(b"read1/1", b"AAAA", left_quality)
                    .with_provenance(source.provenance(MateSide::Left));
                let right = RecordView::new(b"read1/2", b"TTTT", right_quality)
                    .with_provenance(source.provenance(MateSide::Right));
                let mut stats = ReadStats::default();
                stats.record_seen(left.sequence().len());
                stats.record_seen(right.sequence().len());

                let result = admission.admit_pair(left, right, &mut stats);
                match policy {
                    AdmissionPolicy::Error => {
                        assert!(
                            result.is_err(),
                            "Error should halt on an invalid {invalid_mate} mate"
                        );
                    }
                    AdmissionPolicy::Skip => {
                        assert!(
                            result?.is_none(),
                            "Skip should discard the pair when its {invalid_mate} mate is invalid"
                        );
                    }
                }

                assert_eq!(stats.pairs_seen, 1);
                assert_eq!(stats.invalid_pairs, 1);
                assert_eq!(stats.invalid_reads, 1);
                assert_eq!(
                    stats
                        .admission_event_counts
                        .get("sequence_quality_length_mismatch"),
                    Some(&1)
                );
                let [
                    AdmissionEvent::SequenceQualityLengthMismatch {
                        mate, continued, ..
                    },
                ] = stats.admission_samples.as_slice()
                else {
                    panic!("exactly one length-mismatch event should identify the invalid mate");
                };
                assert_eq!(*mate, invalid_mate);
                assert_eq!(*continued, policy == AdmissionPolicy::Skip);
            }
        }
        Ok(())
    }

    #[test]
    fn malformed_pair_classification_does_not_depend_on_policy() -> Result<()> {
        for policy in [AdmissionPolicy::Error, AdmissionPolicy::Skip] {
            let source = PairedSource::LocalSplit {
                input1: PathBuf::from("reads_1.fastq"),
                input2: PathBuf::from("reads_2.fastq"),
            };
            let admission = RecordAdmission::<PairedEnd>::new(&source, policy);
            let left = RecordView::new(b"read1/1", b"AAAA", b"I")
                .with_provenance(source.provenance(MateSide::Left));
            let right = RecordView::new(b"read1/2", b"TTTT", b"K")
                .with_provenance(source.provenance(MateSide::Right));
            let mut stats = ReadStats::default();
            stats.record_seen(left.sequence().len());
            stats.record_seen(right.sequence().len());

            let result = admission.admit_pair(left, right, &mut stats);
            match policy {
                AdmissionPolicy::Error => {
                    assert!(
                        result.is_err(),
                        "Error policy should halt on the malformed pair"
                    );
                }
                AdmissionPolicy::Skip => {
                    assert!(result?.is_none(), "Skip should discard the malformed pair");
                }
            }

            assert_eq!(stats.pairs_seen, 1);
            assert_eq!(stats.invalid_pairs, 1);
            assert_eq!(stats.invalid_reads, 2);
            assert_eq!(
                stats
                    .admission_event_counts
                    .get("sequence_quality_length_mismatch"),
                Some(&2)
            );
        }
        Ok(())
    }

    fn write_fixture(path: &Path, bytes: &[u8]) -> Result<()> {
        let mut file = File::create(path)?;
        file.write_all(bytes)?;
        Ok(())
    }
}
