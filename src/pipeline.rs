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
        self, InputOrigin, PairConstructionError, PairedInputLayout, Result, UnavailableInputError,
        UsageError,
    },
    filter::{MaxNsFilter, MinEntropyFilter, MinLengthFilter, MinMeanQualityFilter},
    observer::{InvalidInputEvent, InvalidInputReport, RunObserver},
    output::{OutputArgs, PairedOutputHandle, SingleOutputHandle, UnitOutput},
    pair_merge::MergePairsTransform,
    plan::{BuildPlan, Execute, Execution, Logical, OrphanPolicy, Plan, TransformArena},
    progress::ProgressReporter,
    quality::{QualityBinCount, QualityBinTransform, QualityTrimTransform},
    record::{InputSource, LeftMate, MateRecord, MateSide, RecordPair, RecordView, RightMate},
    report::{self, RunContext as RunSummaryContext, RunLayout},
};

struct SingleEnd;

struct PairedEnd;

trait RecordLayout {
    type Source;
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

impl SingleSource {
    fn input_label(&self) -> String {
        match self {
            Self::Ena { accession } => format!("ena:{accession}"),
            Self::Local { input } => format!("local:{}", input.display()),
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

    fn input_origin(&self) -> InputOrigin<'_> {
        match self {
            Self::Ena { accession } => InputOrigin::Ena(accession),
            Self::Local { input } => InputOrigin::Local(input),
        }
    }

    fn input_source(&self) -> InputSource<'_> {
        match self {
            Self::Ena { accession } => InputSource::Ena {
                accession: accession.as_str(),
            },
            Self::Local { input } => InputSource::LocalSingle { input },
        }
    }
}

impl PairedSource {
    fn input_label(&self) -> String {
        match self {
            Self::Ena { accession } => format!("ena:{accession}"),
            Self::LocalInterleaved { input } => format!("local-interleaved:{}", input.display()),
            Self::LocalSplit { input1, input2 } => {
                format!("local-paired:{}|{}", input1.display(), input2.display())
            }
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

    fn left_input_origin(&self) -> InputOrigin<'_> {
        match self {
            Self::Ena { accession } => InputOrigin::Ena(accession),
            Self::LocalInterleaved { input } => InputOrigin::Local(input),
            Self::LocalSplit { input1, .. } => InputOrigin::Local(input1),
        }
    }

    fn right_input_origin(&self) -> InputOrigin<'_> {
        match self {
            Self::Ena { accession } => InputOrigin::Ena(accession),
            Self::LocalInterleaved { input } => InputOrigin::Local(input),
            Self::LocalSplit { input2, .. } => InputOrigin::Local(input2),
        }
    }

    const fn input_layout(&self) -> PairedInputLayout {
        match self {
            Self::LocalInterleaved { .. } => PairedInputLayout::Interleaved,
            Self::Ena { .. } | Self::LocalSplit { .. } => PairedInputLayout::Split,
        }
    }

    fn input_source(&self) -> InputSource<'_> {
        match self {
            Self::Ena { accession } => InputSource::Ena {
                accession: accession.as_str(),
            },
            Self::LocalInterleaved { input } => InputSource::LocalInterleavedPaired { input },
            Self::LocalSplit { input1, input2 } => InputSource::LocalPaired { input1, input2 },
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
        let mut observer = run_observer(config, source.input_label())?;
        // build out the mutable state needed to run the application loop
        let mut parser = parse_fastx_reader(reader).map_err(|source_error| {
            error::constructors::single_parser_error(
                source.input_origin(),
                source.input_label(),
                &config.admission_policy,
                &observer,
                source_error,
            )
        })?;
        let mut progress = ProgressReporter::new(ui.progress_mode, config.progress_every);
        let started_at = Instant::now();
        let admission = RecordAdmission::new(&source, config.admission_policy);

        while let Some(next_record) = parser.next() {
            let Some(parsed_record) = (match next_record {
                Ok(record) if record.format() == Format::Fastq => Some(record),
                Ok(record) => {
                    return Err(error::constructors::single_unsupported_format_error(
                        source.input_origin(),
                        source.input_label(),
                        record.format(),
                    ));
                }
                Err(error)
                    if config.admission_policy == AdmissionPolicy::Skip
                        && error.kind == ParseErrorKind::UnequalLengths =>
                {
                    observer.record_unparsed_seen();
                    observer.recoverable_single_parser_failure(&error)?;
                    None
                }
                Err(error) => {
                    return single_parser_failure(
                        &source,
                        config.admission_policy,
                        &mut observer,
                        error,
                    );
                }
            }) else {
                continue;
            };
            let Some(record) = admission.single(&parsed_record, &mut observer)? else {
                continue;
            };
            arena.reset();

            let outcome = plan.execute(record, &mut arena, &mut observer)?;
            output.write_outcome(&outcome, &mut observer)?;

            progress.maybe_report(&observer);
        }

        progress.finish();
        let output_result = output.finish();
        let report_result = observer.finish_invalid_input_report();
        output_result?;
        report_result?;
        let summary = report::RunSummary::from_observer(
            source.summary_context(),
            &observer,
            started_at.elapsed(),
        );
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
        let mut observer = run_observer(config, source.input_label())?;
        let mut progress = ProgressReporter::new(ui.progress_mode, config.progress_every);
        let started_at = Instant::now();
        let admission = RecordAdmission::<PairedEnd>::new(&source, config.admission_policy);

        match readers {
            PairedReaders::Split { left, right } => {
                let mut parser_r1 = parse_fastx_reader(left).map_err(|source_error| {
                    error::constructors::left_parser_error(
                        source.left_input_origin(),
                        source.input_label(),
                        &config.admission_policy,
                        &observer,
                        source_error,
                    )
                })?;
                let mut parser_r2 = parse_fastx_reader(right).map_err(|source_error| {
                    error::constructors::right_parser_error(
                        source.right_input_origin(),
                        source.input_label(),
                        &config.admission_policy,
                        &observer,
                        source_error,
                    )
                })?;

                loop {
                    let next_r1 = parser_r1.next();
                    let next_r2 = parser_r2.next();

                    match (next_r1, next_r2) {
                        (Some(record_r1), Some(record_r2)) => {
                            observer.pairs_seen += 1;
                            let (parsed_r1, left_failure) = match record_r1 {
                                Ok(record) if record.format() == Format::Fastq => {
                                    observer.record_seen(record.raw_seq().len());
                                    (Some(record), None)
                                }
                                Ok(record) => {
                                    return Err(
                                        error::constructors::left_unsupported_format_error(
                                            source.left_input_origin(),
                                            source.input_label(),
                                            record.format(),
                                        ),
                                    );
                                }
                                Err(error)
                                    if config.admission_policy == AdmissionPolicy::Skip
                                        && error.kind == ParseErrorKind::UnequalLengths =>
                                {
                                    observer.record_unparsed_seen();
                                    (None, Some(error))
                                }
                                Err(error) => {
                                    return left_parser_failure(
                                        &source,
                                        config.admission_policy,
                                        &mut observer,
                                        error,
                                    );
                                }
                            };
                            let (parsed_r2, right_failure) = match record_r2 {
                                Ok(record) if record.format() == Format::Fastq => {
                                    observer.record_seen(record.raw_seq().len());
                                    (Some(record), None)
                                }
                                Ok(record) => {
                                    return Err(
                                        error::constructors::right_unsupported_format_error(
                                            source.right_input_origin(),
                                            source.input_label(),
                                            record.format(),
                                        ),
                                    );
                                }
                                Err(error)
                                    if config.admission_policy == AdmissionPolicy::Skip
                                        && error.kind == ParseErrorKind::UnequalLengths =>
                                {
                                    observer.record_unparsed_seen();
                                    (None, Some(error))
                                }
                                Err(error) => {
                                    let right_has_record_slot =
                                        parser_failure_has_record_slot(&error);
                                    if left_failure.is_some() || right_has_record_slot {
                                        observer.invalid_pairs += 1;
                                    }
                                    if right_has_record_slot {
                                        observer.record_unparsed_seen();
                                    }
                                    if let Some(left_failure) = &left_failure {
                                        observer.recoverable_paired_parser_failure(
                                            MateSide::Left,
                                            left_failure,
                                            false,
                                        )?;
                                    }
                                    return finish_right_parser_failure(
                                        &source,
                                        config.admission_policy,
                                        &mut observer,
                                        error,
                                    );
                                }
                            };

                            if let Some(left_failure) = &left_failure {
                                observer.recoverable_paired_parser_failure(
                                    MateSide::Left,
                                    left_failure,
                                    true,
                                )?;
                            }
                            if let Some(right_failure) = &right_failure {
                                observer.recoverable_paired_parser_failure(
                                    MateSide::Right,
                                    right_failure,
                                    true,
                                )?;
                            }

                            let left = match parsed_r1.as_ref() {
                                Some(parsed) => {
                                    let record = admission.left_record(parsed, &observer)?;
                                    Some(record)
                                }
                                None => None,
                            };
                            let right = match parsed_r2.as_ref() {
                                Some(parsed) => {
                                    let record = admission.right_record(parsed, &observer)?;
                                    Some(record)
                                }
                                None => None,
                            };

                            let (Some(left), Some(right)) = (left, right) else {
                                observer.invalid_pairs += 1;
                                continue;
                            };

                            let Some(pair) = admission.admit_pair(left, right, &mut observer)?
                            else {
                                continue;
                            };

                            arena.reset();
                            let outcome = plan.execute(pair, &mut arena, &mut observer)?;
                            output.write_outcome(&outcome, &mut observer)?;
                            progress.maybe_report(&observer);
                        }
                        (None, None) => break,
                        (Some(record_r1), None) => {
                            let parsed = match record_r1 {
                                Ok(record) if record.format() == Format::Fastq => Some(record),
                                Ok(record) => {
                                    return Err(
                                        error::constructors::left_unsupported_format_error(
                                            source.left_input_origin(),
                                            source.input_label(),
                                            record.format(),
                                        ),
                                    );
                                }
                                Err(error)
                                    if config.admission_policy == AdmissionPolicy::Skip
                                        && error.kind == ParseErrorKind::UnequalLengths =>
                                {
                                    observer.record_unparsed_seen();
                                    observer.recoverable_paired_parser_failure(
                                        MateSide::Left,
                                        &error,
                                        true,
                                    )?;
                                    None
                                }
                                Err(error) => {
                                    return left_parser_failure(
                                        &source,
                                        config.admission_policy,
                                        &mut observer,
                                        error,
                                    );
                                }
                            };
                            if let Some(parsed) = parsed.as_ref() {
                                admission.missing_right(parsed, &mut observer)?;
                            }
                        }
                        (None, Some(record_r2)) => {
                            let parsed = match record_r2 {
                                Ok(record) if record.format() == Format::Fastq => Some(record),
                                Ok(record) => {
                                    return Err(
                                        error::constructors::right_unsupported_format_error(
                                            source.right_input_origin(),
                                            source.input_label(),
                                            record.format(),
                                        ),
                                    );
                                }
                                Err(error)
                                    if config.admission_policy == AdmissionPolicy::Skip
                                        && error.kind == ParseErrorKind::UnequalLengths =>
                                {
                                    observer.record_unparsed_seen();
                                    observer.recoverable_paired_parser_failure(
                                        MateSide::Right,
                                        &error,
                                        true,
                                    )?;
                                    None
                                }
                                Err(error) => {
                                    return right_parser_failure(
                                        &source,
                                        config.admission_policy,
                                        &mut observer,
                                        error,
                                    );
                                }
                            };
                            if let Some(parsed) = parsed.as_ref() {
                                admission.missing_left(parsed, &mut observer)?;
                            }
                        }
                    }
                }
            }
            PairedReaders::Interleaved(reader) => {
                let mut parser = parse_fastx_reader(reader).map_err(|source_error| {
                    error::constructors::left_parser_error(
                        source.left_input_origin(),
                        source.input_label(),
                        &config.admission_policy,
                        &observer,
                        source_error,
                    )
                })?;
                let mut left_buffer = InterleavedLeftBuffer::default();

                loop {
                    let next_left = parser.next();
                    let Some(left_record) = next_left else {
                        break;
                    };

                    let (parsed_left, left_failure) = match left_record {
                        Ok(record) if record.format() == Format::Fastq => (Some(record), None),
                        Ok(record) => {
                            return Err(error::constructors::left_unsupported_format_error(
                                source.left_input_origin(),
                                source.input_label(),
                                record.format(),
                            ));
                        }
                        Err(error)
                            if config.admission_policy == AdmissionPolicy::Skip
                                && error.kind == ParseErrorKind::UnequalLengths =>
                        {
                            observer.record_unparsed_seen();
                            (None, Some(error))
                        }
                        Err(error) => {
                            return left_parser_failure(
                                &source,
                                config.admission_policy,
                                &mut observer,
                                error,
                            );
                        }
                    };
                    let left = match parsed_left {
                        Some(parsed) => {
                            let record = admission.buffered_left_record(
                                &parsed,
                                &observer,
                                &mut left_buffer,
                            )?;
                            observer.record_seen(record.into_record().sequence().len());
                            Some(record)
                        }
                        None => None,
                    };

                    let next_right = parser.next();
                    let Some(right_record) = next_right else {
                        if let Some(left_failure) = &left_failure {
                            observer.recoverable_paired_parser_failure(
                                MateSide::Left,
                                left_failure,
                                true,
                            )?;
                        }
                        if let Some(left) = left {
                            admission.missing_right_record(left, &mut observer)?;
                        }
                        break;
                    };
                    observer.pairs_seen += 1;
                    let (parsed_right, right_failure) = match right_record {
                        Ok(record) if record.format() == Format::Fastq => {
                            observer.record_seen(record.raw_seq().len());
                            (Some(record), None)
                        }
                        Ok(record) => {
                            return Err(error::constructors::right_unsupported_format_error(
                                source.right_input_origin(),
                                source.input_label(),
                                record.format(),
                            ));
                        }
                        Err(error)
                            if config.admission_policy == AdmissionPolicy::Skip
                                && error.kind == ParseErrorKind::UnequalLengths =>
                        {
                            observer.record_unparsed_seen();
                            (None, Some(error))
                        }
                        Err(error) => {
                            let right_has_record_slot = parser_failure_has_record_slot(&error);
                            if left_failure.is_some() || right_has_record_slot {
                                observer.invalid_pairs += 1;
                            }
                            if right_has_record_slot {
                                observer.record_unparsed_seen();
                            }
                            if let Some(left_failure) = &left_failure {
                                observer.recoverable_paired_parser_failure(
                                    MateSide::Left,
                                    left_failure,
                                    false,
                                )?;
                            }
                            return finish_right_parser_failure(
                                &source,
                                config.admission_policy,
                                &mut observer,
                                error,
                            );
                        }
                    };
                    if let Some(left_failure) = &left_failure {
                        observer.recoverable_paired_parser_failure(
                            MateSide::Left,
                            left_failure,
                            true,
                        )?;
                    }
                    if let Some(right_failure) = &right_failure {
                        observer.recoverable_paired_parser_failure(
                            MateSide::Right,
                            right_failure,
                            true,
                        )?;
                    }
                    let right = match parsed_right.as_ref() {
                        Some(parsed) => {
                            let record = admission.right_record(parsed, &observer)?;
                            Some(record)
                        }
                        None => None,
                    };

                    let (Some(left), Some(right)) = (left, right) else {
                        observer.invalid_pairs += 1;
                        continue;
                    };

                    let Some(pair) = admission.admit_pair(left, right, &mut observer)? else {
                        continue;
                    };

                    arena.reset();
                    let outcome = plan.execute(pair, &mut arena, &mut observer)?;
                    output.write_outcome(&outcome, &mut observer)?;
                    progress.maybe_report(&observer);
                }
            }
        }

        progress.finish();
        let output_result = output.finish();
        let report_result = observer.finish_invalid_input_report();
        output_result?;
        report_result?;
        let summary = report::RunSummary::from_observer(
            source.summary_context(),
            &observer,
            started_at.elapsed(),
        );
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

fn single_parser_failure<T>(
    source: &SingleSource,
    policy: AdmissionPolicy,
    observer: &mut RunObserver,
    error: ParseError,
) -> Result<T> {
    let source_label = source.input_label();
    if matches!(
        &error.kind,
        ParseErrorKind::Io | ParseErrorKind::UnknownFormat | ParseErrorKind::EmptyFile
    ) {
        return Err(error::constructors::single_parser_error(
            source.input_origin(),
            source_label,
            &policy,
            observer,
            error,
        ));
    }
    observer.record_unparsed_seen();
    observer.terminal_single_parser_failure(&error, policy)?;
    Err(error::constructors::single_parser_error(
        source.input_origin(),
        source_label,
        &policy,
        observer,
        error,
    ))
}

fn left_parser_failure<T>(
    source: &PairedSource,
    policy: AdmissionPolicy,
    observer: &mut RunObserver,
    error: ParseError,
) -> Result<T> {
    let source_label = source.input_label();
    if matches!(
        &error.kind,
        ParseErrorKind::Io | ParseErrorKind::UnknownFormat | ParseErrorKind::EmptyFile
    ) {
        return Err(error::constructors::left_parser_error(
            source.left_input_origin(),
            source_label,
            &policy,
            observer,
            error,
        ));
    }
    observer.record_unparsed_seen();
    observer.terminal_paired_parser_failure(MateSide::Left, &error, policy)?;
    Err(error::constructors::left_parser_error(
        source.left_input_origin(),
        source_label,
        &policy,
        observer,
        error,
    ))
}

fn right_parser_failure<T>(
    source: &PairedSource,
    policy: AdmissionPolicy,
    observer: &mut RunObserver,
    error: ParseError,
) -> Result<T> {
    if parser_failure_has_record_slot(&error) {
        observer.record_unparsed_seen();
    }
    finish_right_parser_failure(source, policy, observer, error)
}

fn finish_right_parser_failure<T>(
    source: &PairedSource,
    policy: AdmissionPolicy,
    observer: &mut RunObserver,
    error: ParseError,
) -> Result<T> {
    let source_label = source.input_label();
    if matches!(
        &error.kind,
        ParseErrorKind::Io | ParseErrorKind::UnknownFormat | ParseErrorKind::EmptyFile
    ) {
        return Err(error::constructors::right_parser_error(
            source.right_input_origin(),
            source_label,
            &policy,
            observer,
            error,
        ));
    }
    observer.terminal_paired_parser_failure(MateSide::Right, &error, policy)?;
    Err(error::constructors::right_parser_error(
        source.right_input_origin(),
        source_label,
        &policy,
        observer,
        error,
    ))
}

fn parser_failure_has_record_slot(error: &ParseError) -> bool {
    !matches!(
        &error.kind,
        ParseErrorKind::Io | ParseErrorKind::UnknownFormat | ParseErrorKind::EmptyFile
    )
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

    fn admit_record<'record>(
        &self,
        record: RecordView<'record>,
        mate: Option<MateSide>,
        pairs_seen: Option<u64>,
        observer: &mut RunObserver,
    ) -> Result<Option<RecordView<'record>>> {
        if record.sequence().len() == record.quality().len() {
            return Ok(Some(record));
        }

        self.record_length_mismatch(record, mate, pairs_seen, observer)?;

        match self.policy {
            AdmissionPolicy::Error => Err(error::constructors::record_length_error(record, mate)),
            AdmissionPolicy::Skip => {
                Self::warn_length_mismatch(record, mate, observer);
                Ok(None)
            }
        }
    }

    fn record_length_mismatch(
        &self,
        record: RecordView<'_>,
        mate: Option<MateSide>,
        pairs_seen: Option<u64>,
        observer: &mut RunObserver,
    ) -> Result<()> {
        observer.invalid_reads += 1;
        observer.record_invalid_input(InvalidInputEvent::SequenceQualityLengthMismatch {
            source: record.source_display(),
            mate,
            header: String::from_utf8_lossy(record.header()).into_owned(),
            sequence_len: record.sequence().len(),
            quality_len: record.quality().len(),
            reads_seen: observer.reads_seen,
            pairs_seen,
            continued: self.policy == AdmissionPolicy::Skip,
        })?;
        Ok(())
    }

    fn warn_length_mismatch(
        record: RecordView<'_>,
        mate: Option<MateSide>,
        observer: &mut RunObserver,
    ) {
        let mate = match mate {
            Some(MateSide::Left) => "left",
            Some(MateSide::Right) => "right",
            None => "single",
        };
        if observer.should_emit_invalid_input_warning() {
            tracing::warn!(
                source = %record.source_display(),
                mate,
                header = %String::from_utf8_lossy(record.header()),
                sequence_len = record.sequence().len(),
                quality_len = record.quality().len(),
                "skipping input record with mismatched sequence and quality lengths"
            );
        } else if observer.should_emit_invalid_input_suppressed_notice() {
            tracing::warn!("further invalid-input warnings suppressed");
        }
    }
}

impl RecordAdmission<'_, SingleEnd> {
    fn single<'record>(
        &'record self,
        parsed_record: &'record SequenceRecord<'_>,
        observer: &mut RunObserver,
    ) -> Result<Option<RecordView<'record>>> {
        let sequence = parsed_record.raw_seq();
        let quality = parsed_record.qual().ok_or_else(|| {
            error::constructors::missing_single_quality_error(
                self.source.input_origin(),
                self.source.input_label(),
                observer,
            )
        })?;
        let record = RecordView::new(parsed_record.id(), sequence, quality)
            .with_source(self.source.input_source());

        observer.record_seen(sequence.len());

        self.admit_record(record, None, None, observer)
    }
}

impl<'source> RecordAdmission<'source, PairedEnd> {
    fn left_record<'record>(
        &'record self,
        parsed_record: &'record SequenceRecord<'_>,
        observer: &RunObserver,
    ) -> Result<MateRecord<'record, LeftMate>> {
        let sequence = parsed_record.raw_seq();
        let quality = parsed_record.qual().ok_or_else(|| {
            error::constructors::missing_left_quality_error(
                self.source.left_input_origin(),
                self.source.input_label(),
                observer,
            )
        })?;
        Ok(RecordView::new(parsed_record.id(), sequence, quality)
            .with_source(self.source.input_source())
            .into())
    }

    fn right_record<'record>(
        &'record self,
        parsed_record: &'record SequenceRecord<'_>,
        observer: &RunObserver,
    ) -> Result<MateRecord<'record, RightMate>> {
        let sequence = parsed_record.raw_seq();
        let quality = parsed_record.qual().ok_or_else(|| {
            error::constructors::missing_right_quality_error(
                self.source.right_input_origin(),
                self.source.input_label(),
                observer,
            )
        })?;
        Ok(RecordView::new(parsed_record.id(), sequence, quality)
            .with_source(self.source.input_source())
            .into())
    }

    fn missing_right(
        &self,
        parsed_record: &SequenceRecord<'_>,
        observer: &mut RunObserver,
    ) -> Result<()> {
        let record = self.left_record(parsed_record, observer)?;
        observer.record_seen(record.into_record().sequence().len());
        self.missing_right_record(record, observer)
    }

    fn missing_left(
        &self,
        parsed_record: &SequenceRecord<'_>,
        observer: &mut RunObserver,
    ) -> Result<()> {
        let record = self.right_record(parsed_record, observer)?;
        observer.record_seen(record.into_record().sequence().len());
        self.missing_left_record(record, observer)
    }

    fn missing_right_record(
        &self,
        left: MateRecord<'_, LeftMate>,
        observer: &mut RunObserver,
    ) -> Result<()> {
        let record = left.into_record();
        let source = record.source_display();
        let header = String::from_utf8_lossy(record.header()).into_owned();
        let continued = self.policy == AdmissionPolicy::Skip;
        observer.record_invalid_input(InvalidInputEvent::MissingMate {
            source: source.clone(),
            present_mate: MateSide::Left,
            header: header.clone(),
            reads_seen: observer.reads_seen,
            pairs_seen: observer.pairs_seen,
            continued,
        })?;
        match self.policy {
            AdmissionPolicy::Error => Err(error::constructors::left_record_count_error(
                self.source.left_input_origin(),
                self.source.input_layout(),
                source,
                header,
                observer,
            )),
            AdmissionPolicy::Skip => {
                if observer.should_emit_invalid_input_warning() {
                    tracing::warn!(source, header, "skipping left record without a right mate");
                } else if observer.should_emit_invalid_input_suppressed_notice() {
                    tracing::warn!("further invalid-input warnings suppressed");
                }
                Ok(())
            }
        }
    }

    fn missing_left_record(
        &self,
        right: MateRecord<'_, RightMate>,
        observer: &mut RunObserver,
    ) -> Result<()> {
        let record = right.into_record();
        let source = record.source_display();
        let header = String::from_utf8_lossy(record.header()).into_owned();
        let continued = self.policy == AdmissionPolicy::Skip;
        observer.record_invalid_input(InvalidInputEvent::MissingMate {
            source: source.clone(),
            present_mate: MateSide::Right,
            header: header.clone(),
            reads_seen: observer.reads_seen,
            pairs_seen: observer.pairs_seen,
            continued,
        })?;

        match self.policy {
            AdmissionPolicy::Error => Err(error::constructors::right_record_count_error(
                self.source.right_input_origin(),
                self.source.input_layout(),
                source,
                header,
                observer,
            )),
            AdmissionPolicy::Skip => {
                if observer.should_emit_invalid_input_warning() {
                    tracing::warn!(source, header, "skipping right record without a left mate");
                } else if observer.should_emit_invalid_input_suppressed_notice() {
                    tracing::warn!("further invalid-input warnings suppressed");
                }
                Ok(())
            }
        }
    }

    fn buffered_left_record<'record>(
        &'record self,
        parsed_record: &SequenceRecord<'_>,
        observer: &RunObserver,
        buffer: &'record mut InterleavedLeftBuffer,
    ) -> Result<MateRecord<'record, LeftMate>>
    where
        'source: 'record,
    {
        let sequence = parsed_record.raw_seq();
        let quality = parsed_record.qual().ok_or_else(|| {
            error::constructors::missing_left_quality_error(
                self.source.left_input_origin(),
                self.source.input_label(),
                observer,
            )
        })?;

        Ok(buffer
            .copy_from(parsed_record.id(), sequence, quality)
            .with_source(self.source.input_source())
            .into())
    }

    fn admit_pair<'record>(
        &self,
        left: impl Into<MateRecord<'record, LeftMate>>,
        right: impl Into<MateRecord<'record, RightMate>>,
        observer: &mut RunObserver,
    ) -> Result<Option<RecordPair<'record>>> {
        let error = match RecordPair::try_new(left.into(), right.into()) {
            Ok(pair) => return Ok(Some(pair)),
            Err(error) => error,
        };

        observer.invalid_pairs += 1;
        let continued = self.policy == AdmissionPolicy::Skip;
        let source = self.source.input_label();

        match error {
            PairConstructionError::LeftRecordLength { .. }
            | PairConstructionError::RightRecordLength { .. } => observer.invalid_reads += 1,
            PairConstructionError::BothRecordLengths { .. } => observer.invalid_reads += 2,
            PairConstructionError::LeftHeaderClaimsRight { .. }
            | PairConstructionError::RightHeaderClaimsLeft { .. }
            | PairConstructionError::BothHeadersContradictPositions { .. }
            | PairConstructionError::IdentifierMismatch { .. } => {}
        }

        observer.record_invalid_input(InvalidInputEvent::PairConstructionFailure {
            source: source.clone(),
            error,
            reads_seen: observer.reads_seen,
            pairs_seen: observer.pairs_seen,
            continued,
        })?;

        match self.policy {
            AdmissionPolicy::Skip => {
                if observer.should_emit_invalid_input_warning() {
                    tracing::warn!(
                        source,
                        pair_error = %error,
                        "skipping paired input that failed construction"
                    );
                } else if observer.should_emit_invalid_input_suppressed_notice() {
                    tracing::warn!("further invalid-input warnings suppressed");
                }
                Ok(None)
            }
            AdmissionPolicy::Error => Err(error::constructors::pair_construction_error(
                self.source.left_input_origin(),
                source,
                error,
            )),
        }
    }
}

fn run_observer(config: &RunConfig, source_label: String) -> Result<RunObserver> {
    let mut observer = RunObserver::new(source_label);
    if let Some(path) = &config.invalid_input_report {
        observer.set_invalid_input_report(InvalidInputReport::create(path)?);
    }
    Ok(observer)
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
        observer::{InvalidInputEvent, InvalidInputReport, RunObserver},
        output::{
            InterleavedOutput, OutputArgs, OutputEncoding, OutputFormat, PairedRecordOutput,
            SingleOutput, SingleRecordOutput, StreamSink,
        },
        record::{MateSide, RecordView},
    };

    use super::{
        InputOrigin, PairedEnd, PairedEndContext, PairedSource, RecordAdmission, RunConfig,
        SingleEnd, SingleSource, single_parser_failure,
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
            progress_every: 1_234,
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
        let mut observer = RunObserver::new(format!("ena:{accession}"));
        observer.reads_seen = 7;
        observer.pairs_seen = 3;

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
            let error = crate::error::constructors::single_parser_error(
                InputOrigin::Ena(&accession),
                format!("ena:{accession}"),
                &AdmissionPolicy::Error,
                &observer,
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
        let error = crate::error::constructors::single_parser_error(
            InputOrigin::Local(local_path),
            "local:reads.fastq".to_owned(),
            &AdmissionPolicy::Error,
            &observer,
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
            let error = crate::error::constructors::single_parser_error(
                InputOrigin::Local(local_path),
                "local:reads.fastq".to_owned(),
                &AdmissionPolicy::Error,
                &observer,
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
        let mut observer = RunObserver::new("local:reads.fastq".to_owned());
        observer.reads_seen = 7;
        observer.pairs_seen = 3;
        let local_path = Path::new("reads.fastq");

        for source in [
            ParseError::new_unknown_format(b'!'),
            ParseError::new_empty_file(),
        ] {
            let expected_kind = source.kind.clone();
            let error = crate::error::constructors::single_parser_error(
                InputOrigin::Local(local_path),
                "local:reads.fastq".to_owned(),
                &AdmissionPolicy::Error,
                &observer,
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
        let mut observer = RunObserver::new("test".to_owned());
        observer.reads_seen = 3;
        observer.pairs_seen = 1;
        let accession = Accession::new("SRR35939766")?;
        let ena_source = PairedSource::Ena { accession };
        let interleaved_source = PairedSource::LocalInterleaved {
            input: "reads.fastq".into(),
        };
        let split_source = PairedSource::LocalSplit {
            input1: "reads_1.fastq".into(),
            input2: "reads_2.fastq".into(),
        };

        assert!(matches!(
            crate::error::constructors::left_record_count_error(
                ena_source.left_input_origin(),
                ena_source.input_layout(),
                ena_source.input_label(),
                "read2/1".to_owned(),
                &observer,
            ),
            RunError::IndeterminateInput(IndeterminateInputError::Content {
                problem: EnaContentProblem::RecordCount {
                    complete_pairs_seen: 1,
                    present_mate: MateSide::Left,
                    ref header,
                },
                ..
            }) if header == "read2/1"
        ));
        assert!(matches!(
            crate::error::constructors::left_record_count_error(
                interleaved_source.left_input_origin(),
                interleaved_source.input_layout(),
                interleaved_source.input_label(),
                "read2/1".to_owned(),
                &observer,
            ),
            RunError::MalformedInput(error)
                if matches!(*error, MalformedInputError::InterleavedRecordCount { .. })
        ));
        assert!(matches!(
            crate::error::constructors::right_record_count_error(
                split_source.right_input_origin(),
                split_source.input_layout(),
                split_source.input_label(),
                "read2/2".to_owned(),
                &observer,
            ),
            RunError::MalformedInput(error)
                if matches!(*error, MalformedInputError::PairedRecordCount { .. })
        ));
        Ok(())
    }

    #[test]
    fn missing_quality_errors_preserve_input_origin() -> Result<()> {
        let mut observer = RunObserver::new("test".to_owned());
        observer.reads_seen = 7;
        observer.pairs_seen = 3;
        let accession = Accession::new("SRR35939766")?;
        let ena_source = SingleSource::Ena { accession };
        let local_source = SingleSource::Local {
            input: "reads.fastq".into(),
        };

        assert!(matches!(
            crate::error::constructors::missing_single_quality_error(
                ena_source.input_origin(),
                ena_source.input_label(),
                &observer,
            ),
            RunError::IndeterminateInput(IndeterminateInputError::Content {
                problem: EnaContentProblem::MissingQuality { mate: "single", .. },
                ..
            })
        ));
        assert!(matches!(
            crate::error::constructors::missing_single_quality_error(
                local_source.input_origin(),
                local_source.input_label(),
                &observer,
            ),
            RunError::MalformedInput(error)
                if matches!(*error, MalformedInputError::MissingQuality { .. })
        ));
        Ok(())
    }

    #[test]
    fn record_length_errors_preserve_input_origin() -> Result<()> {
        let local_source = SingleSource::Local {
            input: PathBuf::from("reads.fastq"),
        };
        let local_admission =
            RecordAdmission::<SingleEnd>::new(&local_source, AdmissionPolicy::Error);
        let local_record =
            RecordView::new(b"read1", b"ACGT", b"I").with_source(local_source.input_source());
        let mut local_observer = RunObserver::new(local_source.input_label());
        local_observer.record_seen(local_record.sequence().len());

        let Err(local_error) =
            local_admission.admit_record(local_record, None, None, &mut local_observer)
        else {
            panic!("local length mismatch should fail admission");
        };
        assert!(matches!(
            local_error,
            RunError::MalformedInput(error)
                if matches!(*error, MalformedInputError::RecordLength { .. })
        ));
        assert_eq!(local_observer.invalid_reads, 1);
        assert_eq!(
            local_observer
                .invalid_input_event_counts()
                .get("sequence_quality_length_mismatch"),
            Some(&1)
        );

        let ena_source = SingleSource::Ena {
            accession: Accession::new("SRR35939766")?,
        };
        let ena_admission = RecordAdmission::<SingleEnd>::new(&ena_source, AdmissionPolicy::Error);
        let ena_record =
            RecordView::new(b"read1", b"ACGT", b"I").with_source(ena_source.input_source());
        let mut ena_observer = RunObserver::new(ena_source.input_label());
        ena_observer.record_seen(ena_record.sequence().len());

        let Err(ena_error) = ena_admission.admit_record(ena_record, None, None, &mut ena_observer)
        else {
            panic!("ENA length mismatch should fail admission");
        };
        assert!(matches!(
            ena_error,
            RunError::IndeterminateInput(IndeterminateInputError::Content {
                problem: EnaContentProblem::RecordLength { .. },
                ..
            })
        ));
        assert_eq!(ena_observer.invalid_reads, 1);
        Ok(())
    }

    #[test]
    fn pair_identifier_errors_preserve_input_origin() -> Result<()> {
        let local_source = PairedSource::LocalInterleaved {
            input: PathBuf::from("reads.fastq"),
        };
        let local_admission =
            RecordAdmission::<PairedEnd>::new(&local_source, AdmissionPolicy::Error);
        let local_left =
            RecordView::new(b"read1/1", b"A", b"I").with_source(local_source.input_source());
        let local_right =
            RecordView::new(b"read2/2", b"T", b"I").with_source(local_source.input_source());
        let mut local_observer = RunObserver::new(local_source.input_label());
        local_observer.record_seen(1);
        local_observer.record_seen(1);

        let Err(local_error) =
            local_admission.admit_pair(local_left, local_right, &mut local_observer)
        else {
            panic!("local mate mismatch should fail admission");
        };
        assert!(matches!(
            local_error,
            RunError::MalformedInput(error)
                if matches!(*error, MalformedInputError::PairConstruction { .. })
        ));
        assert_eq!(local_observer.invalid_pairs, 1);

        let ena_source = PairedSource::Ena {
            accession: Accession::new("SRR35939766")?,
        };
        let ena_admission = RecordAdmission::<PairedEnd>::new(&ena_source, AdmissionPolicy::Error);
        let ena_left =
            RecordView::new(b"read1/1", b"A", b"I").with_source(ena_source.input_source());
        let ena_right =
            RecordView::new(b"read2/2", b"T", b"I").with_source(ena_source.input_source());
        let mut ena_observer = RunObserver::new(ena_source.input_label());
        ena_observer.record_seen(1);
        ena_observer.record_seen(1);

        let Err(ena_error) = ena_admission.admit_pair(ena_left, ena_right, &mut ena_observer)
        else {
            panic!("ENA mate mismatch should fail admission");
        };
        assert!(matches!(
            ena_error,
            RunError::IndeterminateInput(IndeterminateInputError::Content {
                problem: EnaContentProblem::PairConstruction { .. },
                ..
            })
        ));
        assert_eq!(ena_observer.invalid_pairs, 1);
        Ok(())
    }

    #[test]
    fn parser_io_failure_is_not_an_invalid_input_event() -> Result<()> {
        let temp = tempdir()?;
        let report_path = temp.path().join("invalid-input.jsonl");
        let source = SingleSource::Local {
            input: PathBuf::from("reads.fastq"),
        };
        let mut observer = RunObserver::new(source.input_label());
        observer.set_invalid_input_report(InvalidInputReport::create(&report_path)?);

        let error = single_parser_failure::<()>(
            &source,
            AdmissionPolicy::Skip,
            &mut observer,
            ParseError::from(io::Error::other("synthetic read failure")),
        )
        .expect_err("I/O failure should remain fatal");

        assert!(
            matches!(error, RunError::Io(IoError::LocalFastqRead { .. })),
            "I/O parser failure should be classified as local read I/O"
        );
        assert_eq!(observer.invalid_reads, 0);
        assert!(observer.invalid_input_event_counts().is_empty());
        assert_eq!(std::fs::read_to_string(report_path)?, "");
        Ok(())
    }

    #[test]
    fn parser_source_classification_failures_are_not_invalid_input_events() -> Result<()> {
        let temp = tempdir()?;
        let report_path = temp.path().join("invalid-input.jsonl");
        let source = SingleSource::Local {
            input: PathBuf::from("reads.fastq"),
        };
        let mut observer = RunObserver::new(source.input_label());
        observer.set_invalid_input_report(InvalidInputReport::create(&report_path)?);

        for parse_error in [
            ParseError::new_unknown_format(b'!'),
            ParseError::new_empty_file(),
        ] {
            let error = single_parser_failure::<()>(
                &source,
                AdmissionPolicy::Skip,
                &mut observer,
                parse_error,
            )
            .expect_err("source classification failure should remain fatal");
            assert!(
                error
                    .to_string()
                    .contains("input source did not provide a readable FASTQ stream")
            );
        }

        assert_eq!(observer.invalid_reads, 0);
        assert!(observer.invalid_input_event_counts().is_empty());
        assert_eq!(std::fs::read_to_string(report_path)?, "");
        Ok(())
    }

    #[test]
    fn malformed_pair_handling_is_atomic_for_each_invalid_mate() -> Result<()> {
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
                    .with_source(source.input_source());
                let right = RecordView::new(b"read1/2", b"TTTT", right_quality)
                    .with_source(source.input_source());
                let mut observer = RunObserver::new(source.input_label());
                observer.pairs_seen = 1;
                observer.record_seen(left.sequence().len());
                observer.record_seen(right.sequence().len());

                let result = admission.admit_pair(left, right, &mut observer);
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

                assert_eq!(observer.pairs_seen, 1);
                assert_eq!(observer.invalid_pairs, 1);
                assert_eq!(observer.invalid_reads, 1);
                let expected_code = "sequence_quality_length_mismatch";
                assert_eq!(
                    observer.invalid_input_event_counts().get(expected_code),
                    Some(&1)
                );
                let [
                    InvalidInputEvent::PairConstructionFailure {
                        error, continued, ..
                    },
                ] = observer.invalid_input_samples()
                else {
                    panic!("exactly one length-mismatch event should identify the invalid mate");
                };
                assert_eq!(error.code(), expected_code);
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
            let left =
                RecordView::new(b"read1/1", b"AAAA", b"I").with_source(source.input_source());
            let right =
                RecordView::new(b"read1/2", b"TTTT", b"K").with_source(source.input_source());
            let mut observer = RunObserver::new(source.input_label());
            observer.pairs_seen = 1;
            observer.record_seen(left.sequence().len());
            observer.record_seen(right.sequence().len());

            let result = admission.admit_pair(left, right, &mut observer);
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

            assert_eq!(observer.pairs_seen, 1);
            assert_eq!(observer.invalid_pairs, 1);
            assert_eq!(observer.invalid_reads, 2);
            assert_eq!(
                observer
                    .invalid_input_event_counts()
                    .get("sequence_quality_length_mismatch"),
                Some(&1)
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
