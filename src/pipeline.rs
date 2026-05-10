//! Top-level ingress, parsing, and output orchestration.

use std::{
    any::Any,
    fs::File,
    marker::PhantomData,
    panic::{AssertUnwindSafe, catch_unwind},
    path::PathBuf,
    time::Instant,
};

use color_eyre::eyre::{Result, WrapErr, bail, eyre};
use needletail::{errors::ParseError, parse_fastx_reader, parser::SequenceRecord};

use crate::{
    adapter::{AdapterPreset, TrimAdaptersTransform},
    cli::{Cli, Ingress, InvalidFastqPolicy, UiPolicy},
    ena::{Accession, EnaClient, FastqUrlsByLayout},
    filter::{MaxNsFilter, MinEntropyFilter, MinLengthFilter, MinMeanQualityFilter},
    output::{OutputArgs, PairedOutputHandle, SingleOutputHandle, UnitOutput},
    pair_merge::MergePairsTransform,
    plan::{
        BuildPlan, Execute, Execution, Logical, OrphanPolicy, Plan, RecordPair, TransformArena,
    },
    progress::ProgressReporter,
    quality::QualityTrimTransform,
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
    type Readers = (Box<dyn std::io::Read + Send>, Box<dyn std::io::Read + Send>);
    type Output = PairedOutputHandle;
}

trait RunSource {
    fn input_label(&self) -> String;
    fn summary_context(&self) -> RunSummaryContext;
}

enum SingleSource {
    Ena { accession: String },
    Local { input: PathBuf },
}

enum PairedSource {
    Ena { accession: String },
    Local { input1: PathBuf, input2: PathBuf },
}

impl RunSource for SingleSource {
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
                accession: Some(accession.clone()),
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
                source: InputSource::Ena { accession },
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
            Self::Local { input1, input2 } => {
                format!("local-paired:{}|{}", input1.display(), input2.display())
            }
        }
    }

    fn summary_context(&self) -> RunSummaryContext {
        match self {
            Self::Ena { accession } => RunSummaryContext {
                ingress_mode: report::IngressMode::Ena,
                layout: RunLayout::Paired,
                accession: Some(accession.clone()),
                input1: None,
                input2: None,
            },
            Self::Local { input1, input2 } => RunSummaryContext {
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
                source: InputSource::Ena { accession },
                mate: Some(mate),
            },
            Self::Local { input1, input2 } => RecordProvenance {
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
    adapter_preset: AdapterPreset,
    merge_pairs: bool,
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
            adapter_preset: cli.adapter_preset,
            merge_pairs: cli.merge_pairs,
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
            bail!(
                "--merge-pairs requires paired-end input\n\
                 help: provide --in1 and --in2, or use an ENA accession with paired FASTQ files"
            );
        }

        Ok(())
    }

    fn build_plan(&self, layout: RunLayout) -> Result<Plan<Execution>> {
        self.validate_layout(layout)?;

        let plan = Plan::<Logical>::new();
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

        Ok(plan
            .quality_trim(self.trim_min_q)
            .min_length(self.min_length)
            .min_mean_q(self.min_mean_q)
            .min_entropy(self.min_entropy)
            .orphan_policy(OrphanPolicy::DropPair)
            .compile())
    }
}

/// Run the CLI-selected ingress path to completion.
///
/// # Errors
///
/// Returns an error when ingress resolution, reader construction, parsing, writing, or output
/// finalization fails.
pub fn run(cli: &Cli) -> Result<()> {
    let config = RunConfig::from(cli);
    let ui = cli.ui_policy();

    match cli.ingress().wrap_err(
        "invalid input selection\nhelp: choose exactly one ingress mode: --ena ACCESSION, --in1 FASTQ, or --in1 FASTQ --in2 FASTQ",
    )? {
        Ingress::LocalSingle { fastq } => {
            config.validate_layout(RunLayout::Single)?;
            SingleEndContext::open_local(fastq, cli.output_args())?.run(&config, &ui)
        }
        Ingress::LocalPaired { r1, r2 } => {
            config.validate_layout(RunLayout::Paired)?;
            PairedEndContext::open_local(r1, r2, cli.output_args())?.run(&config, &ui)
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
    let client =
        EnaClient::new().wrap_err("failed to construct ENA HTTP client for FASTQ streaming")?;
    let layout = client.lookup_fastq_urls(accession).wrap_err_with(|| {
        format!(
            "failed to resolve ENA FASTQ URLs for accession {accession}\n\
             help: confirm this is a run accession with public FASTQ files in ENA"
        )
    })?;

    match layout {
        FastqUrlsByLayout::Single(url) => {
            config.validate_layout(RunLayout::Single)?;
            SingleEndContext::open_ena(accession, client.open_retrying_stream(url), output_args)?
                .run(config, ui)
        }
        FastqUrlsByLayout::Paired(urls) => {
            config.validate_layout(RunLayout::Paired)?;
            let (r1, r2) = client.open_retrying_paired_streams(urls);
            PairedEndContext::open_ena(accession, r1, r2, output_args)?.run(config, ui)
        }
    }
}

impl RunContext<SingleEnd> {
    fn open_local(fastq: PathBuf, output_args: OutputArgs) -> Result<Self> {
        let reader = File::open(&fastq).wrap_err_with(|| {
            format!(
                "failed to open local single-end FASTQ input: {}\n\
                 help: check that --in1 points to a readable FASTQ file on this machine",
                fastq.display()
            )
        })?;
        let output = output_args
            .resolve_single()
            .wrap_err("failed to configure single-end output")?;
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
        let output = output_args
            .resolve_single()
            .wrap_err("failed to configure single-end output")?;
        Ok(Self {
            source: SingleSource::Ena {
                accession: accession.to_string(),
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

        // build out the mutable state needed to run the application loop
        let mut parser = parse_fastx_reader(reader).wrap_err(
        "failed to initialize FASTQ parser for single-end input\nhelp: confirm the input is readable FASTQ, optionally gzip-compressed if the parser supports it",
    )?;
        let mut arena = TransformArena::new();
        let mut stats = read_stats(config)?;
        let mut progress = ProgressReporter::new(ui.progress_mode, config.progress_every);
        let started_at = Instant::now();
        let input_label = source.input_label();
        let admission = FastqAdmission::<SingleEnd>::new(&source, config.invalid_fastq_policy);

        while let Some(next_record) =
            catch_parser_panic(&input_label, "single", &stats, || parser.next())?
        {
            let parsed_record = admission.parse("single", &mut stats, next_record)?;
            let Some(record) = admission.single(&parsed_record, &mut stats)? else {
                continue;
            };
            arena.reset();

            let outcome = plan.execute(record, &mut arena, &mut stats)?;
            output.write_outcome(&outcome, &mut stats).wrap_err(
                "failed to write single-end output record\nhelp: check downstream pipe or output filesystem health",
            )?;

            progress.maybe_report(&stats);
        }

        progress.finish();
        output.finish().wrap_err(
            "failed to finalize single-end output\nhelp: for gzip output this can indicate a truncated destination or broken downstream pipe",
        )?;
        let summary =
            report::RunSummary::from_stats(source.summary_context(), &stats, started_at.elapsed());
        if ui.show_summary {
            report::print_summary(&summary);
        }
        if let Some(path) = &config.summary {
            report::write_summary_json(path, &summary).wrap_err_with(|| {
            format!(
                "failed to write JSON run summary: {}\nhelp: check that the parent directory exists and is writable",
                path.display()
            )
        })?;
        }
        Ok(())
    }
}

impl RunContext<PairedEnd> {
    fn open_local(r1: PathBuf, r2: PathBuf, output_args: OutputArgs) -> Result<Self> {
        let reader1 = File::open(&r1).wrap_err_with(|| {
            format!(
                "failed to open local paired FASTQ input for read 1: {}\n\
                 help: check that --in1 is readable from the current execution environment",
                r1.display()
            )
        })?;
        let reader2 = File::open(&r2).wrap_err_with(|| {
            format!(
                "failed to open local paired FASTQ input for read 2: {}\n\
                 help: check that --in2 is readable from the current execution environment",
                r2.display()
            )
        })?;
        let output = output_args
            .resolve_paired()
            .wrap_err("failed to configure paired-end output")?;
        Ok(Self {
            source: PairedSource::Local {
                input1: r1,
                input2: r2,
            },
            readers: (Box::new(reader1), Box::new(reader2)),
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
        let output = output_args
            .resolve_paired()
            .wrap_err("failed to configure paired-end output")?;
        Ok(Self {
            source: PairedSource::Ena {
                accession: accession.to_string(),
            },
            readers: (Box::new(r1), Box::new(r2)),
            output,
            _layout: PhantomData,
        })
    }

    fn run(self, config: &RunConfig, ui: &UiPolicy) -> Result<()> {
        let Self {
            source,
            readers: (r1, r2),
            mut output,
            _layout,
        } = self;
        let mut plan = config.build_plan(RunLayout::Paired)?;
        let mut parser_r1 = parse_fastx_reader(r1).wrap_err(
        "failed to initialize FASTQ parser for read 1\nhelp: confirm --in1 is readable FASTQ and has the expected compression",
    )?;
        let mut parser_r2 = parse_fastx_reader(r2).wrap_err(
        "failed to initialize FASTQ parser for read 2\nhelp: confirm --in2 is readable FASTQ and has the expected compression",
    )?;
        let mut arena = TransformArena::new();
        let mut stats = read_stats(config)?;
        let mut progress = ProgressReporter::new(ui.progress_mode, config.progress_every);
        let started_at = Instant::now();
        let input_label = source.input_label();
        let admission = FastqAdmission::<PairedEnd>::new(&source, config.invalid_fastq_policy);

        loop {
            let next_r1 = catch_parser_panic(&input_label, "left", &stats, || parser_r1.next())?;
            let next_r2 = catch_parser_panic(&input_label, "right", &stats, || parser_r2.next())?;

            match (next_r1, next_r2) {
                (Some(record_r1), Some(record_r2)) => {
                    let parsed_r1 = admission.parse("left", &mut stats, record_r1)?;
                    let parsed_r2 = admission.parse("right", &mut stats, record_r2)?;
                    let Some(pair) = admission.pair(&parsed_r1, &parsed_r2, &mut stats)? else {
                        continue;
                    };

                    arena.reset();

                    let outcome = plan.execute(pair, &mut arena, &mut stats)?;
                    output.write_outcome(&outcome, &mut stats).wrap_err(
                        "failed to write paired output record group\nhelp: check downstream pipe or output filesystem health",
                    )?;

                    progress.maybe_report(&stats);
                }
                (None, None) => break,
                _ => bail!(
                    "paired FASTQ inputs have different record counts\n\
                 source: {}\n\
                 complete_pairs_seen: {}\n\
                 reads_seen_before_failure: {}\n\
                 help: confirm --in1 and --in2 are mates from the same run and were not independently filtered or truncated",
                    input_label,
                    stats.pairs_seen,
                    stats.reads_seen,
                ),
            }
        }

        progress.finish();
        output.finish().wrap_err(
        "failed to finalize paired output\nhelp: for gzip output this can indicate a truncated destination or broken downstream pipe",
    )?;
        let summary =
            report::RunSummary::from_stats(source.summary_context(), &stats, started_at.elapsed());
        if ui.show_summary {
            report::print_summary(&summary);
        }
        if let Some(path) = &config.summary {
            report::write_summary_json(path, &summary).wrap_err_with(|| {
            format!(
                "failed to write JSON run summary: {}\nhelp: check that the parent directory exists and is writable",
                path.display()
            )
        })?;
        }
        Ok(())
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
        next_record: Result<SequenceRecord<'record>, ParseError>,
    ) -> Result<SequenceRecord<'record>> {
        match next_record {
            Ok(record) => Ok(record),
            Err(error) => self.parser_error(mate, stats, &error),
        }
    }

    fn parser_error<T>(
        &self,
        mate: &'static str,
        stats: &mut ReadStats,
        error: &ParseError,
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

        bail!(
            "FASTQ parser rejected malformed input while reading source={source} mate={mate}\n\
             invalid_fastq_policy={}\n\
             reads_seen={} pairs_seen={} invalid_reads={} invalid_pairs={}\n\
             parser_error={parser_error_message}\n\
             parser_error_kind={parser_error_kind}\n\
             help: malformed FASTQ parse errors are not currently recoverable; check input integrity and retry ENA-backed reads if the stream may have been interrupted",
            self.policy,
            stats.reads_seen,
            stats.pairs_seen,
            stats.invalid_reads,
            stats.invalid_pairs,
        )
    }
}

impl FastqAdmission<'_, SingleEnd> {
    fn single<'record>(
        &'record self,
        parsed_record: &'record SequenceRecord<'_>,
        stats: &mut ReadStats,
    ) -> Result<Option<RecordView<'record>>> {
        let source = self.source.input_label();
        let sequence = parsed_record.raw_seq();
        let quality = parsed_record
            .qual()
            .ok_or_else(|| missing_quality_error(&source, "single", stats))?;
        let record = RecordView::new(parsed_record.id(), sequence, quality)
            .with_provenance(self.source.provenance());

        stats.record_seen(sequence.len());

        record.validate(self.policy, stats)
    }
}

impl FastqAdmission<'_, PairedEnd> {
    fn pair<'record>(
        &'record self,
        parsed_r1: &'record SequenceRecord<'_>,
        parsed_r2: &'record SequenceRecord<'_>,
        stats: &mut ReadStats,
    ) -> Result<Option<RecordPair<'record>>> {
        let source = self.source.input_label();
        let sequence_r1 = parsed_r1.raw_seq();
        let sequence_r2 = parsed_r2.raw_seq();
        let quality_r1 = parsed_r1
            .qual()
            .ok_or_else(|| missing_quality_error(&source, "left", stats))?;
        let quality_r2 = parsed_r2
            .qual()
            .ok_or_else(|| missing_quality_error(&source, "right", stats))?;
        let left = RecordView::new(parsed_r1.id(), sequence_r1, quality_r1)
            .with_provenance(self.source.provenance(MateSide::Left));
        let right = RecordView::new(parsed_r2.id(), sequence_r2, quality_r2)
            .with_provenance(self.source.provenance(MateSide::Right));

        stats.record_seen(sequence_r1.len());
        stats.record_seen(sequence_r2.len());
        stats.pairs_seen += 1;

        left.validate_pair(right, self.policy, stats)
    }
}

fn catch_parser_panic<T>(
    source: &str,
    mate: &str,
    stats: &ReadStats,
    operation: impl FnOnce() -> T,
) -> Result<T> {
    catch_unwind(AssertUnwindSafe(operation)).map_err(|panic| {
        eyre!(
            "FASTQ parser failed while reading source={source} mate={mate}\n\
             reads_seen={} pairs_seen={} invalid_reads={} invalid_pairs={}\n\
             panic={}\n\
             help: the input stream appears desynchronized; retry ENA accessions and inspect the invalid FASTQ report if one was configured",
            stats.reads_seen,
            stats.pairs_seen,
            stats.invalid_reads,
            stats.invalid_pairs,
            panic_message(&panic),
        )
    })
}

fn missing_quality_error(source: &str, mate: &str, stats: &ReadStats) -> color_eyre::Report {
    eyre!(
        "FASTQ parser did not provide quality scores while reading source={source} mate={mate}\n\
         reads_seen={} pairs_seen={} invalid_reads={} invalid_pairs={}\n\
         help: confirm the input is FASTQ rather than FASTA and that parser quality computation is enabled",
        stats.reads_seen,
        stats.pairs_seen,
        stats.invalid_reads,
        stats.invalid_pairs,
    )
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
        stats.set_invalid_fastq_report(InvalidFastqReport::create(path).wrap_err_with(|| {
            format!(
                "failed to create invalid FASTQ JSONL report: {}\nhelp: check that the parent directory exists and is writable",
                path.display()
            )
        })?);
    }
    Ok(stats)
}

#[cfg(test)]
mod tests {
    use std::{fs::File, io::Cursor, marker::PhantomData, path::Path};

    use color_eyre::{Result, eyre::bail};
    use needletail::parse_fastx_reader;
    use tempfile::tempdir;

    use crate::{
        adapter::AdapterPreset,
        cli::{Cli, InvalidFastqPolicy},
        output::{
            InterleavedOutput, OutputArgs, OutputEncoding, OutputFormat, PairedRecordOutput,
            SingleOutput, SingleRecordOutput, StreamSink,
        },
        record::RecordView,
    };

    use super::{PairedEndContext, PairedSource, RunConfig};

    fn single_output_for_vec(format: OutputFormat) -> SingleOutput<StreamSink<Vec<u8>>> {
        SingleOutput::new(StreamSink::new(Vec::new(), format))
    }

    fn interleaved_output_for_vec(format: OutputFormat) -> InterleavedOutput<StreamSink<Vec<u8>>> {
        InterleavedOutput::new(StreamSink::new(Vec::new(), format))
    }

    fn test_cli() -> Cli {
        Cli {
            ena: None,
            in1: None,
            in2: None,
            min_length: 50,
            max_ns: 4,
            min_mean_q: 20.0,
            trim_min_q: 20,
            adapter_preset: AdapterPreset::IlluminaTruSeq,
            merge_pairs: false,
            merge_min_overlap: 10,
            merge_max_mismatch_rate: 0.2,
            merge_min_correction_delta_q: 0,
            min_entropy: 0.0,
            interleaved: true,
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
            true,
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
            source: PairedSource::Local {
                input1: "reads_1.fastq.gz".into(),
                input2: "reads_2.fastq.gz".into(),
            },
            readers: (Box::new(r1), Box::new(r2)),
            output,
            _layout: PhantomData,
        }
        .run(&config, &ui)
        .expect_err("mismatched paired inputs should fail");

        assert!(
            error
                .to_string()
                .contains("paired FASTQ inputs have different record counts")
        );
        Ok(())
    }

    fn write_fixture(path: &Path, bytes: &[u8]) -> Result<()> {
        use std::io::Write as _;

        let mut file = File::create(path)?;
        file.write_all(bytes)?;
        Ok(())
    }
}
