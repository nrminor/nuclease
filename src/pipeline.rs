//! Top-level ingress, parsing, and output orchestration.

use std::{
    any::Any,
    io::Read,
    panic::{AssertUnwindSafe, catch_unwind},
    time::Instant,
};

use color_eyre::eyre::{Result, WrapErr, bail, eyre};
use helicase::{
    Config, FastqParser, HelicaseParser, ParserOptions,
    input::{FromReader, InputData},
};

use crate::{
    adapter::{AdapterCatalog, TrimAdaptersTransform},
    cli::{Cli, IngressHandle},
    filter::{MaxNsFilter, MinEntropyFilter, MinLengthFilter, MinMeanQualityFilter},
    output::{PairedOutputHandle, SingleOutputHandle, SingleRecordOutput},
    plan::{BuildPlan, Execute, Logical, OrphanPolicy, Plan, RecordPair, TransformArena},
    progress::ProgressReporter,
    quality::QualityTrimTransform,
    record::{InputSource, InvalidFastqReport, MateSide, ReadStats, RecordProvenance, RecordView},
    report::{self, RunContext, RunLayout},
};

const FASTQ_CONFIG: Config = ParserOptions::default().compute_quality().config();

/// Run the CLI-selected ingress path to completion.
///
/// # Errors
///
/// Returns an error when ingress resolution, reader construction, parsing, writing, or output
/// finalization fails.
pub fn run(cli: &Cli) -> Result<()> {
    // work out data ingress based on the information provided by the user in the
    // command line arguments
    let ingress = cli.ingress().wrap_err(
        "invalid input selection\nhelp: choose exactly one ingress mode: --ena ACCESSION, --in1 FASTQ, or --in1 FASTQ --in2 FASTQ",
    )?;
    let summary_context = RunContext::from_ingress(&ingress);
    let ui = cli.ui_policy();

    // open reader(s) for the provided FASTQ(s) or ENA accession(s), build correct
    // writers based on their layout, and let it rip
    match ingress
        .open()
        .wrap_err("failed to open selected FASTQ ingress")?
    {
        IngressHandle::Single(reader) => {
            let output = cli
                .output_args()
                .resolve_single()
                .wrap_err("failed to configure single-end output")?;
            let ctx = summary_context.with_layout(RunLayout::Single);
            run_single(reader, output, cli, &ui, ctx)
        }
        IngressHandle::Paired(r1, r2) => {
            let output = cli
                .output_args()
                .resolve_paired()
                .wrap_err("failed to configure paired-end output")?;
            let ctx = summary_context.with_layout(RunLayout::Paired);
            run_paired(r1, r2, output, cli, &ui, ctx)
        }
    }
}

fn run_single(
    reader: impl Read + Send,
    mut output: SingleOutputHandle,
    cli: &Cli,
    ui: &crate::cli::UiPolicy,
    summary_context: RunContext,
) -> Result<()> {
    // build a plan for what to do with the reads based on the user's provided settings
    let plan = Plan::new()
        .max_ns(cli.max_ns)
        .trim_adapters(AdapterCatalog::illumina_truseq())
        .quality_trim(cli.trim_min_q)
        .min_length(cli.min_length)
        .min_mean_q(cli.min_mean_q)
        .min_entropy(cli.min_entropy)
        .orphan_policy(OrphanPolicy::DropPair)
        .compile();

    // build out the mutable state needed to run the application loop
    let mut parser = FastqParser::<FASTQ_CONFIG, _>::from_reader(reader).wrap_err(
        "failed to initialize FASTQ parser for single-end input\nhelp: confirm the input is readable FASTQ, optionally gzip-compressed if the parser supports it",
    )?;
    let mut arena = TransformArena::new();
    let mut stats = read_stats(cli)?;
    let mut progress = ProgressReporter::new(ui.progress_mode, cli.progress_every);
    let started_at = Instant::now();

    while catch_parser_panic(&summary_context.input_label(), "single", &stats, || {
        parser.next()
    })?
    .is_some()
    {
        let Some(record) = read_validated_single(
            &parser,
            &summary_context,
            cli.invalid_fastq_policy,
            &mut stats,
        )?
        else {
            continue;
        };
        arena.reset();

        let outcome = plan.execute(record, &mut arena, &mut stats);
        if outcome.rejection_count() == 0 {
            for record in outcome.emitted() {
                output.write_record(record).wrap_err_with(|| {
                    format!(
                        "failed to write single-end output record\nheader: {}\nhelp: check downstream pipe or output filesystem health",
                        String::from_utf8_lossy(record.header())
                    )
                })?;
                stats.record_emitted(record.sequence().len());
            }
        }

        progress.maybe_report(&stats);
    }

    progress.finish();
    output.finish().wrap_err(
        "failed to finalize single-end output\nhelp: for gzip output this can indicate a truncated destination or broken downstream pipe",
    )?;
    let summary = report::RunSummary::from_stats(summary_context, &stats, started_at.elapsed());
    if ui.show_summary {
        report::print_summary(&summary);
    }
    if let Some(path) = &cli.summary {
        report::write_summary_json(path, &summary).wrap_err_with(|| {
            format!(
                "failed to write JSON run summary: {}\nhelp: check that the parent directory exists and is writable",
                path.display()
            )
        })?;
    }
    Ok(())
}

fn run_paired(
    r1: impl Read + Send,
    r2: impl Read + Send,
    mut output: PairedOutputHandle,
    cli: &Cli,
    ui: &crate::cli::UiPolicy,
    summary_context: RunContext,
) -> Result<()> {
    let plan = Plan::<Logical>::new()
        .max_ns(cli.max_ns)
        .trim_adapters(AdapterCatalog::illumina_truseq())
        .quality_trim(cli.trim_min_q)
        .min_length(cli.min_length)
        .min_mean_q(cli.min_mean_q)
        .min_entropy(cli.min_entropy)
        .orphan_policy(OrphanPolicy::DropPair)
        .compile();
    let mut parser_r1 = FastqParser::<FASTQ_CONFIG, _>::from_reader(r1).wrap_err(
        "failed to initialize FASTQ parser for read 1\nhelp: confirm --in1 is readable FASTQ and has the expected compression",
    )?;
    let mut parser_r2 = FastqParser::<FASTQ_CONFIG, _>::from_reader(r2).wrap_err(
        "failed to initialize FASTQ parser for read 2\nhelp: confirm --in2 is readable FASTQ and has the expected compression",
    )?;
    let mut arena = TransformArena::new();
    let mut stats = read_stats(cli)?;
    let mut progress = ProgressReporter::new(ui.progress_mode, cli.progress_every);
    let started_at = Instant::now();

    loop {
        let next_r1 = catch_parser_panic(&summary_context.input_label(), "left", &stats, || {
            parser_r1.next()
        })?;
        let next_r2 = catch_parser_panic(&summary_context.input_label(), "right", &stats, || {
            parser_r2.next()
        })?;

        match (next_r1, next_r2) {
            (Some(_), Some(_)) => {
                let Some(pair) = read_validated_pair(
                    &parser_r1,
                    &parser_r2,
                    &summary_context,
                    cli.invalid_fastq_policy,
                    &mut stats,
                )?
                else {
                    continue;
                };

                arena.reset();

                let outcome = plan.execute(pair, &mut arena, &mut stats);
                write_paired_outcome(&mut output, &outcome, &mut stats)?;

                progress.maybe_report(&stats);
            }
            (None, None) => break,
            _ => bail!(
                "paired FASTQ inputs have different record counts\n\
                 source: {}\n\
                 complete_pairs_seen: {}\n\
                 reads_seen_before_failure: {}\n\
                 help: confirm --in1 and --in2 are mates from the same run and were not independently filtered or truncated",
                summary_context.input_label(),
                stats.pairs_seen,
                stats.reads_seen,
            ),
        }
    }

    progress.finish();
    output.finish().wrap_err(
        "failed to finalize paired output\nhelp: for gzip output this can indicate a truncated destination or broken downstream pipe",
    )?;
    let summary = report::RunSummary::from_stats(summary_context, &stats, started_at.elapsed());
    if ui.show_summary {
        report::print_summary(&summary);
    }
    if let Some(path) = &cli.summary {
        report::write_summary_json(path, &summary).wrap_err_with(|| {
            format!(
                "failed to write JSON run summary: {}\nhelp: check that the parent directory exists and is writable",
                path.display()
            )
        })?;
    }
    Ok(())
}

fn read_validated_single<'a, 'p, R>(
    parser: &'a FastqParser<'p, FASTQ_CONFIG, R>,
    summary_context: &'a RunContext,
    policy: crate::cli::InvalidFastqPolicy,
    stats: &mut ReadStats,
) -> Result<Option<RecordView<'a>>>
where
    R: InputData<'p>,
{
    let source = summary_context.input_label();
    let sequence = catch_parser_panic(&source, "single", stats, || parser.get_dna_string())?;
    let quality = catch_parser_panic(&source, "single", stats, || parser.get_quality())?
        .ok_or_else(|| missing_quality_error(&source, "single", stats))?;
    let header = catch_parser_panic(&source, "single", stats, || parser.get_header())?;
    let provenance = single_provenance(summary_context);
    let unvalidated_record = RecordView::new(header, sequence, quality).with_provenance(provenance);

    stats.record_seen(sequence.len());

    unvalidated_record.validate(policy, stats)
}

fn read_validated_pair<'a, 'p, R1, R2>(
    parser_r1: &'a FastqParser<'p, FASTQ_CONFIG, R1>,
    parser_r2: &'a FastqParser<'p, FASTQ_CONFIG, R2>,
    summary_context: &'a RunContext,
    policy: crate::cli::InvalidFastqPolicy,
    stats: &mut ReadStats,
) -> Result<Option<RecordPair<'a>>>
where
    R1: InputData<'p>,
    R2: InputData<'p>,
{
    let source = summary_context.input_label();
    let sequence_r1 = catch_parser_panic(&source, "left", stats, || parser_r1.get_dna_string())?;
    let sequence_r2 = catch_parser_panic(&source, "right", stats, || parser_r2.get_dna_string())?;
    let quality_r1 = catch_parser_panic(&source, "left", stats, || parser_r1.get_quality())?
        .ok_or_else(|| missing_quality_error(&source, "left", stats))?;
    let quality_r2 = catch_parser_panic(&source, "right", stats, || parser_r2.get_quality())?
        .ok_or_else(|| missing_quality_error(&source, "right", stats))?;
    let left_provenance = paired_provenance(summary_context, MateSide::Left);
    let right_provenance = paired_provenance(summary_context, MateSide::Right);
    let header_r1 = catch_parser_panic(&source, "left", stats, || parser_r1.get_header())?;
    let header_r2 = catch_parser_panic(&source, "right", stats, || parser_r2.get_header())?;
    let unvalidated_r1 =
        RecordView::new(header_r1, sequence_r1, quality_r1).with_provenance(left_provenance);
    let unvalidated_r2 =
        RecordView::new(header_r2, sequence_r2, quality_r2).with_provenance(right_provenance);

    stats.record_seen(sequence_r1.len());
    stats.record_seen(sequence_r2.len());
    stats.pairs_seen += 1;

    unvalidated_r1.validate_pair(unvalidated_r2, policy, stats)
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

fn write_paired_outcome(
    output: &mut PairedOutputHandle,
    outcome: &crate::plan::ExecutionOutcome<'_>,
    stats: &mut ReadStats,
) -> Result<()> {
    if outcome.is_fully_emitted() {
        output.write_pair_outcome(outcome, stats).wrap_err(
            "failed to write paired output record group\nhelp: check downstream pipe or output filesystem health",
        )?;
    } else if outcome.is_fully_rejected() && outcome.rejection_count() > 0 {
        stats.record_pair_rejected();
    } else if outcome.is_orphan() {
        let orphan_policy = OrphanPolicy::EmitOrphan;
        let rejection_code = outcome.first_rejection_code().unwrap_or("unknown");
        bail!(
            "paired preprocessing produced an orphan read, but the configured output cannot represent orphans yet\n\
             orphan_policy: {orphan_policy:?}\n\
             first_rejection: {rejection_code}\n\
             pairs_seen: {}\n\
             help: this is an internal plan/output mismatch; use DropPair for paired output until orphan output is implemented",
            stats.pairs_seen,
        );
    } else {
        bail!(
            "paired preprocessing produced an unexpected mixed outcome\n\
             pairs_seen: {} reads_seen: {} reads_rejected: {}\n\
             help: please report this with the input layout and preprocessing flags",
            stats.pairs_seen,
            stats.reads_seen,
            stats.reads_rejected,
        );
    }

    Ok(())
}

fn read_stats(cli: &Cli) -> Result<ReadStats> {
    let mut stats = ReadStats::default();
    if let Some(path) = &cli.invalid_fastq_report {
        stats.set_invalid_fastq_report(InvalidFastqReport::create(path).wrap_err_with(|| {
            format!(
                "failed to create invalid FASTQ JSONL report: {}\nhelp: check that the parent directory exists and is writable",
                path.display()
            )
        })?);
    }
    Ok(stats)
}

fn single_provenance(summary_context: &RunContext) -> RecordProvenance<'_> {
    match &summary_context.ingress_mode {
        report::IngressMode::Ena => RecordProvenance {
            source: InputSource::Ena {
                accession: summary_context
                    .accession
                    .as_deref()
                    .expect("ENA context should carry accession"),
            },
            mate: None,
        },
        report::IngressMode::Local => RecordProvenance {
            source: InputSource::LocalSingle {
                input: std::path::Path::new(
                    summary_context
                        .input1
                        .as_deref()
                        .expect("local single context should carry input path"),
                ),
            },
            mate: None,
        },
    }
}

fn paired_provenance(summary_context: &RunContext, mate: MateSide) -> RecordProvenance<'_> {
    match &summary_context.ingress_mode {
        report::IngressMode::Ena => RecordProvenance {
            source: InputSource::Ena {
                accession: summary_context
                    .accession
                    .as_deref()
                    .expect("ENA context should carry accession"),
            },
            mate: Some(mate),
        },
        report::IngressMode::Local => RecordProvenance {
            source: InputSource::LocalPaired {
                input1: std::path::Path::new(
                    summary_context
                        .input1
                        .as_deref()
                        .expect("paired local context should carry input1 path"),
                ),
                input2: std::path::Path::new(
                    summary_context
                        .input2
                        .as_deref()
                        .expect("paired local context should carry input2 path"),
                ),
            },
            mate: Some(mate),
        },
    }
}

#[cfg(test)]
mod tests {
    use std::{fs::File, io::Cursor, path::Path};

    use color_eyre::{Result, eyre::bail};
    use helicase::{Config, FastqParser, HelicaseParser, ParserOptions, input::FromReader};
    use tempfile::tempdir;

    use crate::{
        cli::{Cli, InvalidFastqPolicy},
        output::{
            InterleavedOutput, OutputArgs, OutputEncoding, OutputFormat, PairedRecordOutput,
            SingleOutput, SingleRecordOutput, StreamSink,
        },
        record::RecordView,
        report,
    };

    use super::run_paired;

    const TEST_FASTQ_CONFIG: Config = ParserOptions::default().compute_quality().config();

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
        let mut parser = FastqParser::<TEST_FASTQ_CONFIG, _>::from_reader(reader)?;
        while parser.next().is_some() {
            let record = RecordView::new(
                parser.get_header(),
                parser.get_dna_string(),
                parser
                    .get_quality()
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
        let mut parser = FastqParser::<TEST_FASTQ_CONFIG, _>::from_reader(reader)?;
        while parser.next().is_some() {
            let record = RecordView::new(
                parser.get_header(),
                parser.get_dna_string(),
                parser
                    .get_quality()
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
        let mut parser = FastqParser::<TEST_FASTQ_CONFIG, _>::from_reader(reader)?;
        while parser.next().is_some() {
            let record = RecordView::new(
                parser.get_header(),
                parser.get_dna_string(),
                parser
                    .get_quality()
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
        let mut parser_r1 = FastqParser::<TEST_FASTQ_CONFIG, _>::from_reader(r1)?;
        let mut parser_r2 = FastqParser::<TEST_FASTQ_CONFIG, _>::from_reader(r2)?;
        loop {
            match (parser_r1.next(), parser_r2.next()) {
                (Some(_), Some(_)) => output.write_pair(
                    RecordView::new(
                        parser_r1.get_header(),
                        parser_r1.get_dna_string(),
                        parser_r1
                            .get_quality()
                            .expect("FASTQ parser must provide quality scores"),
                    ),
                    RecordView::new(
                        parser_r2.get_header(),
                        parser_r2.get_dna_string(),
                        parser_r2
                            .get_quality()
                            .expect("FASTQ parser must provide quality scores"),
                    ),
                )?,
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
        let error = run_paired(
            r1,
            r2,
            output,
            &test_cli(),
            &test_cli().ui_policy(),
            report::RunContext {
                ingress_mode: report::IngressMode::Local,
                layout: report::RunLayout::Paired,
                accession: None,
                input1: Some("reads_1.fastq.gz".to_owned()),
                input2: Some("reads_2.fastq.gz".to_owned()),
            },
        )
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
