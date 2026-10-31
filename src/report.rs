//! Human-readable run summaries written after a streaming job completes.

use std::{
    fmt::Write as _,
    fs::File,
    io::{self, BufWriter, Write as _},
    path::Path,
    time::Duration,
};

use serde::Serialize;

use crate::{
    error::{self, IoError, Result},
    observer::{InvalidInputEvent, RunObserver},
};

/// Lightweight run context needed to explain a preprocessing run.
#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct RunContext {
    /// Whether reads arrived from ENA or local user-provided FASTQ.
    pub ingress_mode: IngressMode,
    /// Whether reads were processed in single-end or paired-end layout.
    pub layout: RunLayout,
    /// ENA run accession when ENA mode was used.
    pub accession: Option<String>,
    /// First local FASTQ path when local ingress was used.
    pub input1: Option<String>,
    /// Second local FASTQ path when paired local ingress was used.
    pub input2: Option<String>,
}

/// Ingress origin for a run summary.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum IngressMode {
    /// Reads were streamed from ENA.
    Ena,
    /// Reads were read from local FASTQ files.
    Local,
}

/// Record layout processed by the run.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum RunLayout {
    /// Single-end records.
    Single,
    /// Paired-end records.
    Paired,
}

/// Stable count breakdown for one rejection reason or transform.
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct CountBreakdown {
    /// Stable machine-readable code.
    pub code: String,
    /// Number of times this code was observed.
    pub count: u64,
    /// Fraction relative to all reads seen.
    pub fraction_of_reads_seen: f64,
    /// Fraction relative to the code family total when applicable.
    pub fraction_of_category: f64,
}

/// Compact owned summary for stderr presentation and JSON serialization.
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct RunSummary {
    /// Run context explaining where data came from and how it was processed.
    pub context: RunContext,
    /// Total elapsed runtime in seconds.
    pub elapsed_seconds: f64,
    /// Overall read throughput.
    pub reads_per_second: f64,
    /// Overall base throughput.
    pub bases_per_second: f64,
    /// Total reads observed.
    pub reads_seen: u64,
    /// Total reads emitted.
    pub reads_emitted: u64,
    /// Total reads rejected.
    pub reads_rejected: u64,
    /// Total reads dropped at ingress because the FASTQ record was malformed.
    pub invalid_reads: u64,
    /// Fraction of reads retained.
    pub read_retention_fraction: f64,
    /// Fraction of reads rejected.
    pub read_rejection_fraction: f64,
    /// Fraction of reads dropped at ingress for FASTQ invalidity.
    pub invalid_read_fraction: f64,
    /// Total bases observed.
    pub bases_seen: u64,
    /// Total bases emitted.
    pub bases_emitted: u64,
    /// Fraction of bases retained.
    pub base_retention_fraction: f64,
    /// Total pairs observed when paired input was processed.
    pub pairs_seen: Option<u64>,
    /// Total pairs emitted when paired input was processed.
    pub pairs_emitted: Option<u64>,
    /// Total pairs rejected when paired input was processed.
    pub pairs_rejected: Option<u64>,
    /// Total pairs dropped at ingress because one or both mates were malformed.
    pub invalid_pairs: Option<u64>,
    /// Fraction of pairs retained when paired input was processed.
    pub pair_retention_fraction: Option<f64>,
    /// Fraction of pairs rejected when paired input was processed.
    pub pair_rejection_fraction: Option<f64>,
    /// Fraction of pairs dropped at ingress for FASTQ invalidity.
    pub invalid_pair_fraction: Option<f64>,
    /// First invalid-input events observed during this run.
    #[serde(rename = "admission_samples")]
    pub invalid_input_samples: Vec<InvalidInputEvent>,
    /// Whether invalid-input samples were truncated to the bounded summary limit.
    #[serde(rename = "admission_samples_truncated")]
    pub invalid_input_samples_truncated: bool,
    /// Invalid-input breakdown sorted by descending count.
    #[serde(rename = "admission_event_breakdown")]
    pub invalid_input_event_breakdown: Vec<CountBreakdown>,
    /// Rejection breakdown sorted by descending count.
    pub rejection_breakdown: Vec<CountBreakdown>,
    /// Transform breakdown sorted by descending count.
    pub transform_breakdown: Vec<CountBreakdown>,
}

impl RunSummary {
    /// Build an owned summary from run context, counters, and elapsed wall time.
    pub fn from_observer(context: RunContext, observer: &RunObserver, elapsed: Duration) -> Self {
        let elapsed_seconds = elapsed.as_secs_f64();
        let paired = context.layout == RunLayout::Paired;

        Self {
            context,
            elapsed_seconds,
            reads_per_second: rate(observer.reads_seen, elapsed_seconds),
            bases_per_second: rate(observer.bases_seen, elapsed_seconds),
            reads_seen: observer.reads_seen,
            reads_emitted: observer.reads_emitted,
            reads_rejected: observer.reads_rejected,
            invalid_reads: observer.invalid_reads,
            read_retention_fraction: fraction(observer.reads_emitted, observer.reads_seen),
            read_rejection_fraction: fraction(observer.reads_rejected, observer.reads_seen),
            invalid_read_fraction: fraction(observer.invalid_reads, observer.reads_seen),
            bases_seen: observer.bases_seen,
            bases_emitted: observer.bases_emitted,
            base_retention_fraction: fraction(observer.bases_emitted, observer.bases_seen),
            pairs_seen: paired.then_some(observer.pairs_seen),
            pairs_emitted: paired.then_some(observer.pairs_emitted),
            pairs_rejected: paired.then_some(observer.pairs_rejected),
            invalid_pairs: paired.then_some(observer.invalid_pairs),
            pair_retention_fraction: paired
                .then(|| fraction(observer.pairs_emitted, observer.pairs_seen)),
            pair_rejection_fraction: paired
                .then(|| fraction(observer.pairs_rejected, observer.pairs_seen)),
            invalid_pair_fraction: paired
                .then(|| fraction(observer.invalid_pairs, observer.pairs_seen)),
            invalid_input_samples: observer.invalid_input_samples().to_vec(),
            invalid_input_samples_truncated: observer.invalid_input_samples_truncated(),
            invalid_input_event_breakdown: breakdowns(
                observer.invalid_input_event_counts(),
                observer.reads_seen,
                observer.invalid_input_event_counts().values().sum(),
            ),
            rejection_breakdown: breakdowns(
                &observer.rejection_counts,
                observer.reads_seen,
                observer.reads_rejected,
            ),
            transform_breakdown: breakdowns(
                &observer.transform_counts,
                observer.reads_seen,
                observer.transform_counts.values().sum(),
            ),
        }
    }
}

/// Print a human-readable run summary to stderr.
pub fn print_summary(summary: &RunSummary) {
    let rendered = render_summary(summary);
    let _ = io::stderr().write_all(rendered.as_bytes());
}

/// Write a JSON summary to disk.
///
/// # Errors
///
/// Returns an error when the file cannot be created or the JSON cannot be serialized.
pub fn write_summary_json(path: &Path, summary: &RunSummary) -> Result<()> {
    let file = File::create(path).map_err(|source| IoError::CreateReport {
        report_kind: "run summary",
        path: path.to_path_buf(),
        source,
    })?;
    let mut writer = BufWriter::new(file);
    write_summary_to(&mut writer, path, summary)
}

fn write_summary_to(
    writer: &mut impl std::io::Write,
    path: &Path,
    summary: &RunSummary,
) -> Result<()> {
    serde_json::to_writer_pretty(&mut *writer, summary)
        .map_err(|source| error::constructors::json_report_error("run summary", path, source))?;
    writer.flush().map_err(|source| IoError::FinalizeReport {
        report_kind: "run summary",
        path: path.to_path_buf(),
        source,
    })?;
    Ok(())
}

fn render_summary(summary: &RunSummary) -> String {
    let mut output = String::new();
    output.push('\n');
    output.push_str("nuclease summary\n");
    write_context(&mut output, summary);
    write_totals(&mut output, summary);
    write_pairs(&mut output, summary);
    write_invalid_input_samples(&mut output, summary);
    write_breakdowns(&mut output, summary);
    output
}

fn write_context(output: &mut String, summary: &RunSummary) {
    let _ = writeln!(
        output,
        "  ingress mode:        {}",
        match summary.context.ingress_mode {
            IngressMode::Ena => "ena",
            IngressMode::Local => "local",
        }
    );
    let _ = writeln!(
        output,
        "  layout:              {}",
        match summary.context.layout {
            RunLayout::Single => "single",
            RunLayout::Paired => "paired",
        }
    );

    if let Some(accession) = &summary.context.accession {
        let _ = writeln!(output, "  accession:           {accession}");
    }

    if let Some(input1) = &summary.context.input1 {
        let _ = writeln!(output, "  input 1:             {input1}");
    }

    if let Some(input2) = &summary.context.input2 {
        let _ = writeln!(output, "  input 2:             {input2}");
    }
}

fn write_totals(output: &mut String, summary: &RunSummary) {
    let _ = writeln!(
        output,
        "  elapsed:             {:.2}s",
        summary.elapsed_seconds
    );
    let _ = writeln!(
        output,
        "  throughput:          {:.1} reads/s, {:.1} bases/s",
        summary.reads_per_second, summary.bases_per_second
    );
    let _ = writeln!(output, "  reads seen:          {}", summary.reads_seen);
    let _ = writeln!(
        output,
        "  reads emitted:       {} ({:.2}%)",
        summary.reads_emitted,
        summary.read_retention_fraction * 100.0
    );
    let _ = writeln!(
        output,
        "  reads rejected:      {} ({:.2}%)",
        summary.reads_rejected,
        summary.read_rejection_fraction * 100.0
    );
    let _ = writeln!(
        output,
        "  invalid reads:       {} ({:.2}%)",
        summary.invalid_reads,
        summary.invalid_read_fraction * 100.0
    );
    let _ = writeln!(output, "  bases seen:          {}", summary.bases_seen);
    let _ = writeln!(
        output,
        "  bases emitted:       {} ({:.2}%)",
        summary.bases_emitted,
        summary.base_retention_fraction * 100.0
    );
}

fn write_pairs(output: &mut String, summary: &RunSummary) {
    if let (Some(pairs_seen), Some(pairs_emitted), Some(pairs_rejected)) = (
        summary.pairs_seen,
        summary.pairs_emitted,
        summary.pairs_rejected,
    ) {
        let _ = writeln!(output, "  pairs seen:          {pairs_seen}");
        let _ = writeln!(
            output,
            "  pairs emitted:       {} ({:.2}%)",
            pairs_emitted,
            summary.pair_retention_fraction.unwrap_or_default() * 100.0
        );
        let _ = writeln!(
            output,
            "  pairs rejected:      {} ({:.2}%)",
            pairs_rejected,
            summary.pair_rejection_fraction.unwrap_or_default() * 100.0
        );
        let _ = writeln!(
            output,
            "  invalid pairs:       {} ({:.2}%)",
            summary.invalid_pairs.unwrap_or_default(),
            summary.invalid_pair_fraction.unwrap_or_default() * 100.0
        );
    }
}

fn write_invalid_input_samples(output: &mut String, summary: &RunSummary) {
    if summary.invalid_input_samples.is_empty() {
        return;
    }

    output.push_str("\n  invalid-input samples:\n");
    for event in &summary.invalid_input_samples {
        match event {
            InvalidInputEvent::SequenceQualityLengthMismatch {
                source,
                mate,
                header,
                sequence_len,
                quality_len,
                reads_seen,
                pairs_seen,
                continued,
            } => {
                let mate = mate.map_or_else(|| "single".to_owned(), |mate| mate.to_string());
                let _ = writeln!(
                    output,
                    "    sequence_quality_length_mismatch source={source} mate={mate} header={header} sequence_len={sequence_len} quality_len={quality_len} reads_seen={reads_seen} pairs_seen={} continued={continued}",
                    pairs_seen
                        .map_or_else(|| "n/a".to_owned(), |pairs_seen| pairs_seen.to_string()),
                );
            }
            InvalidInputEvent::PairConstructionFailure {
                source,
                error,
                reads_seen,
                pairs_seen,
                continued,
            } => {
                let _ = writeln!(
                    output,
                    "    {} source={source} error={error} reads_seen={reads_seen} pairs_seen={pairs_seen} continued={continued}",
                    error.code(),
                );
            }
            InvalidInputEvent::MissingMate {
                source,
                present_mate,
                header,
                reads_seen,
                pairs_seen,
                continued,
            } => {
                let _ = writeln!(
                    output,
                    "    missing_mate source={source} present_mate={present_mate} header={header} reads_seen={reads_seen} pairs_seen={pairs_seen} continued={continued}",
                );
            }
            InvalidInputEvent::SingleRecordParseFailure {
                source,
                failure,
                reads_seen,
                continued,
            } => {
                let _ = writeln!(
                    output,
                    "    record_parse_failure source={source} parser_kind={} message={} line={} reads_seen={reads_seen} continued={continued}",
                    failure.parser_kind,
                    failure.message,
                    failure
                        .line
                        .map_or_else(|| "n/a".to_owned(), |line| line.to_string()),
                );
            }
            InvalidInputEvent::PairedRecordParseFailure {
                source,
                mate,
                failure,
                reads_seen,
                pairs_seen,
                continued,
            } => {
                let _ = writeln!(
                    output,
                    "    record_parse_failure source={source} mate={mate} parser_kind={} message={} line={} reads_seen={reads_seen} pairs_seen={pairs_seen} continued={continued}",
                    failure.parser_kind,
                    failure.message,
                    failure
                        .line
                        .map_or_else(|| "n/a".to_owned(), |line| line.to_string()),
                );
            }
        }
    }

    if summary.invalid_input_samples_truncated {
        output.push_str("    ... additional invalid-input events omitted from summary\n");
    }
}

fn write_breakdowns(output: &mut String, summary: &RunSummary) {
    if !summary.invalid_input_event_breakdown.is_empty() {
        output.push_str("\n  invalid-input events:\n");
        for breakdown in &summary.invalid_input_event_breakdown {
            let _ = writeln!(
                output,
                "    {:<32} {:>10} ({:.2}% of reads, {:.2}% of invalid-input events)",
                breakdown.code,
                breakdown.count,
                breakdown.fraction_of_reads_seen * 100.0,
                breakdown.fraction_of_category * 100.0
            );
        }
    }

    if !summary.rejection_breakdown.is_empty() {
        output.push_str("\n  rejection reasons:\n");
        for breakdown in &summary.rejection_breakdown {
            let _ = writeln!(
                output,
                "    {:<20} {:>10} ({:.2}% of reads, {:.2}% of rejected)",
                breakdown.code,
                breakdown.count,
                breakdown.fraction_of_reads_seen * 100.0,
                breakdown.fraction_of_category * 100.0
            );
        }
    }

    if !summary.transform_breakdown.is_empty() {
        output.push_str("\n  transforms applied:\n");
        for breakdown in &summary.transform_breakdown {
            let _ = writeln!(
                output,
                "    {:<20} {:>10} ({:.2}% of reads)",
                breakdown.code,
                breakdown.count,
                breakdown.fraction_of_reads_seen * 100.0,
            );
        }
    }
}

fn breakdowns(
    counts: &std::collections::BTreeMap<&'static str, u64>,
    reads_seen: u64,
    category_total: u64,
) -> Vec<CountBreakdown> {
    let mut breakdowns = counts
        .iter()
        .map(|(code, count)| CountBreakdown {
            code: (*code).to_owned(),
            count: *count,
            fraction_of_reads_seen: fraction(*count, reads_seen),
            fraction_of_category: fraction(*count, category_total),
        })
        .collect::<Vec<_>>();

    breakdowns.sort_by(|left, right| {
        right
            .count
            .cmp(&left.count)
            .then_with(|| left.code.cmp(&right.code))
    });
    breakdowns
}

fn fraction(numerator: u64, denominator: u64) -> f64 {
    if denominator == 0 {
        0.0
    } else {
        u64_to_f64(numerator) / u64_to_f64(denominator)
    }
}

fn rate(total: u64, elapsed_seconds: f64) -> f64 {
    if elapsed_seconds <= f64::EPSILON {
        0.0
    } else {
        u64_to_f64(total) / elapsed_seconds
    }
}

fn u64_to_f64(value: u64) -> f64 {
    value
        .to_string()
        .parse::<f64>()
        .expect("u64 should always parse into f64")
}

#[cfg(test)]
mod tests {
    use std::{
        io::{self, Write},
        time::Duration,
    };

    use tempfile::tempdir;

    use super::{
        IngressMode, RunContext, RunLayout, RunSummary, render_summary, write_summary_json,
        write_summary_to,
    };
    use crate::{
        error::{IoError, ReportWriteError, RunError},
        observer::{InvalidInputEvent, RunObserver},
        record::MateSide,
    };

    struct FailingWriter {
        fail_write: bool,
        fail_flush: bool,
    }

    impl Write for FailingWriter {
        fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
            if self.fail_write {
                Err(io::Error::new(
                    io::ErrorKind::StorageFull,
                    "test write failure",
                ))
            } else {
                Ok(bytes.len())
            }
        }

        fn flush(&mut self) -> io::Result<()> {
            if self.fail_flush {
                Err(io::Error::new(
                    io::ErrorKind::StorageFull,
                    "test flush failure",
                ))
            } else {
                Ok(())
            }
        }
    }

    fn sample_observer() -> RunObserver {
        let mut observer = RunObserver::new("ena:SRR35939766".to_owned());
        observer.reads_seen = 10;
        observer.reads_emitted = 6;
        observer.reads_rejected = 4;
        observer.invalid_reads = 1;
        observer.bases_seen = 100;
        observer.bases_emitted = 72;
        observer.pairs_seen = 5;
        observer.pairs_emitted = 3;
        observer.pairs_rejected = 2;
        observer.invalid_pairs = 1;
        observer.rejection_counts.insert("too_short", 3);
        observer.rejection_counts.insert("too_many_ns", 1);
        observer.transform_counts.insert("trim_adapters", 4);
        observer
            .record_invalid_input(InvalidInputEvent::SequenceQualityLengthMismatch {
                source: "ena:SRR35939766".to_owned(),
                mate: Some(MateSide::Right),
                header: "SRR35939766.42 instrument/2".to_owned(),
                sequence_len: 267,
                quality_len: 20,
                reads_seen: 84,
                pairs_seen: Some(42),
                continued: true,
            })
            .expect("sample invalid-input event should be recorded");
        observer
    }

    fn sample_context() -> RunContext {
        RunContext {
            ingress_mode: IngressMode::Ena,
            layout: RunLayout::Paired,
            accession: Some("SRR35939766".to_owned()),
            input1: None,
            input2: None,
        }
    }

    #[test]
    fn run_summary_sorts_breakdowns_by_count() {
        let summary =
            RunSummary::from_observer(sample_context(), &sample_observer(), Duration::from_secs(2));

        assert_eq!(summary.rejection_breakdown[0].code, "too_short");
        assert_eq!(summary.rejection_breakdown[0].count, 3);
        assert_eq!(summary.transform_breakdown[0].code, "trim_adapters");
        assert!((summary.read_retention_fraction - 0.6).abs() < f64::EPSILON);
    }

    #[test]
    fn render_summary_includes_breakdowns_and_metadata() {
        let summary =
            RunSummary::from_observer(sample_context(), &sample_observer(), Duration::from_secs(2));
        let rendered = render_summary(&summary);

        assert!(rendered.contains("ingress mode:        ena"));
        assert!(rendered.contains("layout:              paired"));
        assert!(rendered.contains("rejection reasons:"));
        assert!(rendered.contains("too_short"));
        assert!(rendered.contains("transforms applied:"));
        assert!(rendered.contains("trim_adapters"));
        assert!(rendered.contains("invalid-input samples:"));
        assert!(rendered.contains("SRR35939766.42 instrument/2"));
    }

    #[test]
    fn write_summary_json_serializes_summary() {
        let temp = tempdir().expect("tempdir should exist");
        let path = temp.path().join("summary.json");
        let summary =
            RunSummary::from_observer(sample_context(), &sample_observer(), Duration::from_secs(2));

        write_summary_json(&path, &summary).expect("json summary should write");

        let written = std::fs::read_to_string(path).expect("json summary should be readable");
        assert!(written.contains("\"accession\": \"SRR35939766\""));
        assert!(written.contains("\"rejection_breakdown\""));
        assert!(written.contains("\"transform_breakdown\""));
        assert!(written.contains("\"admission_samples\""));
        assert!(written.contains("SRR35939766.42 instrument/2"));
    }

    #[test]
    fn summary_write_failure_retains_path_and_json_io_source() {
        let summary =
            RunSummary::from_observer(sample_context(), &sample_observer(), Duration::from_secs(2));
        let path = std::path::Path::new("summary.json");
        let mut writer = FailingWriter {
            fail_write: true,
            fail_flush: false,
        };

        let error = write_summary_to(&mut writer, path, &summary)
            .expect_err("summary writer failure should be required I/O");
        assert!(matches!(
            error,
            RunError::Io(IoError::WriteReport {
                path: observed_path,
                source: ReportWriteError::Json(source),
                ..
            }) if observed_path == path && source.io_error_kind() == Some(io::ErrorKind::StorageFull)
        ));
    }

    #[test]
    fn summary_flush_failure_is_report_finalization() {
        let summary =
            RunSummary::from_observer(sample_context(), &sample_observer(), Duration::from_secs(2));
        let path = std::path::Path::new("summary.json");
        let mut writer = FailingWriter {
            fail_write: false,
            fail_flush: true,
        };

        let error = write_summary_to(&mut writer, path, &summary)
            .expect_err("summary flush failure should be required I/O");
        assert!(matches!(
            error,
            RunError::Io(IoError::FinalizeReport {
                path: observed_path,
                source,
                ..
            }) if observed_path == path && source.kind() == io::ErrorKind::StorageFull
        ));
    }

    #[test]
    fn summary_create_failure_retains_requested_path() {
        let temp = tempdir().expect("tempdir should exist");
        let path = temp.path().join("missing-parent").join("summary.json");
        let summary =
            RunSummary::from_observer(sample_context(), &sample_observer(), Duration::from_secs(2));

        let error = write_summary_json(&path, &summary)
            .expect_err("summary under missing parent should fail to open");
        assert!(matches!(
            error,
            RunError::Io(IoError::CreateReport {
                path: observed_path,
                source,
                ..
            }) if observed_path == path && source.kind() == io::ErrorKind::NotFound
        ));
    }
}
