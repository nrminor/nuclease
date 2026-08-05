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
    error::{InternalError, IoError, ReportWriteError, Result, RunError},
    record::{InvalidFastqEvent, ReadStats},
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
    /// First invalid FASTQ events observed during this run.
    pub invalid_fastq_samples: Vec<InvalidFastqEvent>,
    /// Whether invalid FASTQ samples were truncated to the bounded summary limit.
    pub invalid_fastq_samples_truncated: bool,
    /// Rejection breakdown sorted by descending count.
    pub rejection_breakdown: Vec<CountBreakdown>,
    /// Transform breakdown sorted by descending count.
    pub transform_breakdown: Vec<CountBreakdown>,
}

impl RunSummary {
    /// Build an owned summary from run context, counters, and elapsed wall time.
    pub fn from_stats(context: RunContext, stats: &ReadStats, elapsed: Duration) -> Self {
        let elapsed_seconds = elapsed.as_secs_f64();

        Self {
            context,
            elapsed_seconds,
            reads_per_second: rate(stats.reads_seen, elapsed_seconds),
            bases_per_second: rate(stats.bases_seen, elapsed_seconds),
            reads_seen: stats.reads_seen,
            reads_emitted: stats.reads_emitted,
            reads_rejected: stats.reads_rejected,
            invalid_reads: stats.invalid_reads,
            read_retention_fraction: fraction(stats.reads_emitted, stats.reads_seen),
            read_rejection_fraction: fraction(stats.reads_rejected, stats.reads_seen),
            invalid_read_fraction: fraction(stats.invalid_reads, stats.reads_seen),
            bases_seen: stats.bases_seen,
            bases_emitted: stats.bases_emitted,
            base_retention_fraction: fraction(stats.bases_emitted, stats.bases_seen),
            pairs_seen: (stats.pairs_seen > 0).then_some(stats.pairs_seen),
            pairs_emitted: (stats.pairs_seen > 0).then_some(stats.pairs_emitted),
            pairs_rejected: (stats.pairs_seen > 0).then_some(stats.pairs_rejected),
            invalid_pairs: (stats.pairs_seen > 0).then_some(stats.invalid_pairs),
            pair_retention_fraction: opt_fraction(stats.pairs_emitted, stats.pairs_seen),
            pair_rejection_fraction: opt_fraction(stats.pairs_rejected, stats.pairs_seen),
            invalid_pair_fraction: opt_fraction(stats.invalid_pairs, stats.pairs_seen),
            invalid_fastq_samples: stats.invalid_fastq_samples.clone(),
            invalid_fastq_samples_truncated: stats.invalid_fastq_samples_truncated,
            rejection_breakdown: breakdowns(
                &stats.rejection_counts,
                stats.reads_seen,
                stats.reads_rejected,
            ),
            transform_breakdown: breakdowns(
                &stats.transform_counts,
                stats.reads_seen,
                stats.transform_counts.values().sum(),
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
        .map_err(|source| classify_json_error("run summary", path, source))?;
    writer.flush().map_err(|source| IoError::FinalizeReport {
        report_kind: "run summary",
        path: path.to_path_buf(),
        source,
    })?;
    Ok(())
}

fn classify_json_error(
    report_kind: &'static str,
    path: &Path,
    source: serde_json::Error,
) -> RunError {
    if source.io_error_kind().is_some() {
        IoError::WriteReport {
            report_kind,
            path: path.to_path_buf(),
            source: ReportWriteError::Json(source),
        }
        .into()
    } else {
        InternalError::SerializeReport {
            report_kind,
            source,
        }
        .into()
    }
}

fn render_summary(summary: &RunSummary) -> String {
    let mut output = String::new();
    output.push('\n');
    output.push_str("nuclease summary\n");
    write_context(&mut output, summary);
    write_totals(&mut output, summary);
    write_pairs(&mut output, summary);
    write_invalid_fastq_samples(&mut output, summary);
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

fn write_invalid_fastq_samples(output: &mut String, summary: &RunSummary) {
    if summary.invalid_fastq_samples.is_empty() {
        return;
    }

    output.push_str("\n  invalid FASTQ samples:\n");
    for event in &summary.invalid_fastq_samples {
        match event.kind {
            "sequence_quality_length_mismatch" => {
                let _ = writeln!(
                    output,
                    "    {} source={} mate={} header={} sequence_len={} quality_len={} reads_seen={} pairs_seen={}",
                    event.kind,
                    event.source,
                    event.mate.unwrap_or("unknown"),
                    event.header.as_deref().unwrap_or("<unknown>"),
                    event.sequence_len.unwrap_or_default(),
                    event.quality_len.unwrap_or_default(),
                    event.reads_seen,
                    event
                        .pairs_seen
                        .map_or_else(|| "n/a".to_owned(), |pairs_seen| pairs_seen.to_string()),
                );
            }
            "paired_header_mismatch" => {
                let _ = writeln!(
                    output,
                    "    {} source={} left_header={} right_header={} reads_seen={} pairs_seen={}",
                    event.kind,
                    event.source,
                    event.left_header.as_deref().unwrap_or("<unknown>"),
                    event.right_header.as_deref().unwrap_or("<unknown>"),
                    event.reads_seen,
                    event
                        .pairs_seen
                        .map_or_else(|| "n/a".to_owned(), |pairs_seen| pairs_seen.to_string()),
                );
            }
            _ => {
                let _ = writeln!(
                    output,
                    "    {} source={} reads_seen={}",
                    event.kind, event.source, event.reads_seen,
                );
            }
        }
    }

    if summary.invalid_fastq_samples_truncated {
        output.push_str("    ... additional invalid FASTQ events omitted from summary\n");
    }
}

fn write_breakdowns(output: &mut String, summary: &RunSummary) {
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

fn opt_fraction(numerator: u64, denominator: u64) -> Option<f64> {
    (denominator > 0).then(|| fraction(numerator, denominator))
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
        collections::BTreeMap,
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
        record::{InvalidFastqEvent, ReadStats},
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

    fn sample_stats() -> ReadStats {
        let mut rejection_counts = BTreeMap::new();
        rejection_counts.insert("too_short", 3);
        rejection_counts.insert("too_many_ns", 1);

        let mut transform_counts = BTreeMap::new();
        transform_counts.insert("trim_adapters", 4);

        ReadStats {
            reads_seen: 10,
            reads_emitted: 6,
            reads_rejected: 4,
            invalid_reads: 1,
            bases_seen: 100,
            bases_emitted: 72,
            pairs_seen: 5,
            pairs_emitted: 3,
            pairs_rejected: 2,
            invalid_pairs: 1,
            rejection_counts,
            transform_counts,
            invalid_fastq_warnings_emitted: 0,
            invalid_fastq_warnings_suppressed: false,
            invalid_fastq_samples: vec![InvalidFastqEvent {
                kind: "sequence_quality_length_mismatch",
                source: "ena:SRR35939766".to_owned(),
                mate: Some("right"),
                header: Some("SRR35939766.42 instrument/2".to_owned()),
                sequence_len: Some(267),
                quality_len: Some(20),
                left_mate: None,
                right_mate: None,
                left_header: None,
                right_header: None,
                reads_seen: 84,
                pairs_seen: Some(42),
                policy: "warn_drop".to_owned(),
                recoverable: true,
                fatal: false,
                parser_error_kind: None,
                parser_error_message: None,
                parser_error_line: None,
            }],
            invalid_fastq_samples_truncated: false,
            invalid_fastq_report: None,
        }
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
            RunSummary::from_stats(sample_context(), &sample_stats(), Duration::from_secs(2));

        assert_eq!(summary.rejection_breakdown[0].code, "too_short");
        assert_eq!(summary.rejection_breakdown[0].count, 3);
        assert_eq!(summary.transform_breakdown[0].code, "trim_adapters");
        assert!((summary.read_retention_fraction - 0.6).abs() < f64::EPSILON);
    }

    #[test]
    fn render_summary_includes_breakdowns_and_metadata() {
        let summary =
            RunSummary::from_stats(sample_context(), &sample_stats(), Duration::from_secs(2));
        let rendered = render_summary(&summary);

        assert!(rendered.contains("ingress mode:        ena"));
        assert!(rendered.contains("layout:              paired"));
        assert!(rendered.contains("rejection reasons:"));
        assert!(rendered.contains("too_short"));
        assert!(rendered.contains("transforms applied:"));
        assert!(rendered.contains("trim_adapters"));
        assert!(rendered.contains("invalid FASTQ samples:"));
        assert!(rendered.contains("SRR35939766.42 instrument/2"));
    }

    #[test]
    fn write_summary_json_serializes_summary() {
        let temp = tempdir().expect("tempdir should exist");
        let path = temp.path().join("summary.json");
        let summary =
            RunSummary::from_stats(sample_context(), &sample_stats(), Duration::from_secs(2));

        write_summary_json(&path, &summary).expect("json summary should write");

        let written = std::fs::read_to_string(path).expect("json summary should be readable");
        assert!(written.contains("\"accession\": \"SRR35939766\""));
        assert!(written.contains("\"rejection_breakdown\""));
        assert!(written.contains("\"transform_breakdown\""));
        assert!(written.contains("\"invalid_fastq_samples\""));
        assert!(written.contains("SRR35939766.42 instrument/2"));
    }

    #[test]
    fn summary_write_failure_retains_path_and_json_io_source() {
        let summary =
            RunSummary::from_stats(sample_context(), &sample_stats(), Duration::from_secs(2));
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
            RunSummary::from_stats(sample_context(), &sample_stats(), Duration::from_secs(2));
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
            RunSummary::from_stats(sample_context(), &sample_stats(), Duration::from_secs(2));

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
