//! Core record and accounting types shared across ingress, parsing, and output layers.

use std::{
    collections::BTreeMap,
    fs::File,
    io::BufWriter,
    path::{Path, PathBuf},
};

use serde::Serialize;
use tracing::warn;

use crate::{
    cli::InvalidFastqPolicy,
    error::{
        EnaContentProblem, IndeterminateInputError, InternalError, IoError, MalformedInputError,
        ReportWriteError, Result, RunError,
    },
    plan::RecordPair,
};

const INVALID_FASTQ_SAMPLE_LIMIT: usize = 20;
const INVALID_FASTQ_WARNING_LIMIT: u64 = 5;

/// Bounded, owned trace information about an invalid FASTQ event.
#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct InvalidFastqEvent {
    /// Stable event kind for machine-readable summaries.
    pub kind: &'static str,
    /// Source label such as `ena:SRR...` or a local path label.
    pub source: String,
    /// Mate label for record-level events.
    pub mate: Option<&'static str>,
    /// Header for record-level events.
    pub header: Option<String>,
    /// Sequence length for record-level events.
    pub sequence_len: Option<usize>,
    /// Quality length for record-level events.
    pub quality_len: Option<usize>,
    /// Left mate label for paired events.
    pub left_mate: Option<&'static str>,
    /// Right mate label for paired events.
    pub right_mate: Option<&'static str>,
    /// Left header for paired events.
    pub left_header: Option<String>,
    /// Right header for paired events.
    pub right_header: Option<String>,
    /// Total reads observed when the event was detected.
    pub reads_seen: u64,
    /// Total pairs observed when the event was detected, if paired input is active.
    pub pairs_seen: Option<u64>,
    /// Invalid FASTQ handling policy active when the event was observed.
    pub policy: String,
    /// Whether the event occurred at a known safe record or pair boundary.
    pub recoverable: bool,
    /// Whether the event forces this run to stop.
    pub fatal: bool,
    /// Parser-specific error kind for fatal parser-level FASTQ failures.
    pub parser_error_kind: Option<String>,
    /// Parser-specific error message for fatal parser-level FASTQ failures.
    pub parser_error_message: Option<String>,
    /// Parser-reported line number for fatal parser-level FASTQ failures.
    pub parser_error_line: Option<u64>,
}

/// Newline-delimited JSON writer for invalid FASTQ events.
#[derive(Debug)]
pub struct InvalidFastqReport {
    path: PathBuf,
    writer: BufWriter<File>,
}

impl InvalidFastqReport {
    /// Create a JSONL report at `path`.
    ///
    /// # Errors
    ///
    /// Returns an error when the destination file cannot be created.
    pub fn create(path: &Path) -> Result<Self> {
        let writer = File::create(path).map_err(|source| IoError::CreateReport {
            report_kind: "invalid FASTQ JSONL",
            path: path.to_path_buf(),
            source,
        })?;
        Ok(Self {
            path: path.to_path_buf(),
            writer: BufWriter::new(writer),
        })
    }

    fn write_event(&mut self, event: &InvalidFastqEvent) -> Result<()> {
        write_invalid_fastq_event_to(&mut self.writer, &self.path, event)
    }

    fn finish(&mut self) -> Result<()> {
        finalize_invalid_fastq_report(&mut self.writer, &self.path)
    }
}

fn write_invalid_fastq_event_to(
    writer: &mut impl std::io::Write,
    path: &Path,
    event: &InvalidFastqEvent,
) -> Result<()> {
    serde_json::to_writer(&mut *writer, event)
        .map_err(|source| classify_json_error(path, source))?;
    writer
        .write_all(b"\n")
        .map_err(|source| IoError::WriteReport {
            report_kind: "invalid FASTQ JSONL",
            path: path.to_path_buf(),
            source: ReportWriteError::Bytes(source),
        })?;
    if event.fatal {
        // Fatal parser events explain why the run is about to exit. Flush those diagnostics
        // before returning the error so failed-job artifacts do not depend on normal teardown.
        finalize_invalid_fastq_report(writer, path)?;
    }
    Ok(())
}

fn finalize_invalid_fastq_report(writer: &mut impl std::io::Write, path: &Path) -> Result<()> {
    writer.flush().map_err(|source| {
        IoError::FinalizeReport {
            report_kind: "invalid FASTQ JSONL",
            path: path.to_path_buf(),
            source,
        }
        .into()
    })
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct InvalidFastqContext {
    reads_seen: u64,
    pairs_seen: Option<u64>,
    policy: InvalidFastqPolicy,
}

impl InvalidFastqContext {
    fn length_mismatch(self, record: RecordView<'_>) -> InvalidFastqEvent {
        InvalidFastqEvent {
            kind: "sequence_quality_length_mismatch",
            source: record.source_display(),
            mate: Some(record.mate_display()),
            header: Some(String::from_utf8_lossy(record.header()).into_owned()),
            sequence_len: Some(record.sequence.len()),
            quality_len: Some(record.quality.len()),
            left_mate: None,
            right_mate: None,
            left_header: None,
            right_header: None,
            reads_seen: self.reads_seen,
            pairs_seen: self.pairs_seen,
            policy: self.policy.to_string(),
            recoverable: true,
            fatal: self.policy == InvalidFastqPolicy::Error,
            parser_error_kind: None,
            parser_error_message: None,
            parser_error_line: None,
        }
    }

    fn paired_header_mismatch(
        self,
        left: RecordView<'_>,
        right: RecordView<'_>,
    ) -> InvalidFastqEvent {
        InvalidFastqEvent {
            kind: "paired_header_mismatch",
            source: left.source_display(),
            mate: None,
            header: None,
            sequence_len: None,
            quality_len: None,
            left_mate: Some(left.mate_display()),
            right_mate: Some(right.mate_display()),
            left_header: Some(String::from_utf8_lossy(left.header()).into_owned()),
            right_header: Some(String::from_utf8_lossy(right.header()).into_owned()),
            reads_seen: self.reads_seen,
            pairs_seen: self.pairs_seen,
            policy: self.policy.to_string(),
            recoverable: true,
            fatal: self.policy == InvalidFastqPolicy::Error,
            parser_error_kind: None,
            parser_error_message: None,
            parser_error_line: None,
        }
    }

    pub(crate) fn parse_error(
        self,
        source: &str,
        mate: &'static str,
        parser_error_kind: String,
        parser_error_message: String,
        parser_error_line: Option<u64>,
    ) -> InvalidFastqEvent {
        InvalidFastqEvent {
            kind: "fastq_parse_error",
            source: source.to_owned(),
            mate: Some(mate),
            header: None,
            sequence_len: None,
            quality_len: None,
            left_mate: None,
            right_mate: None,
            left_header: None,
            right_header: None,
            reads_seen: self.reads_seen,
            pairs_seen: self.pairs_seen,
            policy: self.policy.to_string(),
            recoverable: false,
            fatal: true,
            parser_error_kind: Some(parser_error_kind),
            parser_error_message: Some(parser_error_message),
            parser_error_line,
        }
    }
}

/// Borrowed view of a biological sequence record.
///
/// Implementors expose header, sequence, and optional quality bytes without forcing ownership.
pub trait SequenceRecordRef {
    /// Return the raw record header bytes without the leading `@` or `>` marker.
    fn header(&self) -> &[u8];

    /// Return the raw sequence bytes for the record.
    fn sequence(&self) -> &[u8];

    /// Return quality bytes when the record format carries them.
    fn quality(&self) -> Option<&[u8]>;
}

/// Provenance attached to one borrowed record without taking ownership.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct RecordProvenance<'a> {
    /// Upstream data source that yielded this record.
    pub source: InputSource<'a>,
    /// Mate identity when the record originated from paired-end ingress.
    pub mate: Option<MateSide>,
}

/// Upstream source information attached to one borrowed record.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum InputSource<'a> {
    /// Record streamed directly from ENA.
    Ena { accession: &'a str },
    /// Record loaded from one local FASTQ file.
    LocalSingle { input: &'a Path },
    /// Record loaded from one interleaved local paired FASTQ file.
    LocalInterleavedPaired { input: &'a Path },
    /// Record loaded from one of two local paired FASTQ files.
    LocalPaired { input1: &'a Path, input2: &'a Path },
}

/// Mate identity for records originating from paired-end ingress.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum MateSide {
    /// First mate in a pair.
    Left,
    /// Second mate in a pair.
    Right,
}

/// Canonical borrowed record view used by preprocessing plans.
#[derive(Clone, Copy)]
pub struct RecordView<'a> {
    header: &'a [u8],
    sequence: &'a [u8],
    quality: &'a [u8],
    provenance: Option<RecordProvenance<'a>>,
}

impl SequenceRecordRef for RecordView<'_> {
    fn header(&self) -> &[u8] {
        self.header
    }

    fn sequence(&self) -> &[u8] {
        self.sequence
    }

    fn quality(&self) -> Option<&[u8]> {
        Some(self.quality)
    }
}

impl<'a> RecordView<'a> {
    /// Construct a new borrowed FASTQ record view with no attached provenance.
    pub fn new(header: &'a [u8], sequence: &'a [u8], quality: &'a [u8]) -> Self {
        Self {
            header,
            sequence,
            quality,
            provenance: None,
        }
    }

    /// Attach borrowed provenance metadata to this record view.
    pub fn with_provenance(mut self, provenance: RecordProvenance<'a>) -> Self {
        self.provenance = Some(provenance);
        self
    }

    /// Return the raw record header bytes without the leading `@` marker.
    pub fn header(&self) -> &'a [u8] {
        self.header
    }

    /// Return the sequence bytes.
    pub fn sequence(&self) -> &'a [u8] {
        self.sequence
    }

    /// Return the quality bytes.
    pub fn quality(&self) -> &'a [u8] {
        self.quality
    }

    /// Return attached provenance metadata when available.
    pub fn provenance(&self) -> Option<RecordProvenance<'a>> {
        self.provenance
    }

    /// Return a new record view with updated sequence and quality slices while preserving metadata.
    ///
    /// # Errors
    ///
    /// Returns an error when the provided sequence and quality slices are not the same length.
    pub fn with_sequence_and_quality(self, sequence: &'a [u8], quality: &'a [u8]) -> Result<Self> {
        if sequence.len() != quality.len() {
            return Err(InternalError::ReplacementLength {
                header: String::from_utf8_lossy(self.header).into_owned(),
                sequence_len: sequence.len(),
                quality_len: quality.len(),
            }
            .into());
        }

        Ok(Self {
            header: self.header,
            sequence,
            quality,
            provenance: self.provenance,
        })
    }

    /// Validate one FASTQ record and either admit it, drop it, or error according to policy.
    pub fn validate(
        self,
        policy: InvalidFastqPolicy,
        stats: &mut ReadStats,
    ) -> Result<Option<Self>> {
        if self.sequence.len() == self.quality.len() {
            return Ok(Some(self));
        }

        stats.record_invalid_read(policy, |context| context.length_mismatch(self))?;

        match policy {
            InvalidFastqPolicy::Error => Err(self.record_length_error()),
            InvalidFastqPolicy::WarnDrop => {
                self.warn_invalid_record(stats);
                Ok(None)
            }
            InvalidFastqPolicy::SilentDrop => Ok(None),
        }
    }

    /// Validate two FASTQ mates together, including structural checks and paired ID agreement.
    pub fn validate_pair(
        self,
        mate: Self,
        policy: InvalidFastqPolicy,
        stats: &mut ReadStats,
    ) -> Result<Option<RecordPair<'a>>> {
        let left = self.validate(policy, stats)?;
        let right = mate.validate(policy, stats)?;

        let (Some(left), Some(right)) = (left, right) else {
            stats.record_invalid_pair();
            return Ok(None);
        };

        if left.pair_key() == right.pair_key() {
            Ok(Some(RecordPair { left, right }))
        } else {
            stats.record_invalid_pair_with_event(policy, |context| {
                context.paired_header_mismatch(left, right)
            })?;
            match policy {
                InvalidFastqPolicy::Error => Err(left.mate_identifier_error(right)),
                InvalidFastqPolicy::WarnDrop => {
                    left.warn_invalid_pair(right, stats);
                    Ok(None)
                }
                InvalidFastqPolicy::SilentDrop => Ok(None),
            }
        }
    }

    fn warn_invalid_record(self, stats: &mut ReadStats) {
        if stats.should_emit_invalid_fastq_warning(INVALID_FASTQ_WARNING_LIMIT) {
            warn!(
                source = %self.source_display(),
                mate = %self.mate_display(),
                header = %String::from_utf8_lossy(self.header()),
                sequence_len = self.sequence.len(),
                quality_len = self.quality.len(),
                "dropping invalid FASTQ record with mismatched sequence and quality lengths"
            );
        } else if stats.should_emit_invalid_fastq_suppressed_notice() {
            warn!("further invalid FASTQ warnings suppressed");
        }
    }

    fn record_length_error(self) -> RunError {
        if let Some(accession) = self.ena_accession() {
            return IndeterminateInputError::Content {
                accession: accession.to_owned(),
                problem: EnaContentProblem::RecordLength {
                    mate: self.mate_display(),
                    header: String::from_utf8_lossy(self.header()).into_owned(),
                    sequence_len: self.sequence.len(),
                    quality_len: self.quality.len(),
                },
            }
            .into();
        }

        MalformedInputError::RecordLength {
            source_label: self.source_display(),
            mate: self.mate_display(),
            header: String::from_utf8_lossy(self.header()).into_owned(),
            sequence_len: self.sequence.len(),
            quality_len: self.quality.len(),
        }
        .into()
    }

    fn mate_identifier_error(self, mate: Self) -> RunError {
        if let Some(accession) = self.ena_accession() {
            return IndeterminateInputError::Content {
                accession: accession.to_owned(),
                problem: EnaContentProblem::MateIdentifier {
                    left_header: String::from_utf8_lossy(self.header()).into_owned(),
                    right_header: String::from_utf8_lossy(mate.header()).into_owned(),
                },
            }
            .into();
        }

        MalformedInputError::MateIdentifier {
            source_label: self.source_display(),
            left_mate: self.mate_display(),
            right_mate: mate.mate_display(),
            left_header: String::from_utf8_lossy(self.header()).into_owned(),
            right_header: String::from_utf8_lossy(mate.header()).into_owned(),
        }
        .into()
    }

    fn ena_accession(self) -> Option<&'a str> {
        match self.provenance() {
            Some(RecordProvenance {
                source: InputSource::Ena { accession },
                ..
            }) => Some(accession),
            _ => None,
        }
    }

    fn warn_invalid_pair(self, mate: Self, stats: &mut ReadStats) {
        if stats.should_emit_invalid_fastq_warning(INVALID_FASTQ_WARNING_LIMIT) {
            warn!(
                source = %self.source_display(),
                left_mate = %self.mate_display(),
                right_mate = %mate.mate_display(),
                left_header = %String::from_utf8_lossy(self.header()),
                right_header = %String::from_utf8_lossy(mate.header()),
                "dropping invalid FASTQ pair with mismatched mate identifiers"
            );
        } else if stats.should_emit_invalid_fastq_suppressed_notice() {
            warn!("further invalid FASTQ warnings suppressed");
        }
    }

    pub(crate) fn pair_key(&self) -> &'a [u8] {
        let first_token = self
            .header()
            .split(u8::is_ascii_whitespace)
            .next()
            .unwrap_or(self.header());
        match first_token {
            [prefix @ .., b'/', b'1' | b'2'] => prefix,
            _ => first_token,
        }
    }

    fn source_display(&self) -> String {
        match self.provenance() {
            Some(RecordProvenance {
                source: InputSource::Ena { accession },
                ..
            }) => format!("ena:{accession}"),
            Some(RecordProvenance {
                source: InputSource::LocalSingle { input },
                ..
            }) => format!("local:{}", input.display()),
            Some(RecordProvenance {
                source: InputSource::LocalInterleavedPaired { input },
                ..
            }) => format!("local-interleaved:{}", input.display()),
            Some(RecordProvenance {
                source: InputSource::LocalPaired { input1, input2 },
                ..
            }) => format!("local-paired:{}|{}", input1.display(), input2.display()),
            None => "unknown".to_owned(),
        }
    }

    fn mate_display(&self) -> &'static str {
        match self.provenance().and_then(|provenance| provenance.mate) {
            Some(MateSide::Left) => "left",
            Some(MateSide::Right) => "right",
            None => "single",
        }
    }
}

fn classify_json_error(path: &Path, source: serde_json::Error) -> RunError {
    if source.io_error_kind().is_some() {
        IoError::WriteReport {
            report_kind: "invalid FASTQ JSONL",
            path: path.to_path_buf(),
            source: ReportWriteError::Json(source),
        }
        .into()
    } else {
        InternalError::SerializeReport {
            report_kind: "invalid FASTQ JSONL",
            source,
        }
        .into()
    }
}

/// Running counters describing what the tool has observed and emitted so far.
#[derive(Debug, Default)]
pub struct ReadStats {
    /// Number of records observed from ingress.
    pub reads_seen: u64,
    /// Number of records emitted to the configured output.
    pub reads_emitted: u64,
    /// Number of records rejected by preprocessing.
    pub reads_rejected: u64,
    /// Number of bases observed from ingress.
    pub bases_seen: u64,
    /// Number of bases emitted to the configured output.
    pub bases_emitted: u64,
    /// Number of paired record groups observed.
    pub pairs_seen: u64,
    /// Number of paired record groups emitted.
    pub pairs_emitted: u64,
    /// Number of paired record groups rejected.
    pub pairs_rejected: u64,
    /// Number of records dropped at ingress because the FASTQ record was malformed.
    pub invalid_reads: u64,
    /// Number of paired record groups dropped at ingress because one or both mates were malformed.
    pub invalid_pairs: u64,
    /// Per-reason rejection counts keyed by stable rejection code.
    pub rejection_counts: BTreeMap<&'static str, u64>,
    /// Per-transform counts keyed by stable transform code.
    pub transform_counts: BTreeMap<&'static str, u64>,
    /// Number of invalid FASTQ warnings already emitted during this run.
    pub invalid_fastq_warnings_emitted: u64,
    /// Whether the warning-suppressed notice has already been emitted.
    pub invalid_fastq_warnings_suppressed: bool,
    /// First invalid FASTQ events observed during this run.
    pub invalid_fastq_samples: Vec<InvalidFastqEvent>,
    /// Whether invalid FASTQ sample storage hit its bounded in-memory limit.
    pub invalid_fastq_samples_truncated: bool,
    /// Optional JSONL report writer for every invalid FASTQ event.
    pub invalid_fastq_report: Option<InvalidFastqReport>,
}

impl ReadStats {
    /// Increment counters for one observed record of `bases` length.
    pub fn record_seen(&mut self, bases: usize) {
        self.reads_seen += 1;
        self.bases_seen += bases as u64;
    }

    /// Increment counters for one emitted record of `bases` length.
    pub fn record_emitted(&mut self, bases: usize) {
        self.reads_emitted += 1;
        self.bases_emitted += bases as u64;
    }

    /// Increment the rejected-record counter for one stable rejection code.
    pub fn record_rejected(&mut self, code: &'static str) {
        self.reads_rejected += 1;
        *self.rejection_counts.entry(code).or_default() += 1;
    }

    /// Increment the emitted paired-record counter.
    pub fn record_pair_emitted(&mut self) {
        self.pairs_emitted += 1;
    }

    /// Increment the rejected paired-record counter.
    pub fn record_pair_rejected(&mut self) {
        self.pairs_rejected += 1;
    }

    /// Increment the invalid-record counter for one malformed FASTQ record.
    pub fn record_invalid_read(
        &mut self,
        policy: InvalidFastqPolicy,
        build: impl FnOnce(InvalidFastqContext) -> InvalidFastqEvent,
    ) -> Result<()> {
        self.invalid_reads += 1;
        let context = self.invalid_fastq_context(policy);
        self.record_invalid_fastq_sample(build(context))
    }

    /// Increment the invalid-record counter for one fatal parser-level FASTQ error.
    pub fn record_invalid_parse_error(
        &mut self,
        policy: InvalidFastqPolicy,
        build: impl FnOnce(InvalidFastqContext) -> InvalidFastqEvent,
    ) -> Result<()> {
        self.invalid_reads += 1;
        let context = self.invalid_fastq_context(policy);
        self.record_invalid_fastq_sample(build(context))
    }

    /// Increment the invalid-pair counter for one malformed paired FASTQ record group.
    pub fn record_invalid_pair(&mut self) {
        self.invalid_pairs += 1;
    }

    /// Increment the invalid-pair counter and remember a representative event.
    pub fn record_invalid_pair_with_event(
        &mut self,
        policy: InvalidFastqPolicy,
        build: impl FnOnce(InvalidFastqContext) -> InvalidFastqEvent,
    ) -> Result<()> {
        self.invalid_pairs += 1;
        let context = self.invalid_fastq_context(policy);
        self.record_invalid_fastq_sample(build(context))
    }

    /// Install a JSONL report writer for invalid FASTQ events.
    pub fn set_invalid_fastq_report(&mut self, report: InvalidFastqReport) {
        self.invalid_fastq_report = Some(report);
    }

    /// Flush the optional invalid-FASTQ report before a successful run returns.
    pub fn finish_invalid_fastq_report(&mut self) -> Result<()> {
        if let Some(report) = &mut self.invalid_fastq_report {
            report.finish()?;
        }
        Ok(())
    }

    /// Increment the application count for one transform.
    pub fn record_transform(&mut self, code: &'static str) {
        *self.transform_counts.entry(code).or_default() += 1;
    }

    fn record_invalid_fastq_sample(&mut self, event: InvalidFastqEvent) -> Result<()> {
        if let Some(report) = &mut self.invalid_fastq_report {
            report.write_event(&event)?;
        }

        if self.invalid_fastq_samples.len() < INVALID_FASTQ_SAMPLE_LIMIT {
            self.invalid_fastq_samples.push(event);
        } else {
            self.invalid_fastq_samples_truncated = true;
        }

        Ok(())
    }

    fn invalid_fastq_context(&self, policy: InvalidFastqPolicy) -> InvalidFastqContext {
        InvalidFastqContext {
            reads_seen: self.reads_seen,
            pairs_seen: (self.pairs_seen > 0).then_some(self.pairs_seen),
            policy,
        }
    }

    /// Return `true` when another invalid FASTQ warning may be emitted.
    pub fn should_emit_invalid_fastq_warning(&mut self, limit: u64) -> bool {
        if self.invalid_fastq_warnings_emitted < limit {
            self.invalid_fastq_warnings_emitted += 1;
            true
        } else {
            false
        }
    }

    /// Return `true` only once, when invalid FASTQ warning suppression should be announced.
    pub fn should_emit_invalid_fastq_suppressed_notice(&mut self) -> bool {
        if self.invalid_fastq_warnings_suppressed {
            false
        } else {
            self.invalid_fastq_warnings_suppressed = true;
            true
        }
    }
}

#[cfg(test)]
mod tests {
    use std::{
        fs,
        io::{self, Write},
        path::Path,
    };

    use tempfile::tempdir;

    use super::{
        InputSource, InvalidFastqContext, InvalidFastqReport, MateSide, ReadStats,
        RecordProvenance, RecordView, finalize_invalid_fastq_report, write_invalid_fastq_event_to,
    };
    use crate::{
        cli::InvalidFastqPolicy,
        error::{
            EnaContentProblem, IndeterminateInputError, IoError, MalformedInputError,
            ReportWriteError, RunError,
        },
    };

    struct FailingWriter {
        fail_write: bool,
        fail_flush: bool,
    }

    struct NewlineFailingWriter;

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

    impl Write for NewlineFailingWriter {
        fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
            if bytes == b"\n" {
                Err(io::Error::new(
                    io::ErrorKind::StorageFull,
                    "test newline failure",
                ))
            } else {
                Ok(bytes.len())
            }
        }

        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }

    fn fatal_event() -> super::InvalidFastqEvent {
        InvalidFastqContext {
            reads_seen: 1,
            pairs_seen: None,
            policy: InvalidFastqPolicy::Error,
        }
        .parse_error(
            "local:reads.fastq",
            "single",
            "UnequalLengths".to_owned(),
            "sequence and quality lengths differ".to_owned(),
            Some(1),
        )
    }

    #[test]
    fn pair_key_strips_slash_mate_suffix() {
        let left = RecordView::new(b"read123/1", b"A", b"I");
        let right = RecordView::new(b"read123/2", b"A", b"I");
        let bare = RecordView::new(b"read123", b"A", b"I");

        assert_eq!(left.pair_key(), b"read123");
        assert_eq!(right.pair_key(), b"read123");
        assert_eq!(bare.pair_key(), b"read123");
    }

    #[test]
    fn validate_pair_accepts_matching_mate_ids() {
        let left = RecordView::new(b"read123/1", b"ACGT", b"IIII");
        let right = RecordView::new(b"read123/2", b"TGCA", b"JJJJ");
        let mut stats = ReadStats::default();

        let pair = left
            .validate_pair(right, InvalidFastqPolicy::Error, &mut stats)
            .expect("pair validation should succeed");
        assert!(pair.is_some());
    }

    #[test]
    fn validate_pair_drops_mismatched_ids_under_drop_policy() {
        let left = RecordView::new(b"read123/1", b"ACGT", b"IIII");
        let right = RecordView::new(b"read999/2", b"TGCA", b"JJJJ");
        let mut stats = ReadStats::default();

        let pair = left
            .validate_pair(right, InvalidFastqPolicy::SilentDrop, &mut stats)
            .expect("drop policy should not error");
        assert!(pair.is_none());
        assert_eq!(stats.invalid_pairs, 1);
    }

    #[test]
    fn record_length_error_classification_respects_origin() {
        let mut local_stats = ReadStats::default();
        let local = RecordView::new(b"read1", b"ACGT", b"I").with_provenance(RecordProvenance {
            source: InputSource::LocalSingle {
                input: Path::new("reads.fastq"),
            },
            mate: None,
        });
        let Err(local_error) = local.validate(InvalidFastqPolicy::Error, &mut local_stats) else {
            panic!("local length mismatch should fail");
        };
        assert!(matches!(
            local_error,
            RunError::MalformedInput(error)
                if matches!(*error, MalformedInputError::RecordLength { .. })
        ));

        let mut ena_stats = ReadStats::default();
        let ena = RecordView::new(b"read1", b"ACGT", b"I").with_provenance(RecordProvenance {
            source: InputSource::Ena {
                accession: "SRR35939766",
            },
            mate: None,
        });
        let Err(ena_error) = ena.validate(InvalidFastqPolicy::Error, &mut ena_stats) else {
            panic!("ENA length mismatch should fail");
        };
        assert!(matches!(
            ena_error,
            RunError::IndeterminateInput(IndeterminateInputError::Content {
                problem: EnaContentProblem::RecordLength { .. },
                ..
            })
        ));
    }

    #[test]
    fn mate_identifier_error_classification_respects_origin() {
        let mut local_stats = ReadStats::default();
        let local_provenance = RecordProvenance {
            source: InputSource::LocalInterleavedPaired {
                input: Path::new("reads.fastq"),
            },
            mate: Some(MateSide::Left),
        };
        let local_left = RecordView::new(b"read1/1", b"A", b"I").with_provenance(local_provenance);
        let local_right =
            RecordView::new(b"read2/2", b"T", b"I").with_provenance(RecordProvenance {
                mate: Some(MateSide::Right),
                ..local_provenance
            });
        let Err(local_error) =
            local_left.validate_pair(local_right, InvalidFastqPolicy::Error, &mut local_stats)
        else {
            panic!("local mate mismatch should fail");
        };
        assert!(matches!(
            local_error,
            RunError::MalformedInput(error)
                if matches!(*error, MalformedInputError::MateIdentifier { .. })
        ));

        let mut ena_stats = ReadStats::default();
        let ena_left = RecordView::new(b"read1/1", b"A", b"I").with_provenance(RecordProvenance {
            source: InputSource::Ena {
                accession: "SRR35939766",
            },
            mate: Some(MateSide::Left),
        });
        let ena_right = RecordView::new(b"read2/2", b"T", b"I").with_provenance(RecordProvenance {
            source: InputSource::Ena {
                accession: "SRR35939766",
            },
            mate: Some(MateSide::Right),
        });
        let Err(ena_error) =
            ena_left.validate_pair(ena_right, InvalidFastqPolicy::Error, &mut ena_stats)
        else {
            panic!("ENA mate mismatch should fail");
        };
        assert!(matches!(
            ena_error,
            RunError::IndeterminateInput(IndeterminateInputError::Content {
                problem: EnaContentProblem::MateIdentifier { .. },
                ..
            })
        ));
    }

    #[test]
    fn with_sequence_and_quality_preserves_provenance() {
        let record =
            RecordView::new(b"read123/1", b"ACGT", b"IIII").with_provenance(RecordProvenance {
                source: InputSource::Ena {
                    accession: "ERR000002",
                },
                mate: Some(MateSide::Left),
            });

        let rewritten = record
            .with_sequence_and_quality(b"AC", b"II")
            .expect("replacement slices should be valid");

        assert_eq!(
            rewritten
                .provenance()
                .expect("rewritten record should preserve provenance")
                .mate,
            Some(MateSide::Left)
        );
    }

    #[test]
    fn invalid_fastq_report_flushes_fatal_events_before_drop() {
        let temp = tempdir().expect("tempdir should be created");
        let path = temp.path().join("invalid-fastq.jsonl");
        let mut stats = ReadStats::default();
        stats.set_invalid_fastq_report(
            InvalidFastqReport::create(&path).expect("invalid FASTQ report should be created"),
        );

        stats
            .record_invalid_parse_error(InvalidFastqPolicy::SilentDrop, |context| {
                context.parse_error(
                    "ena:SRR000001",
                    "single",
                    "UnequalLengths".to_owned(),
                    "Unequal length: sequence length is 4 while quality length is 1".to_owned(),
                    Some(1),
                )
            })
            .expect("fatal invalid FASTQ event should be reported");

        let report = fs::read_to_string(path)
            .expect("fatal invalid FASTQ report should be readable before stats is dropped");
        assert!(report.contains("\"kind\":\"fastq_parse_error\""));
        assert!(report.contains("\"fatal\":true"));
        assert!(report.ends_with('\n'));
    }

    #[test]
    fn invalid_fastq_report_flushes_nonfatal_events_on_finish() {
        let temp = tempdir().expect("tempdir should be created");
        let path = temp.path().join("invalid-fastq.jsonl");
        let mut stats = ReadStats::default();
        stats.set_invalid_fastq_report(
            InvalidFastqReport::create(&path).expect("invalid FASTQ report should be created"),
        );

        let dropped = RecordView::new(b"bad", b"AAAA", b"I")
            .validate(InvalidFastqPolicy::SilentDrop, &mut stats)
            .expect("drop policy should retain a safe record boundary");
        assert!(dropped.is_none());
        stats
            .finish_invalid_fastq_report()
            .expect("successful run should flush its invalid FASTQ report");

        let report = fs::read_to_string(path)
            .expect("finished invalid FASTQ report should be immediately readable");
        assert!(report.contains("\"kind\":\"sequence_quality_length_mismatch\""));
        assert!(report.ends_with('\n'));
    }

    #[test]
    fn invalid_fastq_report_write_failure_retains_path_and_json_source() {
        let path = Path::new("invalid-fastq.jsonl");
        let mut writer = FailingWriter {
            fail_write: true,
            fail_flush: false,
        };

        let error = write_invalid_fastq_event_to(&mut writer, path, &fatal_event())
            .expect_err("requested report write failure should be required I/O");
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
    fn invalid_fastq_report_flush_failure_is_finalization() {
        let path = Path::new("invalid-fastq.jsonl");
        let mut writer = FailingWriter {
            fail_write: false,
            fail_flush: true,
        };

        let error = finalize_invalid_fastq_report(&mut writer, path)
            .expect_err("requested report flush failure should be required I/O");
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
    fn invalid_fastq_report_newline_failure_retains_io_source() {
        let path = Path::new("invalid-fastq.jsonl");
        let mut writer = NewlineFailingWriter;

        let error = write_invalid_fastq_event_to(&mut writer, path, &fatal_event())
            .expect_err("JSONL newline failure should be required I/O");
        assert!(matches!(
            error,
            RunError::Io(IoError::WriteReport {
                path: observed_path,
                source: ReportWriteError::Bytes(source),
                ..
            }) if observed_path == path && source.kind() == io::ErrorKind::StorageFull
        ));
    }

    #[test]
    fn invalid_fastq_report_create_failure_retains_requested_path() {
        let temp = tempdir().expect("tempdir should be created");
        let path = temp
            .path()
            .join("missing-parent")
            .join("invalid-fastq.jsonl");

        let error = InvalidFastqReport::create(&path)
            .expect_err("report under missing parent should fail to open");
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
