//! Core record and accounting types shared across ingress, parsing, and output layers.

use std::{
    collections::BTreeMap,
    fs::File,
    io::{BufWriter, Write},
    path::{Path, PathBuf},
};

use serde::Serialize;

use crate::error::{InternalError, IoError, ReportWriteError, Result, RunError};

const ADMISSION_SAMPLE_LIMIT: usize = 20;
const ADMISSION_REPORT_KIND: &str = "invalid-input JSONL";

/// Bounded, owned trace information about one record-admission failure.
#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum AdmissionEvent {
    /// A complete record whose sequence and quality lengths differ.
    SequenceQualityLengthMismatch {
        /// Source label such as `ena:SRR...` or a local path label.
        source: String,
        /// Mate label for this record.
        mate: &'static str,
        /// Record header.
        header: String,
        /// Sequence length.
        sequence_len: usize,
        /// Quality length.
        quality_len: usize,
        /// Total records observed when the event was detected.
        reads_seen: u64,
        /// Total complete pairs observed for paired input.
        pairs_seen: Option<u64>,
        /// Whether processing continued beyond this event.
        continued: bool,
    },
    /// Two complete mate records whose normalized identifiers disagree.
    PairIdentifierMismatch {
        /// Source label.
        source: String,
        /// Left record header.
        left_header: String,
        /// Right record header.
        right_header: String,
        /// Total records observed when the event was detected.
        reads_seen: u64,
        /// Total complete pairs observed when the event was detected.
        pairs_seen: u64,
        /// Whether processing continued beyond this event.
        continued: bool,
    },
    /// A complete record whose corresponding input mate is absent.
    MissingMate {
        /// Source label.
        source: String,
        /// Which mate was present.
        present_mate: MateSide,
        /// Present record's header.
        header: String,
        /// Total records observed when the event was detected.
        reads_seen: u64,
        /// Total complete pairs observed when the event was detected.
        pairs_seen: u64,
        /// Whether processing continued beyond this event.
        continued: bool,
    },
    /// A non-I/O parser failure encountered while attempting the next record.
    RecordParseFailure {
        /// Source label.
        source: String,
        /// Mate label for the parser.
        mate: &'static str,
        /// Stable parser error kind.
        parser_kind: String,
        /// Parser diagnostic.
        message: String,
        /// Parser-reported input line.
        line: Option<u64>,
        /// Total records observed when the event was detected.
        reads_seen: u64,
        /// Total complete pairs observed for paired input.
        pairs_seen: Option<u64>,
        /// Whether processing continued beyond this event.
        continued: bool,
    },
}

impl AdmissionEvent {
    pub(crate) fn kind(&self) -> &'static str {
        match self {
            Self::SequenceQualityLengthMismatch { .. } => "sequence_quality_length_mismatch",
            Self::PairIdentifierMismatch { .. } => "pair_identifier_mismatch",
            Self::MissingMate { .. } => "missing_mate",
            Self::RecordParseFailure { .. } => "record_parse_failure",
        }
    }
}

/// Newline-delimited JSON writer for record-admission failures.
#[derive(Debug)]
pub struct AdmissionReport<W: Write = BufWriter<File>> {
    path: PathBuf,
    writer: W,
}

impl AdmissionReport<BufWriter<File>> {
    /// Create a JSONL report at `path`.
    ///
    /// # Errors
    ///
    /// Returns an error when the destination file cannot be created.
    pub fn create(path: &Path) -> Result<Self> {
        let writer = File::create(path).map_err(|source| IoError::CreateReport {
            report_kind: ADMISSION_REPORT_KIND,
            path: path.to_path_buf(),
            source,
        })?;
        Ok(Self {
            path: path.to_path_buf(),
            writer: BufWriter::new(writer),
        })
    }
}

impl<W: Write> AdmissionReport<W> {
    #[cfg(test)]
    fn from_writer(path: impl Into<PathBuf>, writer: W) -> Self {
        Self {
            path: path.into(),
            writer,
        }
    }

    fn write_event(&mut self, event: &AdmissionEvent) -> Result<()> {
        write_admission_event_to(&mut self.writer, &self.path, event)
    }

    fn finish(&mut self) -> Result<()> {
        finalize_admission_report(&mut self.writer, &self.path)
    }
}

fn write_admission_event_to(
    writer: &mut impl Write,
    path: &Path,
    event: &AdmissionEvent,
) -> Result<()> {
    serde_json::to_writer(&mut *writer, event)
        .map_err(|source| classify_json_error(path, source))?;
    writer
        .write_all(b"\n")
        .map_err(|source| IoError::WriteReport {
            report_kind: ADMISSION_REPORT_KIND,
            path: path.to_path_buf(),
            source: ReportWriteError::Bytes(source),
        })?;
    writer.flush().map_err(|source| IoError::FinalizeReport {
        report_kind: ADMISSION_REPORT_KIND,
        path: path.to_path_buf(),
        source,
    })?;
    Ok(())
}

fn finalize_admission_report(writer: &mut impl Write, path: &Path) -> Result<()> {
    writer.flush().map_err(|source| {
        IoError::FinalizeReport {
            report_kind: ADMISSION_REPORT_KIND,
            path: path.to_path_buf(),
            source,
        }
        .into()
    })
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
#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize)]
#[serde(rename_all = "snake_case")]
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

    pub(crate) fn source_display(&self) -> String {
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

    pub(crate) fn mate_display(&self) -> &'static str {
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
            report_kind: ADMISSION_REPORT_KIND,
            path: path.to_path_buf(),
            source: ReportWriteError::Json(source),
        }
        .into()
    } else {
        InternalError::SerializeReport {
            report_kind: ADMISSION_REPORT_KIND,
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
    /// Per-kind record-admission failure counts.
    pub admission_event_counts: BTreeMap<&'static str, u64>,
    /// Number of record-admission warnings already emitted during this run.
    pub admission_warnings_emitted: u64,
    /// Whether the warning-suppressed notice has already been emitted.
    pub admission_warnings_suppressed: bool,
    /// First record-admission failures observed during this run.
    pub admission_samples: Vec<AdmissionEvent>,
    /// Whether admission sample storage hit its bounded in-memory limit.
    pub admission_samples_truncated: bool,
    /// Optional JSONL report writer for record-admission failures.
    pub admission_report: Option<AdmissionReport>,
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

    /// Install a JSONL report writer for record-admission failures.
    pub fn set_admission_report(&mut self, report: AdmissionReport) {
        self.admission_report = Some(report);
    }

    /// Increment the application count for one transform.
    pub fn record_transform(&mut self, code: &'static str) {
        *self.transform_counts.entry(code).or_default() += 1;
    }

    /// Record one admission failure in counters, the detailed report, and bounded samples.
    pub fn record_admission_event(&mut self, event: AdmissionEvent) -> Result<()> {
        *self.admission_event_counts.entry(event.kind()).or_default() += 1;

        if let Some(report) = &mut self.admission_report {
            report.write_event(&event)?;
        }

        if self.admission_samples.len() < ADMISSION_SAMPLE_LIMIT {
            self.admission_samples.push(event);
        } else {
            self.admission_samples_truncated = true;
        }

        Ok(())
    }

    /// Flush the admission report on successful completion.
    pub fn finish_admission_report(&mut self) -> Result<()> {
        if let Some(report) = &mut self.admission_report {
            report.finish()?;
        }
        Ok(())
    }

    /// Return `true` when another record-admission warning may be emitted.
    pub fn should_emit_admission_warning(&mut self, limit: u64) -> bool {
        if self.admission_warnings_emitted < limit {
            self.admission_warnings_emitted += 1;
            true
        } else {
            false
        }
    }

    /// Return `true` only once, when admission-warning suppression should be announced.
    pub fn should_emit_admission_suppressed_notice(&mut self) -> bool {
        if self.admission_warnings_suppressed {
            false
        } else {
            self.admission_warnings_suppressed = true;
            true
        }
    }
}

#[cfg(test)]
mod tests {
    use std::{
        fs, io,
        io::Write,
        path::{Path, PathBuf},
    };

    use tempfile::tempdir;

    use super::{
        AdmissionEvent, AdmissionReport, InputSource, MateSide, ReadStats, RecordProvenance,
        RecordView, classify_json_error, finalize_admission_report, write_admission_event_to,
    };
    use crate::error::{IoError, ReportWriteError, RunError};

    #[derive(Debug)]
    struct FailingWriter {
        fail_write: bool,
        fail_flush: bool,
    }

    #[derive(Debug)]
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

    fn terminal_event() -> AdmissionEvent {
        AdmissionEvent::RecordParseFailure {
            source: "local:reads.fastq".to_owned(),
            mate: "single",
            parser_kind: "unequal_lengths".to_owned(),
            message: "sequence and quality lengths differ".to_owned(),
            line: Some(1),
            reads_seen: 1,
            pairs_seen: None,
            continued: false,
        }
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
    fn admission_report_flushes_terminal_events_before_drop() {
        let temp = tempdir().expect("tempdir should be created");
        let path = temp.path().join("invalid-input.jsonl");
        let mut stats = ReadStats::default();
        stats.set_admission_report(
            AdmissionReport::create(&path).expect("admission report should be created"),
        );

        stats
            .record_admission_event(AdmissionEvent::RecordParseFailure {
                source: "ena:SRR000001".to_owned(),
                mate: "single",
                parser_kind: "unequal_lengths".to_owned(),
                message: "Sequence length is 4 but quality length is 1".to_owned(),
                line: Some(1),
                reads_seen: 0,
                pairs_seen: None,
                continued: false,
            })
            .expect("terminal admission event should be reported");

        let report = fs::read_to_string(path)
            .expect("terminal admission report should be readable before stats is dropped");
        assert!(report.contains("\"kind\":\"record_parse_failure\""));
        assert!(report.contains("\"continued\":false"));
        assert!(report.ends_with('\n'));
    }

    #[test]
    fn admission_report_flushes_nonterminal_events_immediately() {
        let temp = tempdir().expect("tempdir should be created");
        let path = temp.path().join("invalid-input.jsonl");
        let mut stats = ReadStats::default();
        stats.set_admission_report(
            AdmissionReport::create(&path).expect("admission report should be created"),
        );

        stats
            .record_admission_event(AdmissionEvent::SequenceQualityLengthMismatch {
                source: "local:reads.fastq".to_owned(),
                mate: "single",
                header: "bad".to_owned(),
                sequence_len: 4,
                quality_len: 1,
                reads_seen: 1,
                pairs_seen: None,
                continued: true,
            })
            .expect("nonterminal admission event should be reported");

        let report = fs::read_to_string(path)
            .expect("nonterminal admission report should be readable before finish");
        assert!(report.contains("\"kind\":\"sequence_quality_length_mismatch\""));
        assert!(report.ends_with('\n'));
    }

    #[test]
    fn admission_report_write_failure_retains_path_and_json_source() {
        let path = Path::new("invalid-input.jsonl");
        let mut writer = FailingWriter {
            fail_write: true,
            fail_flush: false,
        };

        let error = write_admission_event_to(&mut writer, path, &terminal_event())
            .expect_err("requested report write failure should be required I/O");
        assert!(matches!(
            error,
            RunError::Io(IoError::WriteReport {
                report_kind: "invalid-input JSONL",
                path: observed_path,
                source: ReportWriteError::Json(source),
            }) if observed_path == path && source.io_error_kind() == Some(io::ErrorKind::StorageFull)
        ));
    }

    #[test]
    fn admission_report_newline_failure_retains_io_source() {
        let path = Path::new("invalid-input.jsonl");
        let mut writer = NewlineFailingWriter;

        let error = write_admission_event_to(&mut writer, path, &terminal_event())
            .expect_err("JSONL newline failure should be required I/O");
        assert!(matches!(
            error,
            RunError::Io(IoError::WriteReport {
                report_kind: "invalid-input JSONL",
                path: observed_path,
                source: ReportWriteError::Bytes(source),
            }) if observed_path == path && source.kind() == io::ErrorKind::StorageFull
        ));
    }

    #[test]
    fn admission_report_event_flush_failure_is_finalization() {
        let path = Path::new("invalid-input.jsonl");
        let mut writer = FailingWriter {
            fail_write: false,
            fail_flush: true,
        };

        let error = write_admission_event_to(&mut writer, path, &terminal_event())
            .expect_err("event flush failure should be required I/O");
        assert!(matches!(
            error,
            RunError::Io(IoError::FinalizeReport {
                report_kind: "invalid-input JSONL",
                path: observed_path,
                source,
            }) if observed_path == path && source.kind() == io::ErrorKind::StorageFull
        ));
    }

    #[test]
    fn admission_report_final_flush_failure_is_finalization() {
        let path = Path::new("invalid-input.jsonl");
        let mut writer = FailingWriter {
            fail_write: false,
            fail_flush: true,
        };

        let error = finalize_admission_report(&mut writer, path)
            .expect_err("requested report flush failure should be required I/O");
        assert!(matches!(
            error,
            RunError::Io(IoError::FinalizeReport {
                report_kind: "invalid-input JSONL",
                path: observed_path,
                source,
            }) if observed_path == path && source.kind() == io::ErrorKind::StorageFull
        ));
    }

    #[test]
    fn admission_report_create_failure_retains_requested_path() {
        let temp = tempdir().expect("tempdir should be created");
        let path = temp
            .path()
            .join("missing-parent")
            .join("invalid-input.jsonl");

        let error = AdmissionReport::create(&path)
            .expect_err("report under missing parent should fail to open");
        assert!(matches!(
            error,
            RunError::Io(IoError::CreateReport {
                report_kind: "invalid-input JSONL",
                path: observed_path,
                source,
            }) if observed_path == path && source.kind() == io::ErrorKind::NotFound
        ));
    }

    #[test]
    fn admission_report_from_writer_uses_fixture_path_for_typed_errors() {
        let path = PathBuf::from("fixture-invalid-input.jsonl");
        let mut report = AdmissionReport::from_writer(
            path.clone(),
            FailingWriter {
                fail_write: false,
                fail_flush: true,
            },
        );

        let error = report
            .write_event(&terminal_event())
            .expect_err("fixture writer flush failure should retain fixture path");
        assert!(matches!(
            error,
            RunError::Io(IoError::FinalizeReport {
                report_kind: "invalid-input JSONL",
                path: observed_path,
                source,
            }) if observed_path == path && source.kind() == io::ErrorKind::StorageFull
        ));
    }

    #[test]
    fn non_io_serialization_defect_is_internal() {
        let source = serde_json::from_str::<serde_json::Value>("{")
            .expect_err("truncated JSON should produce a non-I/O serde_json error");
        let error = classify_json_error(Path::new("invalid-input.jsonl"), source);

        assert!(matches!(error, RunError::Internal(_)));
    }

    #[test]
    fn admission_samples_are_limited_to_twenty_with_deterministic_counts() {
        let mut stats = ReadStats::default();

        for index in 0..25 {
            stats
                .record_admission_event(AdmissionEvent::MissingMate {
                    source: "local-paired:left.fastq|right.fastq".to_owned(),
                    present_mate: MateSide::Left,
                    header: format!("read{index}"),
                    reads_seen: index + 1,
                    pairs_seen: index,
                    continued: true,
                })
                .expect("admission event should be recorded");
        }

        assert_eq!(stats.admission_samples.len(), 20);
        assert!(stats.admission_samples_truncated);
        assert_eq!(stats.admission_event_counts["missing_mate"], 25);
        assert_eq!(stats.invalid_reads, 0);
        assert_eq!(stats.invalid_pairs, 0);
    }
}
