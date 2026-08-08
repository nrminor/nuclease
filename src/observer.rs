//! Run observation, invalid-input reporting, and bounded diagnostic sampling.

use std::{
    collections::BTreeMap,
    fs::File,
    io::{BufWriter, Write},
    path::{Path, PathBuf},
};

use needletail::errors::ParseError;
use serde::{Serialize, Serializer, ser::SerializeStruct};

use crate::{
    cli::AdmissionPolicy,
    error::{self, IoError, PairConstructionError, ReportWriteError, Result},
    record::MateSide,
};

const INVALID_INPUT_SAMPLE_LIMIT: usize = 20;
const INVALID_INPUT_REPORT_KIND: &str = "invalid-input JSONL";

/// Parser-owned evidence captured from one failed physical FASTQ slot.
#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct RecordParseFailure {
    /// Stable parser error kind.
    pub parser_kind: String,
    /// Parser diagnostic.
    pub message: String,
    /// Parser-reported input line.
    pub line: Option<u64>,
}

impl From<&ParseError> for RecordParseFailure {
    fn from(error: &ParseError) -> Self {
        Self {
            parser_kind: error::constructors::parse_error_kind_name(&error.kind).to_owned(),
            message: error.to_string(),
            line: (error.position.line > 0).then_some(error.position.line),
        }
    }
}

/// Bounded, owned trace information about one invalid-input event.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum InvalidInputEvent {
    /// A complete record whose sequence and quality lengths differ.
    SequenceQualityLengthMismatch {
        /// Source label such as `ena:SRR...` or a local path label.
        source: String,
        /// Runtime mate identity for paired input; absent for single input.
        mate: Option<MateSide>,
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
    /// A typed failure to establish a pair from two positional mate claims.
    PairConstructionFailure {
        /// Source label.
        source: String,
        /// Authoritative pair-domain failure.
        error: PairConstructionError,
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
    /// A parser failure while attempting the next single-end record.
    SingleRecordParseFailure {
        /// Source label.
        source: String,
        /// Parser-owned failure evidence.
        failure: RecordParseFailure,
        /// Total records observed when the event was detected.
        reads_seen: u64,
        /// Whether processing continued beyond this event.
        continued: bool,
    },
    /// A parser failure while attempting one mate of the next pair.
    PairedRecordParseFailure {
        /// Source label.
        source: String,
        /// Runtime mate identity.
        mate: MateSide,
        /// Parser-owned failure evidence.
        failure: RecordParseFailure,
        /// Total records observed when the event was detected.
        reads_seen: u64,
        /// Total complete pairs observed when the event was detected.
        pairs_seen: u64,
        /// Whether processing continued beyond this event.
        continued: bool,
    },
}

impl Serialize for InvalidInputEvent {
    fn serialize<S>(&self, serializer: S) -> std::result::Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        match self {
            Self::SequenceQualityLengthMismatch {
                source,
                mate,
                header,
                sequence_len,
                quality_len,
                reads_seen,
                pairs_seen,
                continued,
            } => {
                let field_count =
                    7 + usize::from(mate.is_some()) + usize::from(pairs_seen.is_some());
                let mut event = serializer.serialize_struct("InvalidInputEvent", field_count)?;
                event.serialize_field("kind", self.kind())?;
                event.serialize_field("source", source)?;
                if let Some(mate) = mate {
                    event.serialize_field("mate", mate)?;
                }
                event.serialize_field("header", header)?;
                event.serialize_field("sequence_len", sequence_len)?;
                event.serialize_field("quality_len", quality_len)?;
                event.serialize_field("reads_seen", reads_seen)?;
                if let Some(pairs_seen) = pairs_seen {
                    event.serialize_field("pairs_seen", pairs_seen)?;
                }
                event.serialize_field("continued", continued)?;
                event.end()
            }
            Self::PairConstructionFailure {
                source,
                error,
                reads_seen,
                pairs_seen,
                continued,
            } => {
                let mut event = serializer.serialize_struct("InvalidInputEvent", 6)?;
                event.serialize_field("kind", self.kind())?;
                event.serialize_field("source", source)?;
                event.serialize_field("error", error)?;
                event.serialize_field("reads_seen", reads_seen)?;
                event.serialize_field("pairs_seen", pairs_seen)?;
                event.serialize_field("continued", continued)?;
                event.end()
            }
            Self::MissingMate {
                source,
                present_mate,
                header,
                reads_seen,
                pairs_seen,
                continued,
            } => {
                let mut event = serializer.serialize_struct("InvalidInputEvent", 7)?;
                event.serialize_field("kind", self.kind())?;
                event.serialize_field("source", source)?;
                event.serialize_field("present_mate", present_mate)?;
                event.serialize_field("header", header)?;
                event.serialize_field("reads_seen", reads_seen)?;
                event.serialize_field("pairs_seen", pairs_seen)?;
                event.serialize_field("continued", continued)?;
                event.end()
            }
            Self::SingleRecordParseFailure {
                source,
                failure,
                reads_seen,
                continued,
            } => {
                let mut event = serializer.serialize_struct("InvalidInputEvent", 5)?;
                event.serialize_field("kind", self.kind())?;
                event.serialize_field("source", source)?;
                event.serialize_field("failure", failure)?;
                event.serialize_field("reads_seen", reads_seen)?;
                event.serialize_field("continued", continued)?;
                event.end()
            }
            Self::PairedRecordParseFailure {
                source,
                mate,
                failure,
                reads_seen,
                pairs_seen,
                continued,
            } => {
                let mut event = serializer.serialize_struct("InvalidInputEvent", 7)?;
                event.serialize_field("kind", self.kind())?;
                event.serialize_field("source", source)?;
                event.serialize_field("mate", mate)?;
                event.serialize_field("failure", failure)?;
                event.serialize_field("reads_seen", reads_seen)?;
                event.serialize_field("pairs_seen", pairs_seen)?;
                event.serialize_field("continued", continued)?;
                event.end()
            }
        }
    }
}

impl InvalidInputEvent {
    pub(crate) fn kind(&self) -> &'static str {
        match self {
            Self::SequenceQualityLengthMismatch { .. } => "sequence_quality_length_mismatch",
            Self::PairConstructionFailure { error, .. } => error.code(),
            Self::MissingMate { .. } => "missing_mate",
            Self::SingleRecordParseFailure { .. } | Self::PairedRecordParseFailure { .. } => {
                "record_parse_failure"
            }
        }
    }
}

/// Newline-delimited JSON writer for invalid-input events.
#[derive(Debug)]
pub struct InvalidInputReport<W: Write = BufWriter<File>> {
    path: PathBuf,
    writer: W,
}

impl InvalidInputReport<BufWriter<File>> {
    /// Create a JSONL report at `path`.
    pub fn create(path: &Path) -> Result<Self> {
        let writer = File::create(path).map_err(|source| IoError::CreateReport {
            report_kind: INVALID_INPUT_REPORT_KIND,
            path: path.to_path_buf(),
            source,
        })?;
        Ok(Self {
            path: path.to_path_buf(),
            writer: BufWriter::new(writer),
        })
    }
}

impl<W: Write> InvalidInputReport<W> {
    #[cfg(test)]
    pub(crate) fn from_writer(path: impl Into<PathBuf>, writer: W) -> Self {
        Self {
            path: path.into(),
            writer,
        }
    }

    pub(crate) fn write_event(&mut self, event: &InvalidInputEvent) -> Result<()> {
        write_invalid_input_event_to(&mut self.writer, &self.path, event)
    }

    fn finish(&mut self) -> Result<()> {
        finalize_invalid_input_report(&mut self.writer, &self.path)
    }
}

pub(crate) fn write_invalid_input_event_to(
    writer: &mut impl Write,
    path: &Path,
    event: &InvalidInputEvent,
) -> Result<()> {
    serde_json::to_writer(&mut *writer, event).map_err(|source| {
        error::constructors::json_report_error(INVALID_INPUT_REPORT_KIND, path, source)
    })?;
    writer
        .write_all(b"\n")
        .map_err(|source| IoError::WriteReport {
            report_kind: INVALID_INPUT_REPORT_KIND,
            path: path.to_path_buf(),
            source: ReportWriteError::Bytes(source),
        })?;
    writer.flush().map_err(|source| IoError::FinalizeReport {
        report_kind: INVALID_INPUT_REPORT_KIND,
        path: path.to_path_buf(),
        source,
    })?;
    Ok(())
}

pub(crate) fn finalize_invalid_input_report(writer: &mut impl Write, path: &Path) -> Result<()> {
    writer.flush().map_err(|source| {
        IoError::FinalizeReport {
            report_kind: INVALID_INPUT_REPORT_KIND,
            path: path.to_path_buf(),
            source,
        }
        .into()
    })
}

/// Mutable observation state for one preprocessing run.
#[derive(Debug)]
pub struct RunObserver {
    source_label: String,
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
    /// Number of paired record groups dropped because one or both mates were malformed.
    pub invalid_pairs: u64,
    /// Per-reason rejection counts keyed by stable rejection code.
    pub rejection_counts: BTreeMap<&'static str, u64>,
    /// Per-transform counts keyed by stable transform code.
    pub transform_counts: BTreeMap<&'static str, u64>,
    invalid_input_event_counts: BTreeMap<&'static str, u64>,
    invalid_input_warnings_emitted: u64,
    invalid_input_warnings_suppressed: bool,
    invalid_input_samples: Vec<InvalidInputEvent>,
    invalid_input_samples_truncated: bool,
    invalid_input_report: Option<InvalidInputReport>,
}

impl RunObserver {
    /// Start observing a run from one stable source label.
    pub fn new(source_label: String) -> Self {
        Self {
            source_label,
            reads_seen: 0,
            reads_emitted: 0,
            reads_rejected: 0,
            bases_seen: 0,
            bases_emitted: 0,
            pairs_seen: 0,
            pairs_emitted: 0,
            pairs_rejected: 0,
            invalid_reads: 0,
            invalid_pairs: 0,
            rejection_counts: BTreeMap::new(),
            transform_counts: BTreeMap::new(),
            invalid_input_event_counts: BTreeMap::new(),
            invalid_input_warnings_emitted: 0,
            invalid_input_warnings_suppressed: false,
            invalid_input_samples: Vec::new(),
            invalid_input_samples_truncated: false,
            invalid_input_report: None,
        }
    }

    /// Install a JSONL report writer for invalid-input events.
    pub fn set_invalid_input_report(&mut self, report: InvalidInputReport) {
        self.invalid_input_report = Some(report);
    }

    /// Increment counters for one observed record of `bases` length.
    pub fn record_seen(&mut self, bases: usize) {
        self.reads_seen += 1;
        self.bases_seen += bases as u64;
    }

    /// Count one malformed parser slot whose base count is unavailable.
    pub fn record_unparsed_seen(&mut self) {
        self.reads_seen += 1;
        self.invalid_reads += 1;
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

    /// Increment the application count for one transform.
    pub fn record_transform(&mut self, code: &'static str) {
        *self.transform_counts.entry(code).or_default() += 1;
    }

    /// Record one invalid-input event in counters, the detailed report, and bounded samples.
    pub fn record_invalid_input(&mut self, event: InvalidInputEvent) -> Result<()> {
        *self
            .invalid_input_event_counts
            .entry(event.kind())
            .or_default() += 1;

        if let Some(report) = &mut self.invalid_input_report {
            report.write_event(&event)?;
        }

        if self.invalid_input_samples.len() < INVALID_INPUT_SAMPLE_LIMIT {
            self.invalid_input_samples.push(event);
        } else {
            self.invalid_input_samples_truncated = true;
        }

        Ok(())
    }

    /// Record and present a recoverable single-end parser failure.
    pub fn recoverable_single_parser_failure(&mut self, error: &ParseError) -> Result<()> {
        self.record_invalid_input(InvalidInputEvent::SingleRecordParseFailure {
            source: self.source_label.clone(),
            failure: error.into(),
            reads_seen: self.reads_seen,
            continued: true,
        })?;
        self.warn_recoverable_parser_failure(None, error);
        Ok(())
    }

    /// Record and present a recoverable paired parser failure.
    pub fn recoverable_paired_parser_failure(
        &mut self,
        mate: MateSide,
        error: &ParseError,
        continued: bool,
    ) -> Result<()> {
        self.record_invalid_input(InvalidInputEvent::PairedRecordParseFailure {
            source: self.source_label.clone(),
            mate,
            failure: error.into(),
            reads_seen: self.reads_seen,
            pairs_seen: self.pairs_seen,
            continued,
        })?;
        if continued {
            self.warn_recoverable_parser_failure(Some(mate), error);
        }
        Ok(())
    }

    /// Record and present a terminal single-end parser failure.
    pub fn terminal_single_parser_failure(
        &mut self,
        error: &ParseError,
        policy: AdmissionPolicy,
    ) -> Result<()> {
        self.record_invalid_input(InvalidInputEvent::SingleRecordParseFailure {
            source: self.source_label.clone(),
            failure: error.into(),
            reads_seen: self.reads_seen,
            continued: false,
        })?;
        self.warn_terminal_parser_failure(None, error, policy);
        Ok(())
    }

    /// Record and present a terminal paired parser failure.
    pub fn terminal_paired_parser_failure(
        &mut self,
        mate: MateSide,
        error: &ParseError,
        policy: AdmissionPolicy,
    ) -> Result<()> {
        self.record_invalid_input(InvalidInputEvent::PairedRecordParseFailure {
            source: self.source_label.clone(),
            mate,
            failure: error.into(),
            reads_seen: self.reads_seen,
            pairs_seen: self.pairs_seen,
            continued: false,
        })?;
        self.warn_terminal_parser_failure(Some(mate), error, policy);
        Ok(())
    }

    fn warn_recoverable_parser_failure(&mut self, mate: Option<MateSide>, error: &ParseError) {
        if self.should_emit_invalid_input_warning() {
            if let Some(mate) = mate {
                tracing::warn!(
                    source = %self.source_label,
                    mate = %mate,
                    parser_kind = error::constructors::parse_error_kind_name(&error.kind),
                    parser_error = %error,
                    "skipping FASTQ record rejected for unequal sequence and quality lengths"
                );
            } else {
                tracing::warn!(
                    source = %self.source_label,
                    parser_kind = error::constructors::parse_error_kind_name(&error.kind),
                    parser_error = %error,
                    "skipping FASTQ record rejected for unequal sequence and quality lengths"
                );
            }
        } else if self.should_emit_invalid_input_suppressed_notice() {
            tracing::warn!("further invalid-input warnings suppressed");
        }
    }

    fn warn_terminal_parser_failure(
        &self,
        mate: Option<MateSide>,
        error: &ParseError,
        policy: AdmissionPolicy,
    ) {
        if policy == AdmissionPolicy::Skip {
            if let Some(mate) = mate {
                tracing::warn!(
                    source = %self.source_label,
                    mate = %mate,
                    parser_kind = error::constructors::parse_error_kind_name(&error.kind),
                    parser_error = %error,
                    "record parser error is not recoverable; stopping instead of skipping and continuing"
                );
            } else {
                tracing::warn!(
                    source = %self.source_label,
                    parser_kind = error::constructors::parse_error_kind_name(&error.kind),
                    parser_error = %error,
                    "record parser error is not recoverable; stopping instead of skipping and continuing"
                );
            }
        }
    }

    /// Return `true` when another invalid-input warning may be emitted.
    pub fn should_emit_invalid_input_warning(&mut self) -> bool {
        const WARNING_LIMIT: u64 = 5;
        if self.invalid_input_warnings_emitted < WARNING_LIMIT {
            self.invalid_input_warnings_emitted += 1;
            true
        } else {
            false
        }
    }

    /// Return `true` once when invalid-input warning suppression should be announced.
    pub fn should_emit_invalid_input_suppressed_notice(&mut self) -> bool {
        if self.invalid_input_warnings_suppressed {
            false
        } else {
            self.invalid_input_warnings_suppressed = true;
            true
        }
    }

    /// Flush the invalid-input report on successful completion.
    pub fn finish_invalid_input_report(&mut self) -> Result<()> {
        if let Some(report) = &mut self.invalid_input_report {
            report.finish()?;
        }
        Ok(())
    }

    /// Return the bounded invalid-input samples captured for the run summary.
    pub fn invalid_input_samples(&self) -> &[InvalidInputEvent] {
        &self.invalid_input_samples
    }

    /// Return whether invalid-input samples exceeded their in-memory bound.
    pub fn invalid_input_samples_truncated(&self) -> bool {
        self.invalid_input_samples_truncated
    }

    /// Return invalid-input event counts keyed by stable event code.
    pub fn invalid_input_event_counts(&self) -> &BTreeMap<&'static str, u64> {
        &self.invalid_input_event_counts
    }
}
