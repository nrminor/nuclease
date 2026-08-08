//! Typed application errors that determine Nuclease's terminal behavior.

use std::{
    fmt::{self, Write as _},
    io,
    num::ParseIntError,
    path::{Path, PathBuf},
    process::ExitCode,
    str::Utf8Error,
};

use needletail::{errors::ParseError, parser::Format};
use reqwest::StatusCode;
use serde::{Serialize, Serializer};
use thiserror::Error;

use crate::{ena::Accession, record::MateSide};

const HEADER_EXCERPT_BYTES: usize = 96;

/// Application result retaining a typed [`RunError`] through the pipeline boundary.
pub(crate) type Result<T> = std::result::Result<T, RunError>;

/// Physical origin used to classify an input failure without carrying parser policy.
#[derive(Clone, Copy)]
pub(crate) enum InputOrigin<'source> {
    Ena(&'source Accession),
    Local(&'source Path),
}

/// Local paired-input layout used only for record-count diagnostics.
#[derive(Clone, Copy)]
pub(crate) enum PairedInputLayout {
    Interleaved,
    Split,
}

/// Central constructors that classify concrete failures into the application error hierarchy.
pub(crate) mod constructors {
    use std::{fmt, path::Path};

    use needletail::{
        errors::{ParseError, ParseErrorKind},
        parser::Format,
    };

    use super::{
        EnaContentProblem, IndeterminateInputError, InputOrigin, InternalError, IoError,
        MalformedInputError, PairConstructionError, PairedInputLayout, ReportWriteError, RunError,
    };
    use crate::{
        observer::RunObserver,
        record::{InputSource, MateSide, RecordView},
    };

    pub(crate) fn missing_single_quality_error(
        origin: InputOrigin<'_>,
        source_label: String,
        observer: &RunObserver,
    ) -> RunError {
        missing_quality_error(origin, source_label, "single", observer)
    }

    pub(crate) fn missing_left_quality_error(
        origin: InputOrigin<'_>,
        source_label: String,
        observer: &RunObserver,
    ) -> RunError {
        missing_quality_error(origin, source_label, "left", observer)
    }

    pub(crate) fn missing_right_quality_error(
        origin: InputOrigin<'_>,
        source_label: String,
        observer: &RunObserver,
    ) -> RunError {
        missing_quality_error(origin, source_label, "right", observer)
    }

    fn missing_quality_error(
        origin: InputOrigin<'_>,
        source_label: String,
        mate: &'static str,
        observer: &RunObserver,
    ) -> RunError {
        match origin {
            InputOrigin::Ena(accession) => IndeterminateInputError::Content {
                accession: accession.to_string(),
                problem: EnaContentProblem::MissingQuality {
                    mate,
                    reads_seen: observer.reads_seen,
                    pairs_seen: observer.pairs_seen,
                    invalid_reads: observer.invalid_reads,
                    invalid_pairs: observer.invalid_pairs,
                },
            }
            .into(),
            InputOrigin::Local(_) => MalformedInputError::MissingQuality {
                source_label,
                mate,
                reads_seen: observer.reads_seen,
                pairs_seen: observer.pairs_seen,
                invalid_reads: observer.invalid_reads,
                invalid_pairs: observer.invalid_pairs,
            }
            .into(),
        }
    }

    pub(crate) fn single_parser_error(
        origin: InputOrigin<'_>,
        source_label: String,
        policy: &dyn fmt::Display,
        observer: &RunObserver,
        source: ParseError,
    ) -> RunError {
        parser_error(
            origin,
            source_label,
            "single",
            None,
            policy,
            observer,
            source,
        )
    }

    pub(crate) fn left_parser_error(
        origin: InputOrigin<'_>,
        source_label: String,
        policy: &dyn fmt::Display,
        observer: &RunObserver,
        source: ParseError,
    ) -> RunError {
        parser_error(
            origin,
            source_label,
            "left",
            Some(MateSide::Left),
            policy,
            observer,
            source,
        )
    }

    pub(crate) fn right_parser_error(
        origin: InputOrigin<'_>,
        source_label: String,
        policy: &dyn fmt::Display,
        observer: &RunObserver,
        source: ParseError,
    ) -> RunError {
        parser_error(
            origin,
            source_label,
            "right",
            Some(MateSide::Right),
            policy,
            observer,
            source,
        )
    }

    fn parser_error(
        origin: InputOrigin<'_>,
        source_label: String,
        mate: &'static str,
        ena_mate: Option<MateSide>,
        policy: &dyn fmt::Display,
        observer: &RunObserver,
        source: ParseError,
    ) -> RunError {
        match origin {
            InputOrigin::Ena(accession) => IndeterminateInputError::Parser {
                accession: accession.clone(),
                mate: ena_mate,
                reads_seen: observer.reads_seen,
                pairs_seen: observer.pairs_seen,
                source,
            }
            .into(),
            InputOrigin::Local(path) if source.kind == ParseErrorKind::Io => {
                IoError::LocalFastqRead {
                    path: path.to_path_buf(),
                    source,
                }
                .into()
            }
            InputOrigin::Local(_)
                if matches!(
                    &source.kind,
                    ParseErrorKind::UnknownFormat | ParseErrorKind::EmptyFile
                ) =>
            {
                MalformedInputError::UnreadableFastq {
                    source_label,
                    mate,
                    reads_seen: observer.reads_seen,
                    pairs_seen: observer.pairs_seen,
                    parser_error_kind: parse_error_kind_name(&source.kind),
                    source,
                }
                .into()
            }
            InputOrigin::Local(_) => MalformedInputError::LocalParser {
                source_label,
                mate,
                policy: policy.to_string(),
                reads_seen: observer.reads_seen,
                pairs_seen: observer.pairs_seen,
                invalid_reads: observer.invalid_reads,
                invalid_pairs: observer.invalid_pairs,
                parser_error_kind: parse_error_kind_name(&source.kind).to_owned(),
                source,
            }
            .into(),
        }
    }

    pub(crate) fn single_unsupported_format_error(
        origin: InputOrigin<'_>,
        source_label: String,
        format: Format,
    ) -> RunError {
        unsupported_format_error(origin, source_label, "single", format)
    }

    pub(crate) fn left_unsupported_format_error(
        origin: InputOrigin<'_>,
        source_label: String,
        format: Format,
    ) -> RunError {
        unsupported_format_error(origin, source_label, "left", format)
    }

    pub(crate) fn right_unsupported_format_error(
        origin: InputOrigin<'_>,
        source_label: String,
        format: Format,
    ) -> RunError {
        unsupported_format_error(origin, source_label, "right", format)
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

    pub(crate) fn left_record_count_error(
        origin: InputOrigin<'_>,
        layout: PairedInputLayout,
        source_label: String,
        header: String,
        observer: &RunObserver,
    ) -> RunError {
        record_count_error(
            origin,
            layout,
            source_label,
            MateSide::Left,
            header,
            observer,
        )
    }

    pub(crate) fn right_record_count_error(
        origin: InputOrigin<'_>,
        layout: PairedInputLayout,
        source_label: String,
        header: String,
        observer: &RunObserver,
    ) -> RunError {
        record_count_error(
            origin,
            layout,
            source_label,
            MateSide::Right,
            header,
            observer,
        )
    }

    fn record_count_error(
        origin: InputOrigin<'_>,
        layout: PairedInputLayout,
        source_label: String,
        present_mate: MateSide,
        header: String,
        observer: &RunObserver,
    ) -> RunError {
        match origin {
            InputOrigin::Ena(accession) => IndeterminateInputError::Content {
                accession: accession.to_string(),
                problem: EnaContentProblem::RecordCount {
                    complete_pairs_seen: observer.pairs_seen,
                    present_mate,
                    header,
                },
            }
            .into(),
            InputOrigin::Local(_) => match layout {
                PairedInputLayout::Interleaved => MalformedInputError::InterleavedRecordCount {
                    source_label,
                    present_mate,
                    header,
                    complete_pairs_seen: observer.pairs_seen,
                    reads_seen: observer.reads_seen,
                }
                .into(),
                PairedInputLayout::Split => MalformedInputError::PairedRecordCount {
                    source_label,
                    present_mate,
                    header,
                    complete_pairs_seen: observer.pairs_seen,
                    reads_seen: observer.reads_seen,
                }
                .into(),
            },
        }
    }

    pub(crate) fn pair_construction_error(
        origin: InputOrigin<'_>,
        source_label: String,
        source: PairConstructionError,
    ) -> RunError {
        match origin {
            InputOrigin::Ena(accession) => IndeterminateInputError::Content {
                accession: accession.to_string(),
                problem: EnaContentProblem::PairConstruction {
                    source: Box::new(source),
                },
            }
            .into(),
            InputOrigin::Local(_) => MalformedInputError::PairConstruction {
                source_label,
                source,
            }
            .into(),
        }
    }

    pub(crate) fn record_length_error(record: RecordView<'_>, mate: Option<MateSide>) -> RunError {
        let mate = match mate {
            Some(MateSide::Left) => "left",
            Some(MateSide::Right) => "right",
            None => "single",
        };
        match record.source() {
            Some(InputSource::Ena { accession }) => IndeterminateInputError::Content {
                accession: accession.to_owned(),
                problem: EnaContentProblem::RecordLength {
                    mate,
                    header: String::from_utf8_lossy(record.header()).into_owned(),
                    sequence_len: record.sequence().len(),
                    quality_len: record.quality().len(),
                },
            }
            .into(),
            _ => MalformedInputError::RecordLength {
                source_label: record.source_display(),
                mate,
                header: String::from_utf8_lossy(record.header()).into_owned(),
                sequence_len: record.sequence().len(),
                quality_len: record.quality().len(),
            }
            .into(),
        }
    }

    pub(crate) fn invalid_utf8_error(
        record: RecordView<'_>,
        field: &'static str,
        source: std::str::Utf8Error,
    ) -> RunError {
        match record.source() {
            Some(InputSource::Ena { accession }) => IndeterminateInputError::InvalidUtf8 {
                accession: accession.to_owned(),
                field,
                source,
            }
            .into(),
            _ => MalformedInputError::InvalidUtf8 { field, source }.into(),
        }
    }

    pub(crate) fn json_report_error(
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

    pub(crate) const fn parse_error_kind_name(kind: &ParseErrorKind) -> &'static str {
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
}

/// Terminal Nuclease failures grouped by stable process-response category.
#[derive(Debug, Error)]
pub(crate) enum RunError {
    #[error(transparent)]
    Usage(#[from] UsageError),
    #[error(transparent)]
    MalformedInput(Box<MalformedInputError>),
    #[error(transparent)]
    UnavailableInput(#[from] UnavailableInputError),
    #[error(transparent)]
    Io(#[from] IoError),
    #[error(transparent)]
    IndeterminateInput(#[from] IndeterminateInputError),
    #[error(transparent)]
    Internal(#[from] InternalError),
}

impl RunError {
    /// Return the stable process status for this terminal failure category.
    pub(crate) fn exit_code(&self) -> ExitCode {
        match self {
            Self::Usage(_) => ExitCode::from(2),
            Self::MalformedInput(_) => ExitCode::from(65),
            Self::UnavailableInput(_) => ExitCode::from(66),
            Self::Io(_) => ExitCode::from(74),
            Self::IndeterminateInput(_) => ExitCode::from(75),
            Self::Internal(_) => ExitCode::FAILURE,
        }
    }
}

impl From<MalformedInputError> for RunError {
    fn from(error: MalformedInputError) -> Self {
        Self::MalformedInput(Box::new(error))
    }
}

/// Bounded, allocation-free record-header evidence retained by typed errors.
#[derive(Clone, Copy, Eq, PartialEq)]
pub(crate) struct HeaderExcerpt {
    bytes: [u8; HEADER_EXCERPT_BYTES],
    len: u8,
    truncated: bool,
}

impl HeaderExcerpt {
    /// Capture a bounded prefix of one raw record header without allocating.
    pub(crate) fn capture(header: &[u8]) -> Self {
        let copied = header.len().min(HEADER_EXCERPT_BYTES);
        let mut bytes = [0; HEADER_EXCERPT_BYTES];
        bytes[..copied].copy_from_slice(&header[..copied]);
        Self {
            bytes,
            len: u8::try_from(copied).expect("header excerpt bound should fit in u8"),
            truncated: header.len() > copied,
        }
    }
}

impl fmt::Display for HeaderExcerpt {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        for byte in &self.bytes[..usize::from(self.len)] {
            for escaped in std::ascii::escape_default(*byte) {
                formatter.write_char(char::from(escaped))?;
            }
        }
        if self.truncated {
            formatter.write_str("…")?;
        }
        Ok(())
    }
}

impl fmt::Debug for HeaderExcerpt {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "\"{self}\"")
    }
}

impl Serialize for HeaderExcerpt {
    fn serialize<SerializerType>(
        &self,
        serializer: SerializerType,
    ) -> std::result::Result<SerializerType::Ok, SerializerType::Error>
    where
        SerializerType: Serializer,
    {
        serializer.collect_str(self)
    }
}

/// Pair-domain contradictions discovered while establishing two positional mate claims.
#[derive(Clone, Copy, Debug, Eq, Error, PartialEq, Serialize)]
#[serde(tag = "reason", rename_all = "snake_case")]
pub(crate) enum PairConstructionError {
    #[error(
        "left record has unequal sequence and quality lengths header={header} sequence_len={sequence_len} quality_len={quality_len}"
    )]
    LeftRecordLength {
        header: HeaderExcerpt,
        sequence_len: usize,
        quality_len: usize,
    },
    #[error(
        "right record has unequal sequence and quality lengths header={header} sequence_len={sequence_len} quality_len={quality_len}"
    )]
    RightRecordLength {
        header: HeaderExcerpt,
        sequence_len: usize,
        quality_len: usize,
    },
    #[error(
        "both records have unequal sequence and quality lengths left_header={left_header} left_sequence_len={left_sequence_len} left_quality_len={left_quality_len} right_header={right_header} right_sequence_len={right_sequence_len} right_quality_len={right_quality_len}"
    )]
    BothRecordLengths {
        left_header: HeaderExcerpt,
        left_sequence_len: usize,
        left_quality_len: usize,
        right_header: HeaderExcerpt,
        right_sequence_len: usize,
        right_quality_len: usize,
    },
    #[error("left positional header explicitly claims the right mate header={header}")]
    LeftHeaderClaimsRight { header: HeaderExcerpt },
    #[error("right positional header explicitly claims the left mate header={header}")]
    RightHeaderClaimsLeft { header: HeaderExcerpt },
    #[error(
        "both positional headers contradict their mate positions left_header={left_header} right_header={right_header}"
    )]
    BothHeadersContradictPositions {
        left_header: HeaderExcerpt,
        right_header: HeaderExcerpt,
    },
    #[error(
        "paired record identifiers do not agree left_header={left_header} right_header={right_header}"
    )]
    IdentifierMismatch {
        left_header: HeaderExcerpt,
        right_header: HeaderExcerpt,
    },
}

impl PairConstructionError {
    /// Return the stable admission-report category for this construction failure.
    pub(crate) const fn code(&self) -> &'static str {
        match self {
            Self::LeftRecordLength { .. }
            | Self::RightRecordLength { .. }
            | Self::BothRecordLengths { .. } => "sequence_quality_length_mismatch",
            Self::LeftHeaderClaimsRight { .. }
            | Self::RightHeaderClaimsLeft { .. }
            | Self::BothHeadersContradictPositions { .. } => "pair_mate_designation_mismatch",
            Self::IdentifierMismatch { .. } => "pair_identifier_mismatch",
        }
    }
}

/// Error returned when an ENA accession does not match the CLI value grammar.
#[derive(Debug, Error)]
#[error(
    "invalid ENA run accession {value}: {problem}\n\
     help: use an SRR, ERR, or DRR run accession followed by digits"
)]
pub struct AccessionParseError {
    value: String,
    problem: AccessionSyntaxProblem,
}

impl AccessionParseError {
    pub(crate) fn too_short(value: String) -> Self {
        Self {
            value,
            problem: AccessionSyntaxProblem::TooShort,
        }
    }

    pub(crate) fn unsupported_prefix(value: String) -> Self {
        Self {
            value,
            problem: AccessionSyntaxProblem::UnsupportedPrefix,
        }
    }

    pub(crate) fn non_numeric_suffix(value: String) -> Self {
        Self {
            value,
            problem: AccessionSyntaxProblem::NonNumericSuffix,
        }
    }
}

#[derive(Debug, Error)]
enum AccessionSyntaxProblem {
    #[error("the value is too short")]
    TooShort,
    #[error("the prefix is not SRR, ERR, or DRR")]
    UnsupportedPrefix,
    #[error("the suffix is not numeric")]
    NonNumericSuffix,
}

/// Semantically invalid selections discovered after successful CLI parsing or input resolution.
#[derive(Debug, Error)]
pub(crate) enum UsageError {
    #[error(
        "paired output paths imply different encodings\n\
         out1: {out1} ({encoding1})\n\
         out2: {out2} ({encoding2})\n\
         help: use matching suffixes or specify --output-encoding explicitly"
    )]
    PairedEncodingMismatch {
        out1: PathBuf,
        encoding1: &'static str,
        out2: PathBuf,
        encoding2: &'static str,
    },

    #[error(
        "single-end output accepts either stdout or --out only\n\
         help: remove --out1/--out2 for single-end input"
    )]
    SingleOutputDestination,

    #[error(
        "--merge-pairs requires paired-end input\n\
         help: provide --in1 and --in2, use --in with --paired, or use an ENA accession with paired FASTQ files"
    )]
    MergeRequiresPairedInput,
}

/// Local input failures whose structure is demonstrably malformed.
#[derive(Debug, Error)]
pub(crate) enum MalformedInputError {
    #[error(
        "record parser rejected malformed input while reading source={source_label} mate={mate}\n\
         admission_policy={policy}\n\
         reads_seen={reads_seen} pairs_seen={pairs_seen} invalid_reads={invalid_reads} invalid_pairs={invalid_pairs}\n\
         parser_error_kind={parser_error_kind}\n\
         help: inspect the local FASTQ structure and compression before retrying"
    )]
    LocalParser {
        source_label: String,
        mate: &'static str,
        policy: String,
        reads_seen: u64,
        pairs_seen: u64,
        invalid_reads: u64,
        invalid_pairs: u64,
        parser_error_kind: String,
        #[source]
        source: ParseError,
    },

    #[error(
        "input source did not provide a readable FASTQ stream source={source_label} mate={mate}\n\
         reads_seen={reads_seen} pairs_seen={pairs_seen}\n\
         parser_error_kind={parser_error_kind}\n\
         help: confirm the local input is a non-empty FASTQ stream and that compression was detected correctly"
    )]
    UnreadableFastq {
        source_label: String,
        mate: &'static str,
        reads_seen: u64,
        pairs_seen: u64,
        parser_error_kind: &'static str,
        #[source]
        source: ParseError,
    },

    #[error(
        "unsupported input format source={source_label} mate={mate} format={format:?}\n\
         help: this release accepts FASTQ input; use --output-format fasta only to convert admitted FASTQ records"
    )]
    UnsupportedFormat {
        source_label: String,
        mate: &'static str,
        format: Format,
    },

    #[error(
        "record failed admission source={source_label} mate={mate} header={header} sequence_len={sequence_len} quality_len={quality_len}"
    )]
    RecordLength {
        source_label: String,
        mate: &'static str,
        header: String,
        sequence_len: usize,
        quality_len: usize,
    },

    #[error(
        "paired input could not be established source={source_label}: {source}\n\
         help: confirm both inputs contain aligned mates in left/right order"
    )]
    PairConstruction {
        source_label: String,
        #[source]
        source: PairConstructionError,
    },

    #[error(
        "paired FASTQ inputs have different record counts\n\
         source: {source_label}\n\
         present_mate: {present_mate}\n\
         header: {header}\n\
         complete_pairs_seen: {complete_pairs_seen}\n\
         reads_seen_before_failure: {reads_seen}\n\
         help: confirm both inputs are complete mates from the same run"
    )]
    PairedRecordCount {
        source_label: String,
        present_mate: MateSide,
        header: String,
        complete_pairs_seen: u64,
        reads_seen: u64,
    },

    #[error(
        "interleaved paired FASTQ ended with an unpaired read\n\
         source: {source_label}\n\
         present_mate: {present_mate}\n\
         header: {header}\n\
         complete_pairs_seen: {complete_pairs_seen}\n\
         reads_seen_before_failure: {reads_seen}\n\
         help: confirm the interleaved input contains adjacent read pairs and was not truncated"
    )]
    InterleavedRecordCount {
        source_label: String,
        present_mate: MateSide,
        header: String,
        complete_pairs_seen: u64,
        reads_seen: u64,
    },

    #[error(
        "{field} must be UTF-8 for paired-read assembly\n\
         help: inspect the affected local FASTQ record or disable --merge-pairs"
    )]
    InvalidUtf8 {
        field: &'static str,
        #[source]
        source: Utf8Error,
    },

    #[error(
        "FASTQ parser did not provide quality scores while reading source={source_label} mate={mate}\n\
         reads_seen={reads_seen} pairs_seen={pairs_seen} invalid_reads={invalid_reads} invalid_pairs={invalid_pairs}\n\
         help: confirm the input is FASTQ rather than FASTA and that parser quality computation is enabled"
    )]
    MissingQuality {
        source_label: String,
        mate: &'static str,
        reads_seen: u64,
        pairs_seen: u64,
        invalid_reads: u64,
        invalid_pairs: u64,
    },
}

/// Requested inputs that cannot be opened or resolved in a supported form.
#[derive(Debug, Error)]
pub(crate) enum UnavailableInputError {
    #[error(
        "failed to open local FASTQ input {path}\n\
         help: check that the path exists and is readable from this environment"
    )]
    OpenLocalFastq {
        path: PathBuf,
        #[source]
        source: io::Error,
    },

    #[error(
        "ENA metadata for accession {accession} supplied no supported FASTQ input: {problem}\n\
         help: confirm the accession resolves to a supported run-level FASTQ layout"
    )]
    EnaMetadata {
        accession: Accession,
        problem: EnaAvailabilityProblem,
    },
}

/// Authoritative ENA metadata outcomes that supply no supported input.
#[derive(Debug, Error)]
pub(crate) enum EnaAvailabilityProblem {
    #[error("no data row")]
    NoDataRow,
    #[error("the FASTQ location is empty")]
    EmptyFastqLocation,
    #[error("the FASTQ layout or cardinality is unsupported")]
    UnsupportedLayout,
}

/// Resolved destination of one required scientific output stream.
#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) enum OutputDestination {
    Stdout,
    File(PathBuf),
    MateFile { mate: MateSide, path: PathBuf },
}

impl fmt::Display for OutputDestination {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Stdout => formatter.write_str("stdout"),
            Self::File(path) => write!(formatter, "{}", path.display()),
            Self::MateFile { mate, path } => {
                write!(formatter, "{mate} mate at {}", path.display())
            }
        }
    }
}

/// External local read, required-output, and requested-report failures.
#[derive(Debug, Error)]
pub(crate) enum IoError {
    #[error(
        "failed while reading local FASTQ input {path}\n\
         help: check filesystem and storage health, then retry the read"
    )]
    LocalFastqRead {
        path: PathBuf,
        #[source]
        source: ParseError,
    },

    #[error(
        "failed to create required output {destination}\n\
         encoding: {encoding}\n\
         help: check that the parent directory exists and is writable"
    )]
    CreateOutput {
        destination: OutputDestination,
        encoding: &'static str,
        #[source]
        source: io::Error,
    },

    #[error(
        "required output {destination} was closed before nuclease finished writing\n\
         help: ensure the downstream process consumes the complete output"
    )]
    BrokenPipe {
        destination: OutputDestination,
        #[source]
        source: io::Error,
    },

    #[error(
        "failed to write required output {destination}\n\
         help: check the output filesystem or downstream process"
    )]
    WriteOutput {
        destination: OutputDestination,
        #[source]
        source: io::Error,
    },

    #[error(
        "failed to finalize required output {destination}\n\
         help: check the output filesystem or downstream process"
    )]
    FinalizeOutput {
        destination: OutputDestination,
        #[source]
        source: io::Error,
    },

    #[error(
        "failed to create requested {report_kind} report at {path}\n\
         help: check that the parent directory exists and is writable"
    )]
    CreateReport {
        report_kind: &'static str,
        path: PathBuf,
        #[source]
        source: io::Error,
    },

    #[error(
        "failed to write requested {report_kind} report at {path}\n\
         help: check that the destination is writable and has sufficient space"
    )]
    WriteReport {
        report_kind: &'static str,
        path: PathBuf,
        #[source]
        source: ReportWriteError,
    },

    #[error(
        "failed to finalize requested {report_kind} report at {path}\n\
         help: check that the destination is writable and has sufficient space"
    )]
    FinalizeReport {
        report_kind: &'static str,
        path: PathBuf,
        #[source]
        source: io::Error,
    },
}

/// Concrete mechanism failure encountered while writing a requested report.
#[derive(Debug, Error)]
pub(crate) enum ReportWriteError {
    #[error(transparent)]
    Json(serde_json::Error),
    #[error(transparent)]
    Bytes(io::Error),
}

/// ENA failures for which this run could not establish trustworthy input evidence.
#[derive(Debug, Error)]
pub(crate) enum IndeterminateInputError {
    #[error(
        "ENA metadata request failed for accession {accession}\n\
         help: check network access to ebi.ac.uk and consider a bounded whole-run retry"
    )]
    MetadataTransport {
        accession: Accession,
        #[source]
        source: reqwest::Error,
    },

    #[error(
        "ENA returned HTTP status {status} for accession {accession}\n\
         help: consider a bounded whole-run retry and retain diagnostics if the response recurs"
    )]
    Status {
        accession: Accession,
        status: StatusCode,
    },

    #[error(
        "ENA metadata for accession {accession} could not be trusted: {problem}\n\
         help: consider a bounded whole-run retry and retain the response details if it recurs"
    )]
    Metadata {
        accession: Accession,
        #[source]
        problem: EnaMetadataProblem,
    },

    #[error(
        "failed to open and validate ENA FASTQ stream\n\
         accession: {accession}\n\
         mate: {mate:?}\n\
         help: consider a bounded retry of the complete ENA-backed run"
    )]
    Stream {
        accession: Accession,
        mate: Option<MateSide>,
        #[source]
        source: io::Error,
    },

    #[error(
        "FASTQ parsing stopped before ENA input could be trusted\n\
         accession: {accession}\n\
         mate: {mate:?}\n\
         reads_seen: {reads_seen}\n\
         pairs_seen: {pairs_seen}\n\
         help: consider a bounded retry of the complete ENA-backed run"
    )]
    Parser {
        accession: Accession,
        mate: Option<MateSide>,
        reads_seen: u64,
        pairs_seen: u64,
        #[source]
        source: ParseError,
    },

    #[error(
        "ENA FASTQ content for accession {accession} could not be trusted: {problem}\n\
         help: consider a bounded retry of the complete ENA-backed run"
    )]
    Content {
        accession: String,
        #[source]
        problem: EnaContentProblem,
    },

    #[error(
        "ENA FASTQ content for accession {accession} could not be trusted: {field} was not valid UTF-8\n\
         help: consider a bounded retry of the complete ENA-backed run"
    )]
    InvalidUtf8 {
        accession: String,
        field: &'static str,
        #[source]
        source: Utf8Error,
    },
}

/// Closed set of malformed or contradictory ENA file-report states.
#[derive(Debug, Error)]
pub(crate) enum EnaMetadataProblem {
    #[error("the file-report response was empty")]
    EmptyResponse,
    #[error("the file-report returned more than one run row")]
    MultipleRows,
    #[error("the file-report row shape did not match its header")]
    RowShape,
    #[error("the file-report did not include required field {field}")]
    MissingField { field: &'static str },
    #[error("required file-report field {field} was empty")]
    EmptyField { field: &'static str },
    #[error("the file-report returned accession {returned} instead of the requested accession")]
    AccessionMismatch { returned: String },
    #[error("FASTQ URL, byte-count, and MD5 cardinalities differ")]
    CardinalityMismatch,
    #[error("paired FASTQ URLs were not distinct")]
    DuplicatePairedUrls,
    #[error("FASTQ URL used a non-HTTPS scheme: {value}")]
    NonHttpsUrl { value: String },
    #[error("FASTQ URL path did not end with .fastq.gz: {value}")]
    UnexpectedUrlSuffix { value: String },
    #[error("fastq_md5 did not contain 32 ASCII hexadecimal characters: {value}")]
    InvalidMd5Shape { value: String },
    #[error("fastq_bytes value was not numeric: {value}")]
    InvalidByteCount {
        value: String,
        #[source]
        source: ParseIntError,
    },
    #[error("FASTQ URL was invalid: {value}")]
    InvalidUrl {
        value: String,
        #[source]
        source: url::ParseError,
    },
    #[error("fastq_md5 value was not hexadecimal: {value}")]
    InvalidMd5 {
        value: String,
        #[source]
        source: ParseIntError,
    },
}

/// Closed set of record-level defects observed in untrusted ENA content.
#[derive(Debug, Error)]
pub(crate) enum EnaContentProblem {
    #[error(
        "sequence and quality lengths differ for mate={mate} header={header} sequence_len={sequence_len} quality_len={quality_len}"
    )]
    RecordLength {
        mate: &'static str,
        header: String,
        sequence_len: usize,
        quality_len: usize,
    },
    #[error("paired-record construction failed: {source}")]
    PairConstruction {
        #[source]
        source: Box<PairConstructionError>,
    },
    #[error(
        "FASTQ record for mate={mate} did not include quality scores after reads_seen={reads_seen} pairs_seen={pairs_seen} invalid_reads={invalid_reads} invalid_pairs={invalid_pairs}"
    )]
    MissingQuality {
        mate: &'static str,
        reads_seen: u64,
        pairs_seen: u64,
        invalid_reads: u64,
        invalid_pairs: u64,
    },
    #[error("input provided unsupported format {format:?} for mate={mate}")]
    UnsupportedFormat { mate: &'static str, format: Format },
    #[error(
        "paired FASTQ inputs ended after {complete_pairs_seen} complete pairs with unmatched {present_mate} record {header}"
    )]
    RecordCount {
        complete_pairs_seen: u64,
        present_mate: MateSide,
        header: String,
    },
}

/// Software defects, initialization failures, and violated internal expectations.
#[derive(Debug, Error)]
pub(crate) enum InternalError {
    #[error(
        "failed to initialize tracing subscriber\n\
         help: report this diagnostic if it persists in a clean environment"
    )]
    Tracing {
        #[source]
        source: Box<dyn std::error::Error + Send + Sync>,
    },

    #[error(
        "failed to initialize the ENA HTTP client\n\
         help: report this diagnostic if it persists in a clean environment"
    )]
    HttpClientInitialization {
        #[source]
        source: reqwest::Error,
    },

    #[error(
        "paired-read assembly failed unexpectedly\n\
         help: report this diagnostic with the triggering input and options"
    )]
    PairAssembly {
        #[source]
        source: libpairassembly::Error,
    },

    #[error(
        "Nuclease reached an impossible parsed CLI state: {detail}\n\
         help: report this as a Nuclease bug"
    )]
    CliInvariant { detail: String },

    #[error(
        "Nuclease reached an unsupported execution-plan state: {detail}\n\
         help: report this as a Nuclease bug"
    )]
    PlanInvariant { detail: String },

    #[error(
        "replacement sequence and quality lengths differ header={header} sequence_len={sequence_len} quality_len={quality_len}\n\
         help: report this as a Nuclease bug"
    )]
    ReplacementLength {
        header: String,
        sequence_len: usize,
        quality_len: usize,
    },

    #[error(
        "FASTQ output requires quality scores\n\
         help: report this as a Nuclease bug"
    )]
    MissingOutputQuality,

    #[error(
        "failed to serialize internally constructed {report_kind} report\n\
         help: report this as a Nuclease bug"
    )]
    SerializeReport {
        report_kind: &'static str,
        #[source]
        source: serde_json::Error,
    },
}

#[cfg(test)]
mod tests {
    use std::{error::Error as _, io, path::PathBuf, process::ExitCode};

    use super::{
        EnaAvailabilityProblem, HeaderExcerpt, IndeterminateInputError, InputOrigin, InternalError,
        IoError, MalformedInputError, OutputDestination, PairConstructionError, RunError,
        UnavailableInputError, UsageError, constructors,
    };
    use crate::ena::Accession;

    #[test]
    fn run_error_categories_map_to_stable_exit_codes() {
        let accession = Accession::new("SRR35939766").expect("test accession should be valid");
        let cases = [
            (RunError::from(UsageError::MergeRequiresPairedInput), 2),
            (
                RunError::from(MalformedInputError::MissingQuality {
                    source_label: "local:reads.fastq".to_owned(),
                    mate: "single",
                    reads_seen: 0,
                    pairs_seen: 0,
                    invalid_reads: 0,
                    invalid_pairs: 0,
                }),
                65,
            ),
            (
                RunError::from(UnavailableInputError::EnaMetadata {
                    accession: accession.clone(),
                    problem: EnaAvailabilityProblem::NoDataRow,
                }),
                66,
            ),
            (
                RunError::from(IoError::WriteOutput {
                    destination: OutputDestination::File(PathBuf::from("results.fastq")),
                    source: io::Error::other("test output failure"),
                }),
                74,
            ),
            (
                RunError::from(IndeterminateInputError::Status {
                    accession,
                    status: reqwest::StatusCode::SERVICE_UNAVAILABLE,
                }),
                75,
            ),
            (
                RunError::from(InternalError::CliInvariant {
                    detail: "test invariant".to_owned(),
                }),
                1,
            ),
        ];

        for (error, expected) in cases {
            assert_eq!(error.exit_code(), ExitCode::from(expected));
        }
    }

    fn identifier_mismatch() -> PairConstructionError {
        PairConstructionError::IdentifierMismatch {
            left_header: HeaderExcerpt::capture(b"left/1"),
            right_header: HeaderExcerpt::capture(b"right/2"),
        }
    }

    #[test]
    fn local_pair_construction_error_preserves_typed_source_chain() {
        let error = constructors::pair_construction_error(
            InputOrigin::Local(std::path::Path::new("reads.fastq")),
            "local:reads.fastq".to_owned(),
            identifier_mismatch(),
        );

        assert_eq!(error.exit_code(), ExitCode::from(65));
        let RunError::MalformedInput(error) = &error else {
            panic!("local pair failure should be malformed input");
        };
        assert!(matches!(
            &**error,
            MalformedInputError::PairConstruction {
                source: PairConstructionError::IdentifierMismatch { .. },
                ..
            }
        ));
        assert!(error.source().is_some());
    }

    #[test]
    fn ena_pair_construction_error_preserves_typed_source_chain() {
        let accession = Accession::new("SRR35939766").expect("test accession should be valid");
        let error = constructors::pair_construction_error(
            InputOrigin::Ena(&accession),
            format!("ena:{accession}"),
            identifier_mismatch(),
        );

        assert_eq!(error.exit_code(), ExitCode::from(75));
        let RunError::IndeterminateInput(IndeterminateInputError::Content { problem, .. }) = &error
        else {
            panic!("ENA pair failure should be indeterminate input");
        };
        assert!(matches!(
            problem,
            super::EnaContentProblem::PairConstruction { source }
                if matches!(
                    &**source,
                    PairConstructionError::IdentifierMismatch { .. }
                )
        ));
        assert!(problem.source().is_some());
    }

    #[test]
    fn header_excerpt_escapes_and_bounds_arbitrary_bytes() {
        let mut header = vec![b'A'; 120];
        header[1] = 0xff;
        header[2] = b'\n';
        header[3] = 0x01;

        let rendered = HeaderExcerpt::capture(&header).to_string();

        assert!(rendered.starts_with("A\\xff\\n\\x01"));
        assert!(rendered.ends_with('…'));
        assert!(!rendered.contains('\n'));
        assert!(rendered.len() < header.len() + 16);
    }
}
