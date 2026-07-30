//! Typed application errors that determine Nuclease's terminal behavior.

use std::{fmt, io, num::ParseIntError, path::PathBuf, str::Utf8Error};

use needletail::errors::ParseError;
use reqwest::StatusCode;
use thiserror::Error;

use crate::{ena::Accession, record::MateSide};

/// Application result retaining a typed [`RunError`] through the pipeline boundary.
pub(crate) type Result<T> = std::result::Result<T, RunError>;

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

impl From<MalformedInputError> for RunError {
    fn from(error: MalformedInputError) -> Self {
        Self::MalformedInput(Box::new(error))
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
        "FASTQ parser rejected malformed input while reading source={source_label} mate={mate}\n\
         invalid_fastq_policy={policy}\n\
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
        "invalid FASTQ record source={source_label} mate={mate} header={header} sequence_len={sequence_len} quality_len={quality_len}"
    )]
    RecordLength {
        source_label: String,
        mate: &'static str,
        header: String,
        sequence_len: usize,
        quality_len: usize,
    },

    #[error(
        "paired FASTQ headers do not agree source={source_label} left_mate={left_mate} right_mate={right_mate} left_header={left_header} right_header={right_header}\n\
         help: confirm both files contain aligned mates in the same order"
    )]
    MateIdentifier {
        source_label: String,
        left_mate: &'static str,
        right_mate: &'static str,
        left_header: String,
        right_header: String,
    },

    #[error(
        "paired FASTQ inputs have different record counts\n\
         source: {source_label}\n\
         complete_pairs_seen: {complete_pairs_seen}\n\
         reads_seen_before_failure: {reads_seen}\n\
         help: confirm both inputs are complete mates from the same run"
    )]
    PairedRecordCount {
        source_label: String,
        complete_pairs_seen: u64,
        reads_seen: u64,
    },

    #[error(
        "interleaved paired FASTQ ended with an unpaired read\n\
         source: {source_label}\n\
         complete_pairs_seen: {complete_pairs_seen}\n\
         reads_seen_before_failure: {reads_seen}\n\
         help: confirm the interleaved input contains adjacent read pairs and was not truncated"
    )]
    InterleavedRecordCount {
        source_label: String,
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
         reads_seen={reads_seen} pairs_seen={pairs_seen}\n\
         help: confirm the input is FASTQ rather than FASTA and that parser quality computation is enabled"
    )]
    MissingQuality {
        source_label: String,
        mate: &'static str,
        reads_seen: u64,
        pairs_seen: u64,
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
                write!(
                    formatter,
                    "{} mate at {}",
                    mate_label(*mate),
                    path.display()
                )
            }
        }
    }
}

const fn mate_label(mate: MateSide) -> &'static str {
    match mate {
        MateSide::Left => "left",
        MateSide::Right => "right",
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

    #[error(
        "Needletail panicked while parsing unverified ENA input {accession}\n\
         mate: {mate}\n\
         reads_seen: {reads_seen}\n\
         pairs_seen: {pairs_seen}\n\
         panic: {panic}\n\
         help: retry the complete run and report this diagnostic if it is reproducible"
    )]
    ParserPanic {
        accession: Accession,
        mate: &'static str,
        reads_seen: u64,
        pairs_seen: u64,
        panic: String,
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
    #[error("paired identifiers differ: left_header={left_header} right_header={right_header}")]
    MateIdentifier {
        left_header: String,
        right_header: String,
    },
    #[error("FASTQ record for mate={mate} did not include quality scores")]
    MissingQuality { mate: &'static str },
    #[error("paired FASTQ inputs ended after {complete_pairs_seen} complete pairs")]
    RecordCount { complete_pairs_seen: u64 },
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

    #[error(
        "Needletail panicked while parsing local input {source_label}\n\
         mate: {mate}\n\
         reads_seen: {reads_seen}\n\
         pairs_seen: {pairs_seen}\n\
         panic: {panic}\n\
         help: report this diagnostic with the triggering input if possible"
    )]
    LocalParserPanic {
        source_label: String,
        mate: &'static str,
        reads_seen: u64,
        pairs_seen: u64,
        panic: String,
    },
}
