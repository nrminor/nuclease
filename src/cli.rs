//! Command-line interface and ingress selection types.

use std::{
    fmt,
    fs::File,
    io::{self, IsTerminal, Read},
    path::PathBuf,
};

use clap::{
    ArgAction, ArgGroup, Parser, ValueEnum,
    builder::{
        Styles,
        styling::{AnsiColor, Effects},
    },
};
use color_eyre::eyre::{Result, WrapErr, bail, eyre};
use tracing::level_filters::LevelFilter;
use tracing_subscriber::EnvFilter;

use crate::{
    adapter::AdapterPreset,
    ena::{Accession, EnaClient, FastqUrlsByLayout},
    output::{OutputArgs, OutputEncoding, OutputFormat},
    progress::ProgressMode,
};

pub const INFO: &str = r"
   ▄     ▄   ▄█▄    █     ▄███▄   ██      ▄▄▄▄▄   ▄███▄   
    █     █  █▀ ▀▄  █     █▀   ▀  █ █    █     ▀▄ █▀   ▀  
██   █ █   █ █   ▀  █     ██▄▄    █▄▄█ ▄  ▀▀▀▀▄   ██▄▄    
█ █  █ █   █ █▄  ▄▀ ███▄  █▄   ▄▀ █  █  ▀▄▄▄▄▀    █▄   ▄▀ 
█  █ █ █▄ ▄█ ▀███▀      ▀ ▀███▀      █            ▀███▀   
█   ██  ▀▀▀                         █                     
                                   ▀                      
==================================================

`nuclease` streams local or ENA FASTQ through a compact preprocessing plan
and emits cleaned FASTQ/FASTA for downstream tools. Defaults are chosen to
be sensible, while additional flags compose into stricter preprocessing when
needed.
";

const AFTER_HELP: &str = "\
Examples:
  nuclease --in1 reads.fastq.gz > cleaned.fastq

  nuclease --in1 reads_1.fastq.gz --in2 reads_2.fastq.gz \
    --out1 cleaned_1.fastq.gz --out2 cleaned_2.fastq.gz

  nuclease --ena SRR35939766 --summary run-summary.json > cleaned.fastq

  nuclease --in1 reads.fastq.gz --trim-min-q 20 --min-length 75 --min-entropy 1.2 \
    > stricter.fastq

  nuclease --in1 reads.fastq.gz | sourmash scripts singlesketch - --output reads.sig
";

pub const STYLES: Styles = Styles::styled()
    .header(AnsiColor::Blue.on_default().effects(Effects::BOLD))
    .usage(AnsiColor::Blue.on_default().effects(Effects::BOLD))
    .literal(AnsiColor::Cyan.on_default().effects(Effects::BOLD))
    .placeholder(AnsiColor::Cyan.on_default())
    .error(AnsiColor::Red.on_default().effects(Effects::BOLD))
    .valid(AnsiColor::Green.on_default().effects(Effects::BOLD))
    .invalid(AnsiColor::Yellow.on_default().effects(Effects::BOLD));

/// Rendering and logging policy derived from `-v`/`-q` flags and terminal capabilities.
pub(crate) struct UiPolicy {
    /// Maximum tracing level to emit, or `None` when tracing should stay disabled.
    pub log_level: Option<LevelFilter>,
    /// Whether to print the final human-readable summary.
    pub show_summary: bool,
    /// Progress rendering mode to use during long-running work.
    pub progress_mode: ProgressMode,
}

/// Policy for handling invalid FASTQ input.
///
/// Some invalid FASTQ events are recoverable because `nuclease` still has a safe record or pair
/// boundary. Parser-level or stream-level corruption is still managed by this policy for reporting
/// and diagnostics, but remains fatal when the stream cannot be safely resynchronized.
#[derive(Clone, Copy, Debug, Eq, PartialEq, ValueEnum)]
pub enum InvalidFastqPolicy {
    /// Fail immediately when invalid FASTQ input is observed.
    Error,
    /// Warn and drop invalid FASTQ records or pairs when recovery is safe.
    WarnDrop,
    /// Drop invalid FASTQ records or pairs without warning when recovery is safe.
    SilentDrop,
}

impl fmt::Display for InvalidFastqPolicy {
    /// Format the policy as its stable machine-readable token.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::Error => "error",
            Self::WarnDrop => "warn_drop",
            Self::SilentDrop => "silent_drop",
        })
    }
}

/// Top-level CLI arguments for `nuclease`.
#[derive(Clone, Debug, Parser)]
#[command(
    name = "nuclease",
    version,
    about = INFO,
    styles = STYLES,
    after_help = AFTER_HELP,
    arg_required_else_help = true,
)]
#[command(group(
    ArgGroup::new("ingress")
        .required(true)
        .args(["ena", "in1"])
))]
pub struct Cli {
    /// ENA run accession to resolve and stream.
    #[arg(long, help_heading = "Inputs", help = "ENA run accession to stream")]
    pub ena: Option<Accession>,

    /// Local FASTQ input for read 1.
    #[arg(long, help_heading = "Inputs", help = "Local FASTQ input for read 1")]
    pub in1: Option<PathBuf>,

    /// Local FASTQ input for read 2.
    #[arg(long, help_heading = "Inputs", help = "Local FASTQ input for read 2")]
    pub in2: Option<PathBuf>,

    /// Minimum post-trim read length.
    #[arg(
        long,
        default_value_t = 50,
        help_heading = "Preprocessing",
        help = "Minimum read length after trimming"
    )]
    pub min_length: usize,

    /// Maximum number of ambiguous `N` bases allowed per record.
    #[arg(
        long,
        default_value_t = 4,
        help_heading = "Preprocessing",
        help = "Maximum number of Ns allowed per read"
    )]
    pub max_ns: usize,

    /// Minimum mean Phred quality allowed per record.
    #[arg(
        long,
        default_value_t = 20.0,
        help_heading = "Preprocessing",
        help = "Minimum mean Phred quality after trimming"
    )]
    pub min_mean_q: f64,

    /// Minimum post-transform Shannon entropy allowed per record.
    #[arg(
        long,
        default_value_t = 0.0,
        help_heading = "Preprocessing",
        help = "Minimum Shannon entropy after trimming"
    )]
    pub min_entropy: f64,

    /// 3' quality trimming cutoff in Phred units.
    #[arg(
        long = "trim-min-q",
        default_value_t = 20,
        help_heading = "Preprocessing",
        help = "3' quality trimming cutoff in Phred units"
    )]
    pub trim_min_q: u8,

    /// Adapter trimming preset to apply.
    #[arg(
        long,
        value_enum,
        default_value_t = AdapterPreset::IlluminaTruSeq,
        help_heading = "Preprocessing",
        help = "Adapter trimming preset to apply"
    )]
    pub adapter_preset: AdapterPreset,

    /// How to handle invalid FASTQ input.
    #[arg(
        long,
        value_enum,
        default_value_t = InvalidFastqPolicy::Error,
        help_heading = "Preprocessing",
        help = "How to handle invalid FASTQ input; unrecoverable parser or stream errors are reported and remain fatal"
    )]
    pub invalid_fastq_policy: InvalidFastqPolicy,

    /// Emit paired records as an interleaved output stream.
    #[arg(
        long,
        help_heading = "Outputs",
        help = "Write paired output as one interleaved stream"
    )]
    pub interleaved: bool,

    /// Output record format.
    #[arg(
        long,
        value_enum,
        default_value_t = OutputFormat::Fastq,
        help_heading = "Outputs",
        help = "Output record format"
    )]
    pub output_format: OutputFormat,

    /// Output encoding.
    #[arg(
        long,
        value_enum,
        help_heading = "Outputs",
        help = "Output encoding (defaults to gzip for .gz outputs, otherwise plain)"
    )]
    pub output_encoding: Option<OutputEncoding>,

    /// Single output path for single-end or interleaved paired output.
    #[arg(
        long,
        conflicts_with_all = ["out1", "out2"],
        help_heading = "Outputs",
        help = "Write output to a single file instead of stdout"
    )]
    pub out: Option<PathBuf>,

    /// First output path for split paired output.
    #[arg(
        long,
        requires = "out2",
        conflicts_with = "out",
        help_heading = "Outputs",
        help = "Write paired read 1 output to this file"
    )]
    pub out1: Option<PathBuf>,

    /// Second output path for split paired output.
    #[arg(
        long,
        requires = "out1",
        conflicts_with = "out",
        help_heading = "Outputs",
        help = "Write paired read 2 output to this file"
    )]
    pub out2: Option<PathBuf>,

    /// Number of records between progress log updates.
    #[arg(
        long,
        default_value_t = 100_000,
        help_heading = "Reporting",
        help = "Reads between progress log updates"
    )]
    pub progress_every: u64,

    /// Write a JSON run summary to this path.
    #[arg(
        long,
        help_heading = "Reporting",
        help = "Write a JSON run summary to this path"
    )]
    pub summary: Option<PathBuf>,

    /// Write every invalid FASTQ event as newline-delimited JSON to this path.
    #[arg(
        long,
        help_heading = "Reporting",
        help = "Write all invalid FASTQ events to this JSONL path"
    )]
    pub invalid_fastq_report: Option<PathBuf>,

    /// Increase tracing verbosity (`-v` = INFO, `-vv` = DEBUG, `-vvv` = TRACE).
    #[arg(short = 'v', long = "verbose", action = ArgAction::Count, help_heading = "Execution")]
    pub verbose: u8,

    /// Reduce logs and reporting (`-q` = ERROR only, `-qq` = no logs or summary, `-qqq` = fully silent).
    #[arg(short = 'q', long = "quiet", action = ArgAction::Count, help_heading = "Execution")]
    pub quiet: u8,
}

/// Abstract ingress choice selected by the user before opening concrete readers.
#[derive(Clone, Debug)]
pub enum Ingress {
    /// ENA-backed ingress driven by a validated run accession.
    Ena { accession: Accession },
    /// Local single-end FASTQ ingress.
    LocalSingle { fastq: PathBuf },
    /// Local paired-end FASTQ ingress.
    LocalPaired { r1: PathBuf, r2: PathBuf },
}

/// Opened ingress readers ready to feed the parser path.
pub enum IngressHandle {
    /// Single-reader ingress.
    Single(Box<dyn Read + Send>),
    /// Paired-reader ingress.
    Paired(Box<dyn Read + Send>, Box<dyn Read + Send>),
}

impl Ingress {
    /// Open concrete readers for the selected ingress.
    ///
    /// # Errors
    ///
    /// Returns an error when local files cannot be opened or ENA lookup/opening fails.
    pub fn open(self) -> Result<IngressHandle> {
        match self {
            Self::LocalSingle { fastq } => {
                let reader = File::open(&fastq).wrap_err_with(|| {
                    format!(
                        "failed to open local single-end FASTQ input: {}\n\
                         help: check that --in1 points to a readable FASTQ file on this machine",
                        fastq.display()
                    )
                })?;
                Ok(IngressHandle::Single(Box::new(reader)))
            }
            Self::LocalPaired { r1, r2 } => Ok(IngressHandle::Paired(
                Box::new(File::open(&r1).wrap_err_with(|| {
                    format!(
                        "failed to open local paired FASTQ input for read 1: {}\n\
                         help: check that --in1 is readable from the current execution environment",
                        r1.display()
                    )
                })?),
                Box::new(File::open(&r2).wrap_err_with(|| {
                    format!(
                        "failed to open local paired FASTQ input for read 2: {}\n\
                         help: check that --in2 is readable from the current execution environment",
                        r2.display()
                    )
                })?),
            )),
            Self::Ena { accession } => {
                let client = EnaClient::new()
                    .wrap_err("failed to construct ENA HTTP client for FASTQ streaming")?;
                match client.lookup_fastq_urls(&accession).wrap_err_with(|| {
                    format!(
                        "failed to resolve ENA FASTQ URLs for accession {accession}\n\
                         help: confirm this is a run accession with public FASTQ files in ENA"
                    )
                })? {
                    FastqUrlsByLayout::Single(url) => {
                        let r = client.open_retrying_stream(url);
                        Ok(IngressHandle::Single(Box::new(r)))
                    }
                    FastqUrlsByLayout::Paired(urls) => {
                        let (r1, r2) = client.open_retrying_paired_streams(urls);
                        Ok(IngressHandle::Paired(Box::new(r1), Box::new(r2)))
                    }
                }
            }
        }
    }
}

impl Cli {
    /// Initialize stderr-backed tracing using the selected verbosity and any `RUST_LOG` override.
    pub fn init_tracing(&self) -> Result<()> {
        let Some(level) = self.ui_policy().log_level else {
            return Ok(());
        };

        let filter = EnvFilter::builder()
            .with_default_directive(level.into())
            .from_env_lossy();

        tracing_subscriber::fmt()
            .with_env_filter(filter)
            .with_writer(io::stderr)
            .try_init()
            .map_err(|error| eyre!("failed to initialize tracing subscriber: {error}"))
    }

    /// Derive logging, summary, and progress behavior from quiet/verbose flags and terminal state.
    pub fn ui_policy(&self) -> UiPolicy {
        let stderr_is_tty = io::stderr().is_terminal();

        let (log_level, show_summary, show_progress) = match self.quiet {
            0 => (
                Some(match self.verbose {
                    0 => LevelFilter::WARN,
                    1 => LevelFilter::INFO,
                    2 => LevelFilter::DEBUG,
                    _ => LevelFilter::TRACE,
                }),
                true,
                true,
            ),
            1 => (Some(LevelFilter::ERROR), true, true),
            2 => (None, false, true),
            _ => (None, false, false),
        };

        let progress_mode = if !show_progress {
            ProgressMode::Off
        } else if stderr_is_tty {
            ProgressMode::Live
        } else {
            ProgressMode::Plain
        };

        UiPolicy {
            log_level,
            show_summary,
            progress_mode,
        }
    }

    /// Resolve mutually exclusive ingress arguments into a typed ingress choice.
    ///
    /// # Errors
    ///
    /// Returns an error when the CLI-selected ingress arguments do not describe a supported
    /// combination.
    pub fn ingress(&self) -> Result<Ingress> {
        match (&self.ena, &self.in1, &self.in2) {
            (Some(accession), None, None) => Ok(Ingress::Ena {
                accession: accession.clone(),
            }),
            (None, Some(in1), Some(in2)) => Ok(Ingress::LocalPaired {
                r1: in1.clone(),
                r2: in2.clone(),
            }),
            (None, Some(in1), None) => Ok(Ingress::LocalSingle { fastq: in1.clone() }),
            (None, None, Some(_)) => bail!("--in2 requires --in1"),
            _ => bail!("choose either --ena or local FASTQ input"),
        }
    }

    /// Convert CLI output flags into raw output arguments for typestate resolution.
    pub fn output_args(&self) -> OutputArgs {
        OutputArgs::new(
            self.interleaved,
            self.output_format,
            self.output_encoding,
            self.out.clone(),
            self.out1.clone(),
            self.out2.clone(),
        )
    }
}

#[cfg(test)]
mod tests {
    use std::path::PathBuf;

    use color_eyre::Result;
    use tracing::level_filters::LevelFilter;

    use super::{Cli, Ingress, InvalidFastqPolicy, UiPolicy};
    use crate::{
        adapter::AdapterPreset,
        ena::Accession,
        output::{OutputEncoding, OutputFormat},
        progress::ProgressMode,
    };

    fn base_cli() -> Cli {
        Cli {
            ena: None,
            in1: Some(PathBuf::from("reads.fastq.gz")),
            in2: None,
            min_length: 50,
            max_ns: 4,
            min_mean_q: 20.0,
            min_entropy: 0.0,
            trim_min_q: 20,
            adapter_preset: AdapterPreset::IlluminaTruSeq,
            invalid_fastq_policy: InvalidFastqPolicy::Error,
            interleaved: false,
            output_format: OutputFormat::Fastq,
            output_encoding: None,
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

    fn live_policy(cli: &Cli) -> UiPolicy {
        let mut policy = cli.ui_policy();
        policy.progress_mode = ProgressMode::Live;
        policy
    }

    #[test]
    fn ingress_returns_ena_variant_for_accession_input() -> Result<()> {
        let cli = Cli {
            ena: Some(Accession::new("SRR35939766")?),
            in1: None,
            in2: None,
            min_length: 50,
            max_ns: 4,
            min_mean_q: 20.0,
            trim_min_q: 20,
            min_entropy: 0.0,
            adapter_preset: AdapterPreset::IlluminaTruSeq,
            invalid_fastq_policy: InvalidFastqPolicy::Error,
            interleaved: false,
            output_format: OutputFormat::Fastq,
            output_encoding: Some(OutputEncoding::Plain),
            out: None,
            out1: None,
            out2: None,
            progress_every: 100_000,
            summary: None,
            invalid_fastq_report: None,
            verbose: 0,
            quiet: 0,
        };

        let ingress = cli.ingress()?;
        assert!(matches!(ingress, Ingress::Ena { .. }));
        Ok(())
    }

    #[test]
    fn ingress_returns_local_paired_variant_when_both_fastqs_are_present() -> Result<()> {
        let cli = Cli {
            ena: None,
            in1: Some(PathBuf::from("reads_1.fastq.gz")),
            in2: Some(PathBuf::from("reads_2.fastq.gz")),
            min_length: 50,
            max_ns: 4,
            min_mean_q: 20.0,
            trim_min_q: 20,
            min_entropy: 0.0,
            adapter_preset: AdapterPreset::IlluminaTruSeq,
            invalid_fastq_policy: InvalidFastqPolicy::Error,
            interleaved: true,
            output_format: OutputFormat::Fastq,
            output_encoding: Some(OutputEncoding::Plain),
            out: None,
            out1: None,
            out2: None,
            progress_every: 100_000,
            summary: None,
            invalid_fastq_report: None,
            verbose: 0,
            quiet: 0,
        };

        let ingress = cli.ingress()?;
        assert!(matches!(ingress, Ingress::LocalPaired { .. }));
        Ok(())
    }

    #[test]
    fn ui_policy_defaults_to_warn_logs_with_live_progress() {
        let cli = base_cli();
        let policy = live_policy(&cli);

        assert_eq!(policy.log_level, Some(LevelFilter::WARN));
        assert!(policy.show_summary);
        assert_eq!(policy.progress_mode, ProgressMode::Live);
    }

    #[test]
    fn ui_policy_maps_verbose_flags_upward() {
        let mut cli = base_cli();
        cli.verbose = 2;
        let policy = live_policy(&cli);

        assert_eq!(policy.log_level, Some(LevelFilter::DEBUG));
        assert!(policy.show_summary);
    }

    #[test]
    fn ui_policy_maps_single_quiet_to_error_only() {
        let mut cli = base_cli();
        cli.quiet = 1;
        let policy = live_policy(&cli);

        assert_eq!(policy.log_level, Some(LevelFilter::ERROR));
        assert!(policy.show_summary);
        assert_eq!(policy.progress_mode, ProgressMode::Live);
    }

    #[test]
    fn ui_policy_maps_double_quiet_to_progress_only() {
        let mut cli = base_cli();
        cli.quiet = 2;
        let policy = live_policy(&cli);

        assert_eq!(policy.log_level, None);
        assert!(!policy.show_summary);
        assert_eq!(policy.progress_mode, ProgressMode::Live);
    }

    #[test]
    fn ui_policy_maps_triple_quiet_to_full_silence() {
        let mut cli = base_cli();
        cli.quiet = 3;
        let policy = cli.ui_policy();

        assert_eq!(policy.log_level, None);
        assert!(!policy.show_summary);
        assert_eq!(policy.progress_mode, ProgressMode::Off);
    }
}
