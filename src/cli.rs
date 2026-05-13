//! Command-line interface and ingress selection types.

use std::{
    fmt,
    io::{self, IsTerminal},
    path::PathBuf,
};

use clap::{
    ArgAction, ArgGroup, Parser, ValueEnum,
    builder::{
        Styles,
        styling::{AnsiColor, Effects},
    },
};
use color_eyre::eyre::{Result, bail, eyre};
use tracing::level_filters::LevelFilter;
use tracing_subscriber::EnvFilter;

use crate::{
    adapter::AdapterPreset,
    ena::Accession,
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
=========================================================

`nuclease` is a fast and resource-frugal sequencing read preprocessor. It streams local or ENA FASTQ through a processing plan and emits cleaned FASTQ/FASTA for downstream tools. It starts with sensible defaults that users can override with additional flags when stricter filtering or transformation is needed. See below for more.
";

const AFTER_HELP: &str = "\
Examples:
  nuclease --in reads.fastq.gz > cleaned.fastq

  nuclease --in1 reads_1.fastq.gz --in2 reads_2.fastq.gz \
    --out1 cleaned_1.fastq.gz --out2 cleaned_2.fastq.gz

  nuclease --in reads.interleaved.fastq.gz --paired --out cleaned.interleaved.fastq.gz

  nuclease --ena SRR35939766 --summary run-summary.json > cleaned.fastq

  nuclease --in reads.fastq.gz --trim-min-q 20 --min-length 75 --min-entropy 1.2 \
    > stricter.fastq

  nuclease --in reads.fastq.gz | sourmash scripts singlesketch - --output reads.sig
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
        .args(["ena", "input", "in1"])
))]
pub struct Cli {
    /// ENA run accession to resolve and stream.
    #[arg(long, help_heading = "Inputs", help = "ENA run accession to stream")]
    pub ena: Option<Accession>,

    /// Local FASTQ input. Treated as single-end unless `--paired` is also passed.
    #[arg(long = "in", help_heading = "Inputs", help = "Local FASTQ input")]
    pub input: Option<PathBuf>,

    /// Treat `--in` as interleaved paired-end FASTQ.
    #[arg(
        long,
        requires = "input",
        help_heading = "Inputs",
        help = "Treat --in as interleaved paired-end FASTQ"
    )]
    pub paired: bool,

    /// Local FASTQ input for split paired read 1.
    #[arg(
        long,
        help_heading = "Inputs",
        help = "Local FASTQ input for split paired read 1"
    )]
    pub in1: Option<PathBuf>,

    /// Local FASTQ input for split paired read 2.
    #[arg(
        long,
        help_heading = "Inputs",
        help = "Local FASTQ input for split paired read 2"
    )]
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

    /// Attempt to merge paired-end reads before per-record filters.
    #[arg(
        long,
        help_heading = "Preprocessing",
        help = "Attempt to merge paired-end reads before trimming and filtering; requires paired input and single-stream output"
    )]
    pub merge_pairs: bool,

    /// Minimum overlap length required before paired reads can merge.
    #[arg(
        long,
        default_value_t = 10,
        help_heading = "Preprocessing",
        help = "Minimum overlap length required before paired reads can merge"
    )]
    pub merge_min_overlap: usize,

    /// Maximum mismatch fraction allowed in candidate paired-read overlaps.
    #[arg(
        long,
        default_value_t = 0.2,
        help_heading = "Preprocessing",
        help = "Maximum mismatch fraction allowed in candidate paired-read overlaps"
    )]
    pub merge_max_mismatch_rate: f32,

    /// Minimum Phred-quality gap required before overlap correction rewrites disagreeing bases.
    #[arg(
        long,
        default_value_t = 0,
        help_heading = "Preprocessing",
        help = "Minimum Phred-quality gap required before overlap correction rewrites disagreeing bases"
    )]
    pub merge_min_correction_delta_q: u8,

    /// How to handle invalid FASTQ input.
    #[arg(
        long,
        value_enum,
        default_value_t = InvalidFastqPolicy::Error,
        help_heading = "Preprocessing",
        help = "How to handle invalid FASTQ input; unrecoverable parser or stream errors are reported and remain fatal"
    )]
    pub invalid_fastq_policy: InvalidFastqPolicy,

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

    /// Write the run's natural single output stream to this path.
    #[arg(
        long,
        conflicts_with_all = ["out1", "out2"],
        help_heading = "Outputs",
        help = "Write the natural single output stream to this file"
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
    /// Local interleaved paired-end FASTQ ingress.
    LocalInterleavedPaired { fastq: PathBuf },
    /// Local split paired-end FASTQ ingress.
    LocalSplitPaired { r1: PathBuf, r2: PathBuf },
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
        if self.ena.is_some() && self.paired {
            bail!(
                "--paired applies only to local --in input; ENA layout is detected automatically"
            );
        }

        if self.paired && self.input.is_none() {
            bail!("--paired requires --in because split --in1/--in2 input is already paired");
        }

        if self.ena.is_some() && (self.input.is_some() || self.in1.is_some() || self.in2.is_some())
        {
            bail!("choose either --ena or local FASTQ input, not both");
        }

        if self.input.is_some() && (self.in1.is_some() || self.in2.is_some()) {
            bail!("--in cannot be combined with --in1/--in2");
        }

        match (&self.ena, &self.input, &self.in1, &self.in2, self.paired) {
            (Some(accession), None, None, None, false) => Ok(Ingress::Ena {
                accession: accession.clone(),
            }),
            (None, Some(input), None, None, true) => Ok(Ingress::LocalInterleavedPaired {
                fastq: input.clone(),
            }),
            (None, Some(input), None, None, false) => Ok(Ingress::LocalSingle {
                fastq: input.clone(),
            }),
            (None, None, Some(in1), Some(in2), _) => Ok(Ingress::LocalSplitPaired {
                r1: in1.clone(),
                r2: in2.clone(),
            }),
            (None, None, Some(_), None, _) => {
                bail!("--in1 requires --in2; use --in for single-end input")
            }
            (None, None, None, Some(_), _) => bail!("--in2 requires --in1"),
            _ => bail!("choose either --ena, --in, or --in1/--in2"),
        }
    }

    /// Convert CLI output flags into raw output arguments for typestate resolution.
    pub fn output_args(&self) -> OutputArgs {
        OutputArgs::new(
            self.output_format,
            self.output_encoding,
            self.out.clone(),
            self.out1.clone(),
            self.out2.clone(),
        )
        .with_merge_pairs(self.merge_pairs)
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
            input: Some(PathBuf::from("reads.fastq.gz")),
            paired: false,
            in1: None,
            in2: None,
            min_length: 50,
            max_ns: 4,
            min_mean_q: 20.0,
            min_entropy: 0.0,
            trim_min_q: 20,
            adapter_preset: AdapterPreset::IlluminaTruSeq,
            merge_pairs: false,
            merge_min_overlap: 10,
            merge_max_mismatch_rate: 0.2,
            merge_min_correction_delta_q: 0,
            invalid_fastq_policy: InvalidFastqPolicy::Error,
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
            input: None,
            paired: false,
            in1: None,
            in2: None,
            min_length: 50,
            max_ns: 4,
            min_mean_q: 20.0,
            trim_min_q: 20,
            min_entropy: 0.0,
            adapter_preset: AdapterPreset::IlluminaTruSeq,
            merge_pairs: false,
            merge_min_overlap: 10,
            merge_max_mismatch_rate: 0.2,
            merge_min_correction_delta_q: 0,
            invalid_fastq_policy: InvalidFastqPolicy::Error,
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
            input: None,
            paired: false,
            in1: Some(PathBuf::from("reads_1.fastq.gz")),
            in2: Some(PathBuf::from("reads_2.fastq.gz")),
            min_length: 50,
            max_ns: 4,
            min_mean_q: 20.0,
            trim_min_q: 20,
            min_entropy: 0.0,
            adapter_preset: AdapterPreset::IlluminaTruSeq,
            merge_pairs: false,
            merge_min_overlap: 10,
            merge_max_mismatch_rate: 0.2,
            merge_min_correction_delta_q: 0,
            invalid_fastq_policy: InvalidFastqPolicy::Error,
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
        assert!(matches!(ingress, Ingress::LocalSplitPaired { .. }));
        Ok(())
    }

    #[test]
    fn ingress_returns_local_single_variant_for_in_without_paired() -> Result<()> {
        let ingress = base_cli().ingress()?;
        assert!(matches!(ingress, Ingress::LocalSingle { .. }));
        Ok(())
    }

    #[test]
    fn ingress_returns_interleaved_paired_variant_for_in_with_paired() -> Result<()> {
        let mut cli = base_cli();
        cli.paired = true;

        let ingress = cli.ingress()?;
        assert!(matches!(ingress, Ingress::LocalInterleavedPaired { .. }));
        Ok(())
    }

    #[test]
    fn ingress_rejects_in1_without_in2() {
        let mut cli = base_cli();
        cli.input = None;
        cli.in1 = Some(PathBuf::from("reads_1.fastq.gz"));

        let error = cli
            .ingress()
            .expect_err("--in1 alone should no longer mean single-end input");
        assert!(error.to_string().contains("--in1 requires --in2"));
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
