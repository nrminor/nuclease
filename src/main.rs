//! Streaming FASTQ ingress and preprocessing for `nuclease`.
//!
//! The binary owns local or ENA ingress, emits sequence records in a caller-selected
//! layout/format/encoding, and is intended to compose directly with downstream tools such as
//! `sourmash scripts singlesketch`.

use std::{
    fmt,
    io::{self, Write as _},
    process::ExitCode,
};

mod adapter;
mod cli;
mod ena;
mod error;
mod filter;
mod output;
mod pair_merge;
mod pipeline;
mod plan;
mod progress;
mod quality;
mod record;
mod report;

use clap::Parser;

/// Parse CLI arguments, initialize tracing, and execute the selected pipeline.
fn main() -> ExitCode {
    if let Err(error) = color_eyre::install() {
        stderr_best_effort(format_args!("failed to install error reporting: {error}"));
        return ExitCode::FAILURE;
    }

    let cli = cli::Cli::parse();
    match pipeline::run(&cli) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            let exit_code = error.exit_code();
            let report = color_eyre::Report::new(error);
            stderr_best_effort(format_args!("{report:?}"));
            exit_code
        }
    }
}

fn stderr_best_effort(arguments: fmt::Arguments<'_>) {
    let _ = writeln!(io::stderr(), "{arguments}");
}
