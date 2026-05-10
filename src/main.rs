//! Streaming FASTQ ingress and preprocessing for `nuclease`.
//!
//! The binary owns local or ENA ingress, emits sequence records in a caller-selected
//! layout/format/encoding, and is intended to compose directly with downstream tools such as
//! `sourmash scripts singlesketch`.

mod adapter;
mod cli;
mod ena;
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
use color_eyre::eyre::Result;

/// Parse CLI arguments, initialize tracing, and execute the selected pipeline.
///
/// # Errors
///
/// Returns an error when ingress selection, reader construction, parsing, transport, or output
/// finalization fails.
fn main() -> Result<()> {
    color_eyre::install()?;
    let cli = cli::Cli::parse();
    cli.init_tracing()?;
    tracing::info!(version = env!("CARGO_PKG_VERSION"), "starting nuclease");
    pipeline::run(&cli)
}
