# nuclease

[![CI](https://github.com/nrminor/nuclease/actions/workflows/ci.yml/badge.svg)](https://github.com/nrminor/nuclease/actions/workflows/ci.yml)
[![Latest release](https://img.shields.io/github/v/release/nrminor/nuclease?label=release)](https://github.com/nrminor/nuclease/releases/latest)

`nuclease` streams local or ENA FASTQ through a compact preprocessing plan and emits cleaned FASTQ or FASTA for downstream tools.

The initial codebase was extracted from the `nt-somm` CLI crate in [`nrminor/nt-terroir`](https://github.com/nrminor/nt-terroir), where its earlier development history lives under `crates/nt-somm-cli`.

## Installation

The recommended installation method is the conda-aware install script:

```sh
curl -fsSL https://raw.githubusercontent.com/nrminor/nuclease/main/INSTALL.sh | bash
```

This downloads a pre-built binary for your platform when one is available. If no binary is available, it falls back to building from source with Cargo. When a conda, mamba, or pixi environment is active, the installer places the binary in that environment's `bin` directory.

## Usage

```sh
nuclease --in1 reads.fastq.gz > cleaned.fastq
nuclease --ena SRR35939766 --summary run-summary.json > cleaned.fastq
```

Run `nuclease --help` for the full CLI.

## Development

```sh
just check
just smoke
```
