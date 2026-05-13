# nuclease - a streaming FASTQ preprocessor

[![CI](https://github.com/nrminor/nuclease/actions/workflows/ci.yml/badge.svg)](https://github.com/nrminor/nuclease/actions/workflows/ci.yml)
[![Latest release](https://img.shields.io/github/v/release/nrminor/nuclease?label=release)](https://github.com/nrminor/nuclease/releases/latest)

## Overview

`nuclease` is a preprocessor for large sequencing read datasets, akin to [fastp](https://github.com/OpenGene/fastp) or [VSEARCH](https://github.com/torognes/vsearch) (though for now with fewer features). It was initially created to allow reads to be streamed directly from ENA and trimmed _in-process_, meaning that preprocessed reads can be streamed directly from ENA into whatever tool. For large projects with lots of big FASTQs, this means fewer temporary FASTQs on disk and thus smaller Nextflow work directories.

Here's what I mean:

```bash
nuclease \
--ena SRR35946892 \
--interleaved-out \
--min-length 100 \
--max-ns 0 \
--min-mean-q 20 \
--min-entropy 0.7 \
--trim-min-q 10 \
--adapter-preset illumina-truseq \
--summary summary.json \
| sourmash sketch dna -p k=31 - -o sample.zip
```

Here, `nuclease` is used to stream and preprocess reads from ENA and then pass them along into kmer sketching with Sourmash. At no point in this preprocessing will the full read dataset be materialized on disk. Instead, reads are streamed and operated on as they arrive from ENA.

`nuclease`'s second big feature is developer-focused: it is highly extensible. See the [Contributing to nuclease](#contributing-to-nuclease) section for more.

## Installation

A simple curl install is the recommended installation method:

```sh
curl -fsSL https://raw.githubusercontent.com/nrminor/nuclease/main/INSTALL.sh | bash
```

This downloads a pre-built binary for your platform when one is available. If no binary is available, it falls back to building from source with Cargo. When a conda, mamba, or pixi environment is active, the installer places the binary in that environment's `bin` directory.

## Supported Preprocessing Operations

`nuclease` currently has a small built-in set of filters and transforms. Filters decide whether a read should continue through the stream. Transforms may change the sequence and quality data before later filters are evaluated. Most operations act on one record at a time, but paired-end input can also run pair-aware preprocessing before the ordinary per-record steps.

Supported filters are:

- Maximum ambiguous bases (`--max-ns`): rejects reads with more than the configured number of `N` bases.
- Minimum read length (`--min-length`): rejects post-transform reads shorter than the configured length.
- Minimum mean Phred quality (`--min-mean-q`): rejects post-transform reads whose mean quality is below the configured threshold.
- Minimum Shannon entropy (`--min-entropy`): rejects post-transform reads whose `A`/`C`/`G`/`T` composition is below the configured entropy threshold.

Supported transforms are:

- Paired-read overlap merging (`--merge-pairs`): attempts to merge paired-end reads before per-record trimming and filtering. Merging is powered by [`libpairassembly`](https://github.com/nrminor/pairassembler), and can be tuned with `--merge-min-overlap`, `--merge-max-mismatch-rate`, and `--merge-min-correction-delta-q`.
- Illumina TruSeq adapter trimming (`--adapter-preset illumina-truseq`): trims detected 3' adapter overlap using the built-in TruSeq R1/R2 adapter catalog. Use `--adapter-preset none` to disable adapter trimming.
- 3' quality trimming (`--trim-min-q`): trims low-quality suffixes using the configured Phred cutoff.

For paired-end input, the CLI currently uses the `DropPair` orphan policy during paired processing: if one mate is rejected, the surviving mate is suppressed rather than emitted. With `--merge-pairs`, successfully merged pairs become single records and continue through the rest of the pipeline as single records. Pairs that cannot be merged continue through paired processing as the original mate pair. Because a merge-enabled run may contain a mix of merged single records and unmerged pairs, `--merge-pairs` requires paired input and single-stream output rather than split `--out1`/`--out2` files. The execution engine also has an internal `EmitOrphan` policy, but orphan emission is not exposed as a supported CLI/output mode yet.

Local input can be single-end (`--in reads.fastq.gz`), interleaved paired-end (`--in reads.fastq.gz --paired`), or split paired-end (`--in1 reads_1.fastq.gz --in2 reads_2.fastq.gz`). Paired output can be split with `--out1`/`--out2` or written as one stream with `--interleaved-out`.

Invalid FASTQ handling is controlled with `--invalid-fastq-policy`. Recoverable invalid records or pairs, such as mate ID mismatches, can fail immediately, warn and drop, or silently drop depending on the selected policy. Parser-level corruption that cannot be safely resynchronized is still reported through the same policy and optional `--invalid-fastq-report` JSONL file, but remains fatal.

## Contributing to nuclease

Nuclease exists in part because I like developing tools in Rust with future hacking in mind. That ethos is pervasive in nuclease, with clear pathways for adding new operations to be run on FASTQs via the `ReadFilter` and `ReadTransform` traits.

> [!NOTE]
> The following demos are quite nice in my opinion...but they will only seem that way if you know Rust!

### Extending nuclease with new filters

In `nuclease`, most operations are _filters_, where some computation is run to decide whether a read should be passed on the in the stream or not. These operations are represented with the `ReadFilter` trait.

To see how this works, let's implement a new filter for GC-contents. Filters each start as a Rust struct that takes responsibility for filtering:

```rust
pub(crate) struct GcRange {
    min_gc: f64,
    max_gc: f64,
}

impl GcRange {
    pub(crate) const fn new(min_gc: f64, max_gc: f64) -> Self {
        Self { min_gc, max_gc }
    }
}
```

`nuclease` needs to know why this struct may reject a read. We can use a marker struct to carry this information (you'll see why in a moment):

```rust
use crate::plan::RejectionReason;

#[derive(Debug)]
pub(crate) struct GcOutOfRange;

impl RejectionReason for GcOutOfRange {
    fn code(&self) -> &'static str {
        "gc_out_of_range"
    }
}
```

Let's actually implement the filter now:

```rust
use crate::{plan::ReadFilter, record::RecordView};

impl ReadFilter for GcRange {
    type Reason = GcOutOfRange;

    fn evaluate(&self, record: &RecordView<'_>) -> Result<(), Self::Reason> {
        let seq = record.sequence();

        if seq.is_empty() {
            return Ok(());
        }

        let gc = seq
            .iter()
            .filter(|&&base| matches!(base, b'G' | b'g' | b'C' | b'c'))
            .count() as f64
            / seq.len() as f64;

        if gc < self.min_gc || gc > self.max_gc {
            Err(GcOutOfRange)
        } else {
            Ok(())
        }
    }
}
```

That's technically all we need to use the filter in `nuclease`. But let's go one step further and bring it into our composable fluent interface:

```rust
use crate::plan::{BuildPlan, ExecutionStep, FilterStep, IntoExecutionStep};

impl IntoExecutionStep for GcRange {
    fn into_execution_step(self) -> Box<dyn ExecutionStep> {
        Box::new(FilterStep(self))
    }
}

pub(crate) trait GcRangeFilter: BuildPlan {
    fn gc_range(self, min_gc: f64, max_gc: f64) -> Self {
        self.step(GcRange::new(min_gc, max_gc))
    }
}

impl<T> GcRangeFilter for T where T: BuildPlan {}
```

This teaches the execution planner that `GcRange` is a filter step. And then making the GcRangeFilter trait and blanket implementing it for all `BuildPlan` types means we can now include our new operation in the `nuclease` plan builder:

```rust
use crate::filter::{
    GcRangeFilter, MaxNsFilter, MinEntropyFilter, MinLengthFilter, MinMeanQualityFilter,
};

let plan = Plan::new()
    // our new filter is now available!
    .gc_range(cli.min_gc, cli.max_gc)

    // from there, we can chain on any number of other desired steps to the plan
    .max_ns(cli.max_ns)
    .trim_adapters(AdapterCatalog::illumina_truseq())
    .quality_trim(cli.trim_min_q)
    .min_length(cli.min_length)
    .min_mean_q(cli.min_mean_q)
    .min_entropy(cli.min_entropy)
    .orphan_policy(OrphanPolicy::DropPair)
    .compile();
```

Transforms that produce mutated versions of the input reads work exactly the same way.

### Extending nuclease with new transforms

Transforms produce a possibly changed `RecordView` rather than accepting or rejecting the read. They are represented with the `ReadTransform` trait. A simple example is trimming a fixed number of bases from the 5' end of each read.

The transform itself is another small Rust struct:

```rust
pub(crate) struct TrimPrefix {
    bases: usize,
}

impl TrimPrefix {
    pub(crate) const fn new(bases: usize) -> Self {
        Self { bases }
    }
}
```

The `ReadTransform` implementation receives a borrowed `RecordView` and a reusable `TransformArena`. If the transform changes the read, it places the changed sequence and quality slices in the arena and returns a new view:

```rust
use crate::{
    plan::{ReadTransform, TransformArena, TransformResult},
    record::RecordView,
};

impl ReadTransform for TrimPrefix {
    fn code(&self) -> &'static str {
        "trim_prefix"
    }

    fn apply<'a>(&self, record: RecordView<'a>, arena: &'a TransformArena) -> TransformResult<'a> {
        assert_eq!(
            record.sequence().len(),
            record.quality().len(),
            "prefix trimming requires equal sequence and quality lengths"
        );

        let trim_start = self.bases.min(record.sequence().len());

        if trim_start == 0 {
            return TransformResult {
                record,
                applied: false,
            };
        }

        TransformResult {
            record: record
                .with_sequence_and_quality(
                    arena.alloc_slice_copy(&record.sequence()[trim_start..]),
                    arena.alloc_slice_copy(&record.quality()[trim_start..]),
                )
                .expect("prefix trimming should preserve equal sequence and quality lengths"),
            applied: true,
        }
    }
}
```

As with filters, bridge the transform into the executor and expose a fluent method:

```rust
use crate::plan::{BuildPlan, ExecutionStep, IntoExecutionStep, TransformStep};

impl IntoExecutionStep for TrimPrefix {
    fn into_execution_step(self) -> Box<dyn ExecutionStep> {
        Box::new(TransformStep(self))
    }
}

pub(crate) trait TrimPrefixTransform: BuildPlan {
    fn trim_prefix(self, bases: usize) -> Self {
        self.step(TrimPrefix::new(bases))
    }
}

impl<T> TrimPrefixTransform for T where T: BuildPlan {}
```

If the transform should be user-configurable, add a CLI field:

```rust
#[arg(long, default_value_t = 0, help_heading = "Preprocessing")]
pub trim_prefix: usize,
```

And compose it into the same plan chain:

```rust
use crate::quality::{QualityTrimTransform, TrimPrefixTransform};

let plan = Plan::new()
    .max_ns(cli.max_ns)
    .trim_adapters(AdapterCatalog::illumina_truseq())
    .trim_prefix(cli.trim_prefix) // here's the new transform
    .quality_trim(cli.trim_min_q)
    .min_length(cli.min_length)
    .min_mean_q(cli.min_mean_q)
    .min_entropy(cli.min_entropy)
    .orphan_policy(OrphanPolicy::DropPair)
    .compile();
```

The transform reports whether it was actually applied. When `applied` is true, the executor records the transform code (`trim_prefix`) in run statistics without the pipeline needing transform-specific accounting.

In both cases, the following steps are all that's needed to do new things in `nuclease`:

1. Make a new struct to handle the transform or filter.
2. Make a struct that records the reason for rejection, if a filter. This isn't necessary for transforms.
3. Implement `ReadFilter` or `ReadTransform`.
4. Implement `IntoExecutionStep` for the new filter/transform struct.
5. Package the struct in a trait and blanket implement it for all plans.

### Roadmap

If you're guessing I have plans to make use of this extensibility, you'd be right! These plans currently include:

- [ ] Internal integration with [Deacon](https://github.com/bede/deacon) to allow decontamination or enrichment with Deacon indices.
- [x] Support for paired-read overlap merging.
- [ ] Support for memory-frugal sorting of reads by sequence homology to improve compression.
- [ ] More adapter trimming support.
- [ ] Read demultiplexing.
