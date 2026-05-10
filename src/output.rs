//! Output selection, construction, and record-writing abstractions.

use std::{
    fs::File,
    io::{self, BufWriter, Write},
    marker::PhantomData,
    path::{Path, PathBuf},
};

use clap::ValueEnum;
use color_eyre::eyre::{Result, WrapErr, bail, eyre};
use flate2::{Compression, write::GzEncoder};

use crate::{
    plan::{EmittedUnit, ExecutionOutcome},
    record::{ReadStats, RecordView, SequenceRecordRef},
};

/// Raw output arguments selected by the CLI before they are resolved into valid output handles.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct OutputArgs {
    interleaved: bool,
    format: OutputFormat,
    encoding: Option<OutputEncoding>,
    out: Option<PathBuf>,
    out1: Option<PathBuf>,
    out2: Option<PathBuf>,
    merge_pairs: bool,
}

impl OutputArgs {
    /// Construct raw output arguments from CLI-provided selections.
    pub fn new(
        interleaved: bool,
        format: OutputFormat,
        encoding: Option<OutputEncoding>,
        out: Option<PathBuf>,
        out1: Option<PathBuf>,
        out2: Option<PathBuf>,
    ) -> Self {
        Self {
            interleaved,
            format,
            encoding,
            out,
            out1,
            out2,
            merge_pairs: false,
        }
    }

    /// Preserve the simple paired-output resolution call while carrying plan cardinality context.
    pub(crate) fn with_merge_pairs(mut self, merge_pairs: bool) -> Self {
        self.merge_pairs = merge_pairs;
        self
    }

    fn resolved_single_encoding(&self) -> OutputEncoding {
        if let Some(encoding) = self.encoding {
            return encoding;
        }

        self.out
            .as_ref()
            .map_or(OutputEncoding::Plain, |path| infer_encoding_from_path(path))
    }

    fn resolved_paired_encoding(&self) -> Result<OutputEncoding> {
        if let Some(encoding) = self.encoding {
            return Ok(encoding);
        }

        match (&self.out, &self.out1, &self.out2) {
            (Some(path), None, None) => Ok(infer_encoding_from_path(path)),
            (None, Some(r1), Some(r2)) => {
                let e1 = infer_encoding_from_path(r1);
                let e2 = infer_encoding_from_path(r2);
                if e1 != e2 {
                    bail!(
                        "paired output paths imply different encodings\n\
                         out1: {} ({e1:?})\n\
                         out2: {} ({e2:?})\n\
                         help: use matching suffixes or specify --output-encoding explicitly",
                        r1.display(),
                        r2.display(),
                    )
                }
                Ok(e1)
            }
            _ => Ok(OutputEncoding::Plain),
        }
    }

    fn validate_paired_merge_output(&self) -> Result<()> {
        if self.merge_pairs && (self.out1.is_some() || self.out2.is_some()) {
            bail!(
                "--merge-pairs cannot be used with split paired output yet\n\
                 help: merged and unmerged reads require one output stream; use --interleaved or --out"
            );
        }

        Ok(())
    }

    /// Resolve runtime-selected single-end output arguments into an opened output handle.
    ///
    /// Missing format and encoding selections are filled in by builder defaults when this method
    /// drives the typestate builder.
    ///
    /// # Errors
    ///
    /// Returns an error when the raw output arguments describe an invalid single-end output
    /// combination or when the selected sink cannot be opened.
    pub fn resolve_single(self) -> Result<SingleOutputHandle> {
        let encoding = self.resolved_single_encoding();
        match (&self.out, &self.out1, &self.out2, self.interleaved) {
            (_, _, _, true) => bail!(
                "single-end ingress cannot be emitted as interleaved output\n\
                 help: remove --interleaved, or provide paired input with --in1 and --in2"
            ),
            (None, None, None, false) => match (self.format, encoding) {
                (OutputFormat::Fastq, OutputEncoding::Plain) => {
                    OutputBuilder::new().single().stdout().build()
                }
                (OutputFormat::Fastq, OutputEncoding::Gzip) => {
                    OutputBuilder::new().single().stdout().gzip().build()
                }
                (OutputFormat::Fasta, OutputEncoding::Plain) => {
                    OutputBuilder::new().single().stdout().fasta().build()
                }
                (OutputFormat::Fasta, OutputEncoding::Gzip) => OutputBuilder::new()
                    .single()
                    .stdout()
                    .fasta()
                    .gzip()
                    .build(),
            },
            (Some(path), None, None, false) => match (self.format, encoding) {
                (OutputFormat::Fastq, OutputEncoding::Plain) => {
                    OutputBuilder::new().single().file(path.clone()).build()
                }
                (OutputFormat::Fastq, OutputEncoding::Gzip) => OutputBuilder::new()
                    .single()
                    .file(path.clone())
                    .gzip()
                    .build(),
                (OutputFormat::Fasta, OutputEncoding::Plain) => OutputBuilder::new()
                    .single()
                    .file(path.clone())
                    .fasta()
                    .build(),
                (OutputFormat::Fasta, OutputEncoding::Gzip) => OutputBuilder::new()
                    .single()
                    .file(path.clone())
                    .fasta()
                    .gzip()
                    .build(),
            },
            _ => bail!(
                "single-end output accepts either stdout or --out only\n\
                 help: remove --out1/--out2 for single-end input"
            ),
        }
    }

    /// Resolve runtime-selected paired-end output arguments into an opened output handle.
    ///
    /// Missing format and encoding selections are filled in by builder defaults when this method
    /// drives the typestate builder.
    ///
    /// # Errors
    ///
    /// Returns an error when the raw output arguments describe an invalid paired-end output
    /// combination or when the selected sinks cannot be opened.
    pub fn resolve_paired(self) -> Result<PairedOutputHandle> {
        self.validate_paired_merge_output()?;
        let encoding = self.resolved_paired_encoding()?;
        match (self.interleaved, &self.out, &self.out1, &self.out2) {
            (true, None, Some(_), Some(_)) => {
                bail!(
                    "interleaved paired output uses stdout or --out, not --out1/--out2\n\
                     help: either remove --interleaved for split paired files, or replace --out1/--out2 with --out"
                )
            }
            (false, None, None, None) => bail!(
                "split paired output requires --out1 and --out2\n\
                 help: provide both split output paths, or add --interleaved to write paired reads to stdout"
            ),
            (false, Some(_), _, _) => bail!(
                "split paired output cannot use --out\n\
                 help: use --out1/--out2 for split paired output, or add --interleaved when using --out"
            ),
            (true, None, None, None) => match (self.format, encoding) {
                (OutputFormat::Fastq, OutputEncoding::Plain) => {
                    OutputBuilder::new().paired().interleaved_stdout().build()
                }
                (OutputFormat::Fastq, OutputEncoding::Gzip) => OutputBuilder::new()
                    .paired()
                    .interleaved_stdout()
                    .gzip()
                    .build(),
                (OutputFormat::Fasta, OutputEncoding::Plain) => OutputBuilder::new()
                    .paired()
                    .interleaved_stdout()
                    .fasta()
                    .build(),
                (OutputFormat::Fasta, OutputEncoding::Gzip) => OutputBuilder::new()
                    .paired()
                    .interleaved_stdout()
                    .fasta()
                    .gzip()
                    .build(),
            },
            (true, Some(path), None, None) => match (self.format, encoding) {
                (OutputFormat::Fastq, OutputEncoding::Plain) => OutputBuilder::new()
                    .paired()
                    .interleaved_file(path.clone())
                    .build(),
                (OutputFormat::Fastq, OutputEncoding::Gzip) => OutputBuilder::new()
                    .paired()
                    .interleaved_file(path.clone())
                    .gzip()
                    .build(),
                (OutputFormat::Fasta, OutputEncoding::Plain) => OutputBuilder::new()
                    .paired()
                    .interleaved_file(path.clone())
                    .fasta()
                    .build(),
                (OutputFormat::Fasta, OutputEncoding::Gzip) => OutputBuilder::new()
                    .paired()
                    .interleaved_file(path.clone())
                    .fasta()
                    .gzip()
                    .build(),
            },
            (false, None, Some(r1), Some(r2)) => match (self.format, encoding) {
                (OutputFormat::Fastq, OutputEncoding::Plain) => OutputBuilder::new()
                    .paired()
                    .split_files(r1.clone(), r2.clone())
                    .build(),
                (OutputFormat::Fastq, OutputEncoding::Gzip) => OutputBuilder::new()
                    .paired()
                    .split_files(r1.clone(), r2.clone())
                    .gzip()
                    .build(),
                (OutputFormat::Fasta, OutputEncoding::Plain) => OutputBuilder::new()
                    .paired()
                    .split_files(r1.clone(), r2.clone())
                    .fasta()
                    .build(),
                (OutputFormat::Fasta, OutputEncoding::Gzip) => OutputBuilder::new()
                    .paired()
                    .split_files(r1.clone(), r2.clone())
                    .fasta()
                    .gzip()
                    .build(),
            },
            _ => bail!(
                "invalid paired output destination combination\n\
                 help: paired output must be either split (--out1 and --out2) or interleaved (--interleaved with stdout or --out)"
            ),
        }
    }
}

fn infer_encoding_from_path(path: &Path) -> OutputEncoding {
    if path.extension().is_some_and(|ext| ext == "gz") {
        OutputEncoding::Gzip
    } else {
        OutputEncoding::Plain
    }
}

/// Output record format.
#[derive(Copy, Clone, Debug, Eq, PartialEq, ValueEnum)]
pub enum OutputFormat {
    /// Emit FASTQ records with qualities.
    Fastq,
    /// Emit FASTA records without qualities.
    Fasta,
}

/// Output encoding wrapped around the selected output format.
#[derive(Copy, Clone, Debug, Eq, PartialEq, ValueEnum)]
pub enum OutputEncoding {
    /// Plain uncompressed bytes.
    Plain,
    /// Gzip-compressed bytes.
    Gzip,
}

/// Typestate marker for an unset builder axis.
pub struct Unset;
/// Typestate marker for single-end output.
pub struct Single;
/// Typestate marker for paired-end output.
pub struct Paired;
/// Typestate marker for stdout destination.
pub struct Stdout;
/// Typestate marker for single-file destination.
pub struct FileOutput;
/// Typestate marker for split paired-file destination.
pub struct FilePair;
/// Typestate marker for default output format.
pub struct DefaultFormat;
/// Typestate marker for default output encoding.
pub struct DefaultEncoding;
/// Typestate marker for FASTQ format.
pub struct Fastq;
/// Typestate marker for FASTA format.
pub struct Fasta;
/// Typestate marker for plain output encoding.
pub struct Plain;
/// Typestate marker for gzip output encoding.
pub struct Gzip;

/// Typestate builder for output construction.
pub struct OutputBuilder<L, D, F, E> {
    out: Option<PathBuf>,
    out1: Option<PathBuf>,
    out2: Option<PathBuf>,
    _marker: PhantomData<(L, D, F, E)>,
}

impl OutputBuilder<Unset, Unset, Unset, Unset> {
    /// Construct a new typestate builder with no selections applied yet.
    pub const fn new() -> Self {
        Self {
            out: None,
            out1: None,
            out2: None,
            _marker: PhantomData,
        }
    }

    /// Begin constructing a single-end output.
    pub fn single(self) -> OutputPlan<Single, Unset, DefaultFormat, DefaultEncoding> {
        let OutputBuilder { .. } = self;
        OutputPlan::new(OutputBuilder {
            out: None,
            out1: None,
            out2: None,
            _marker: PhantomData,
        })
    }

    /// Begin constructing a paired-end output.
    pub fn paired(self) -> OutputPlan<Paired, Unset, DefaultFormat, DefaultEncoding> {
        let OutputBuilder { .. } = self;
        OutputPlan::new(OutputBuilder {
            out: None,
            out1: None,
            out2: None,
            _marker: PhantomData,
        })
    }
}

/// Typestate wrapper for a partially or fully specified output configuration.
pub struct OutputPlan<L, D, F, E> {
    state: OutputBuilder<L, D, F, E>,
}

impl<L, D, F, E> OutputPlan<L, D, F, E> {
    fn new(state: OutputBuilder<L, D, F, E>) -> Self {
        Self { state }
    }
}

impl<F, E> OutputPlan<Single, Unset, F, E> {
    /// Direct single-end output to stdout.
    pub fn stdout(self) -> OutputPlan<Single, Stdout, F, E> {
        let OutputPlan { state: _ } = self;
        OutputPlan::new(OutputBuilder {
            out: None,
            out1: None,
            out2: None,
            _marker: PhantomData,
        })
    }

    /// Direct single-end output to one file.
    pub fn file(self, path: PathBuf) -> OutputPlan<Single, FileOutput, F, E> {
        let OutputPlan { state: _ } = self;
        OutputPlan::new(OutputBuilder {
            out: Some(path),
            out1: None,
            out2: None,
            _marker: PhantomData,
        })
    }
}

impl<F, E> OutputPlan<Paired, Unset, F, E> {
    /// Direct paired-end output to one interleaved stdout stream.
    pub fn interleaved_stdout(self) -> OutputPlan<Paired, Stdout, F, E> {
        let OutputPlan { state: _ } = self;
        OutputPlan::new(OutputBuilder {
            out: None,
            out1: None,
            out2: None,
            _marker: PhantomData,
        })
    }

    /// Direct paired-end output to one interleaved file.
    pub fn interleaved_file(self, path: PathBuf) -> OutputPlan<Paired, FileOutput, F, E> {
        let OutputPlan { state: _ } = self;
        OutputPlan::new(OutputBuilder {
            out: Some(path),
            out1: None,
            out2: None,
            _marker: PhantomData,
        })
    }

    /// Direct paired-end output to two split files.
    pub fn split_files(self, r1: PathBuf, r2: PathBuf) -> OutputPlan<Paired, FilePair, F, E> {
        let OutputPlan { state: _ } = self;
        OutputPlan::new(OutputBuilder {
            out: None,
            out1: Some(r1),
            out2: Some(r2),
            _marker: PhantomData,
        })
    }
}

impl<L, D, E> OutputPlan<L, D, DefaultFormat, E> {
    /// Select FASTQ output.
    pub fn fastq(self) -> OutputPlan<L, D, Fastq, E> {
        OutputPlan::new(OutputBuilder {
            out: self.state.out,
            out1: self.state.out1,
            out2: self.state.out2,
            _marker: PhantomData,
        })
    }

    /// Select FASTA output.
    pub fn fasta(self) -> OutputPlan<L, D, Fasta, E> {
        OutputPlan::new(OutputBuilder {
            out: self.state.out,
            out1: self.state.out1,
            out2: self.state.out2,
            _marker: PhantomData,
        })
    }
}

impl<L, D, F> OutputPlan<L, D, F, DefaultEncoding> {
    /// Select plain uncompressed output.
    pub fn plain(self) -> OutputPlan<L, D, F, Plain> {
        OutputPlan::new(OutputBuilder {
            out: self.state.out,
            out1: self.state.out1,
            out2: self.state.out2,
            _marker: PhantomData,
        })
    }

    /// Select gzip-compressed output.
    pub fn gzip(self) -> OutputPlan<L, D, F, Gzip> {
        OutputPlan::new(OutputBuilder {
            out: self.state.out,
            out1: self.state.out1,
            out2: self.state.out2,
            _marker: PhantomData,
        })
    }
}

impl<L, D> OutputPlan<L, D, DefaultFormat, DefaultEncoding>
where
    OutputPlan<L, D, Fastq, Plain>: BuildHandle,
{
    /// Build using the default output format and encoding: FASTQ + plain.
    pub fn build(self) -> Result<<OutputPlan<L, D, Fastq, Plain> as BuildHandle>::Handle> {
        self.fastq().plain().build()
    }
}

impl<L, D, F> OutputPlan<L, D, F, DefaultEncoding>
where
    OutputPlan<L, D, F, Plain>: BuildHandle,
{
    /// Build using the default output encoding for the selected format: plain.
    pub fn build(self) -> Result<<OutputPlan<L, D, F, Plain> as BuildHandle>::Handle> {
        self.plain().build()
    }
}

impl<L, D, E> OutputPlan<L, D, DefaultFormat, E>
where
    OutputPlan<L, D, Fastq, E>: BuildHandle,
{
    /// Build using the default output format for the selected encoding: FASTQ.
    pub fn build(self) -> Result<<OutputPlan<L, D, Fastq, E> as BuildHandle>::Handle> {
        self.fastq().build()
    }
}

/// Trait implemented by fully specified output plans that can open concrete handles.
pub(crate) trait BuildHandle {
    /// The output handle produced by the fully resolved plan.
    type Handle;

    /// Open the selected output sink or sinks.
    ///
    /// # Errors
    ///
    /// Returns an error when the configured sink or sinks cannot be opened.
    fn build(self) -> Result<Self::Handle>;
}

fn build_single_handle(
    format: OutputFormat,
    encoding: OutputEncoding,
    out: Option<PathBuf>,
) -> Result<SingleOutputHandle> {
    let writer = match out {
        Some(path) => finalizable_file_writer(&path, encoding)?,
        None => finalizable_stdout_writer(encoding),
    };

    Ok(SingleOutputHandle {
        inner: SingleOutput::new(StreamSink::new(writer, format)),
    })
}

fn build_interleaved_paired_handle(
    format: OutputFormat,
    encoding: OutputEncoding,
    out: Option<PathBuf>,
) -> Result<PairedOutputHandle> {
    let writer = match out {
        Some(path) => finalizable_file_writer(&path, encoding)?,
        None => finalizable_stdout_writer(encoding),
    };

    Ok(PairedOutputHandle::Interleaved(Box::new(
        InterleavedOutputHandle {
            inner: InterleavedOutput::new(StreamSink::new(writer, format)),
        },
    )))
}

fn build_split_paired_handle(
    format: OutputFormat,
    encoding: OutputEncoding,
    r1: &PathBuf,
    r2: &PathBuf,
) -> Result<PairedOutputHandle> {
    let w1 = finalizable_file_writer(r1, encoding)?;
    let w2 = finalizable_file_writer(r2, encoding)?;

    Ok(PairedOutputHandle::Split(Box::new(SplitOutputHandle {
        inner: SplitOutput::new(StreamSink::new(w1, format), StreamSink::new(w2, format)),
    })))
}

impl BuildHandle for OutputPlan<Single, Stdout, Fastq, Plain> {
    type Handle = SingleOutputHandle;

    fn build(self) -> Result<Self::Handle> {
        build_single_handle(OutputFormat::Fastq, OutputEncoding::Plain, None)
    }
}

impl BuildHandle for OutputPlan<Single, Stdout, Fastq, Gzip> {
    type Handle = SingleOutputHandle;

    fn build(self) -> Result<Self::Handle> {
        build_single_handle(OutputFormat::Fastq, OutputEncoding::Gzip, None)
    }
}

impl BuildHandle for OutputPlan<Single, Stdout, Fasta, Plain> {
    type Handle = SingleOutputHandle;

    fn build(self) -> Result<Self::Handle> {
        build_single_handle(OutputFormat::Fasta, OutputEncoding::Plain, None)
    }
}

impl BuildHandle for OutputPlan<Single, Stdout, Fasta, Gzip> {
    type Handle = SingleOutputHandle;

    fn build(self) -> Result<Self::Handle> {
        build_single_handle(OutputFormat::Fasta, OutputEncoding::Gzip, None)
    }
}

impl BuildHandle for OutputPlan<Single, FileOutput, Fastq, Plain> {
    type Handle = SingleOutputHandle;

    fn build(self) -> Result<Self::Handle> {
        build_single_handle(OutputFormat::Fastq, OutputEncoding::Plain, self.state.out)
    }
}

impl BuildHandle for OutputPlan<Single, FileOutput, Fastq, Gzip> {
    type Handle = SingleOutputHandle;

    fn build(self) -> Result<Self::Handle> {
        build_single_handle(OutputFormat::Fastq, OutputEncoding::Gzip, self.state.out)
    }
}

impl BuildHandle for OutputPlan<Single, FileOutput, Fasta, Plain> {
    type Handle = SingleOutputHandle;

    fn build(self) -> Result<Self::Handle> {
        build_single_handle(OutputFormat::Fasta, OutputEncoding::Plain, self.state.out)
    }
}

impl BuildHandle for OutputPlan<Single, FileOutput, Fasta, Gzip> {
    type Handle = SingleOutputHandle;

    fn build(self) -> Result<Self::Handle> {
        build_single_handle(OutputFormat::Fasta, OutputEncoding::Gzip, self.state.out)
    }
}

impl BuildHandle for OutputPlan<Paired, Stdout, Fastq, Plain> {
    type Handle = PairedOutputHandle;

    fn build(self) -> Result<Self::Handle> {
        build_interleaved_paired_handle(OutputFormat::Fastq, OutputEncoding::Plain, None)
    }
}

impl BuildHandle for OutputPlan<Paired, Stdout, Fastq, Gzip> {
    type Handle = PairedOutputHandle;

    fn build(self) -> Result<Self::Handle> {
        build_interleaved_paired_handle(OutputFormat::Fastq, OutputEncoding::Gzip, None)
    }
}

impl BuildHandle for OutputPlan<Paired, Stdout, Fasta, Plain> {
    type Handle = PairedOutputHandle;

    fn build(self) -> Result<Self::Handle> {
        build_interleaved_paired_handle(OutputFormat::Fasta, OutputEncoding::Plain, None)
    }
}

impl BuildHandle for OutputPlan<Paired, Stdout, Fasta, Gzip> {
    type Handle = PairedOutputHandle;

    fn build(self) -> Result<Self::Handle> {
        build_interleaved_paired_handle(OutputFormat::Fasta, OutputEncoding::Gzip, None)
    }
}

impl BuildHandle for OutputPlan<Paired, FileOutput, Fastq, Plain> {
    type Handle = PairedOutputHandle;

    fn build(self) -> Result<Self::Handle> {
        build_interleaved_paired_handle(OutputFormat::Fastq, OutputEncoding::Plain, self.state.out)
    }
}

impl BuildHandle for OutputPlan<Paired, FileOutput, Fastq, Gzip> {
    type Handle = PairedOutputHandle;

    fn build(self) -> Result<Self::Handle> {
        build_interleaved_paired_handle(OutputFormat::Fastq, OutputEncoding::Gzip, self.state.out)
    }
}

impl BuildHandle for OutputPlan<Paired, FileOutput, Fasta, Plain> {
    type Handle = PairedOutputHandle;

    fn build(self) -> Result<Self::Handle> {
        build_interleaved_paired_handle(OutputFormat::Fasta, OutputEncoding::Plain, self.state.out)
    }
}

impl BuildHandle for OutputPlan<Paired, FileOutput, Fasta, Gzip> {
    type Handle = PairedOutputHandle;

    fn build(self) -> Result<Self::Handle> {
        build_interleaved_paired_handle(OutputFormat::Fasta, OutputEncoding::Gzip, self.state.out)
    }
}

impl BuildHandle for OutputPlan<Paired, FilePair, Fastq, Plain> {
    type Handle = PairedOutputHandle;

    fn build(self) -> Result<Self::Handle> {
        build_split_paired_handle(
            OutputFormat::Fastq,
            OutputEncoding::Plain,
            self.state
                .out1
                .as_ref()
                .expect("split file mate 1 path must be present"),
            self.state
                .out2
                .as_ref()
                .expect("split file mate 2 path must be present"),
        )
    }
}

impl BuildHandle for OutputPlan<Paired, FilePair, Fastq, Gzip> {
    type Handle = PairedOutputHandle;

    fn build(self) -> Result<Self::Handle> {
        build_split_paired_handle(
            OutputFormat::Fastq,
            OutputEncoding::Gzip,
            self.state
                .out1
                .as_ref()
                .expect("split file mate 1 path must be present"),
            self.state
                .out2
                .as_ref()
                .expect("split file mate 2 path must be present"),
        )
    }
}

impl BuildHandle for OutputPlan<Paired, FilePair, Fasta, Plain> {
    type Handle = PairedOutputHandle;

    fn build(self) -> Result<Self::Handle> {
        build_split_paired_handle(
            OutputFormat::Fasta,
            OutputEncoding::Plain,
            self.state
                .out1
                .as_ref()
                .expect("split file mate 1 path must be present"),
            self.state
                .out2
                .as_ref()
                .expect("split file mate 2 path must be present"),
        )
    }
}

impl BuildHandle for OutputPlan<Paired, FilePair, Fasta, Gzip> {
    type Handle = PairedOutputHandle;

    fn build(self) -> Result<Self::Handle> {
        build_split_paired_handle(
            OutputFormat::Fasta,
            OutputEncoding::Gzip,
            self.state
                .out1
                .as_ref()
                .expect("split file mate 1 path must be present"),
            self.state
                .out2
                .as_ref()
                .expect("split file mate 2 path must be present"),
        )
    }
}

/// Generic sink that can write records implementing [`SequenceRecordRef`].
pub trait RecordSink<R> {
    /// Write one logical record to the sink.
    fn write_record(&mut self, record: R) -> Result<()>;
}

/// Record sink over one writable byte stream with a selected output format.
pub struct StreamSink<W> {
    writer: W,
    format: OutputFormat,
}

impl<W> StreamSink<W>
where
    W: Write,
{
    /// Construct a stream sink from a writable byte stream and output format.
    pub fn new(writer: W, format: OutputFormat) -> Self {
        Self { writer, format }
    }

    pub(crate) fn into_inner(self) -> W {
        self.writer
    }
}

impl<W, R> RecordSink<R> for StreamSink<W>
where
    W: Write,
    R: SequenceRecordRef,
{
    fn write_record(&mut self, record: R) -> Result<()> {
        match self.format {
            OutputFormat::Fastq => {
                let quality = record
                    .quality()
                    .ok_or_else(|| eyre!("FASTQ output requires quality scores"))?;

                self.writer.write_all(b"@")?;
                self.writer.write_all(record.header())?;
                self.writer.write_all(b"\n")?;
                self.writer.write_all(record.sequence())?;
                self.writer.write_all(b"\n+\n")?;
                self.writer.write_all(quality)?;
                self.writer.write_all(b"\n")?;
            }
            OutputFormat::Fasta => {
                self.writer.write_all(b">")?;
                self.writer.write_all(record.header())?;
                self.writer.write_all(b"\n")?;
                self.writer.write_all(record.sequence())?;
                self.writer.write_all(b"\n")?;
            }
        }

        Ok(())
    }
}

/// Single-stream output arrangement.
pub struct SingleOutput<S> {
    sink: S,
}

impl<S> SingleOutput<S> {
    /// Construct a single-output arrangement from a sink.
    pub fn new(sink: S) -> Self {
        Self { sink }
    }

    pub(crate) fn into_inner(self) -> S {
        self.sink
    }
}

/// Trait for things that can consume one record at a time.
pub trait SingleRecordOutput {
    /// Write one record.
    fn write_record(&mut self, record: RecordView<'_>) -> Result<()>;
}

/// Trait for outputs that can consume a full plan execution outcome.
pub(crate) trait UnitOutput {
    /// Write the emitted records in one outcome and update output counters.
    fn write_outcome(
        &mut self,
        outcome: &ExecutionOutcome<'_>,
        stats: &mut ReadStats,
    ) -> Result<()>;
}

impl<S> SingleRecordOutput for SingleOutput<S>
where
    S: for<'a> RecordSink<RecordView<'a>>,
{
    fn write_record(&mut self, record: RecordView<'_>) -> Result<()> {
        self.sink.write_record(record)
    }
}

/// Interleaved paired-output arrangement that writes both mates into one sink.
pub struct InterleavedOutput<S> {
    sink: S,
}

impl<S> InterleavedOutput<S> {
    /// Construct an interleaved paired-output arrangement from one sink.
    pub fn new(sink: S) -> Self {
        Self { sink }
    }

    pub(crate) fn into_inner(self) -> S {
        self.sink
    }
}

/// Split paired-output arrangement that writes mates into separate sinks.
pub struct SplitOutput<S1, S2> {
    r1: S1,
    r2: S2,
}

impl<S1, S2> SplitOutput<S1, S2> {
    /// Construct a split paired-output arrangement from mate-specific sinks.
    pub fn new(r1: S1, r2: S2) -> Self {
        Self { r1, r2 }
    }

    pub(crate) fn into_parts(self) -> (S1, S2) {
        (self.r1, self.r2)
    }
}

/// Trait for things that can consume paired records.
pub trait PairedRecordOutput {
    /// Write one logical pair of records.
    fn write_pair(&mut self, r1: RecordView<'_>, r2: RecordView<'_>) -> Result<()>;
}

impl<S> PairedRecordOutput for InterleavedOutput<S>
where
    S: for<'a> RecordSink<RecordView<'a>>,
{
    fn write_pair(&mut self, r1: RecordView<'_>, r2: RecordView<'_>) -> Result<()> {
        self.sink.write_record(r1)?;
        self.sink.write_record(r2)?;
        Ok(())
    }
}

impl<S1, S2> PairedRecordOutput for SplitOutput<S1, S2>
where
    S1: for<'a> RecordSink<RecordView<'a>>,
    S2: for<'a> RecordSink<RecordView<'a>>,
{
    fn write_pair(&mut self, r1: RecordView<'_>, r2: RecordView<'_>) -> Result<()> {
        self.r1.write_record(r1)?;
        self.r2.write_record(r2)?;
        Ok(())
    }
}

/// Internal trait for byte writers that need an explicit end-of-stream finalization step.
trait FinishableWrite: Write {
    /// Finalize the writer and make the output fully readable.
    fn finish(self: Box<Self>) -> Result<()>;
}

impl FinishableWrite for BufWriter<io::Stdout> {
    fn finish(mut self: Box<Self>) -> Result<()> {
        self.flush()?;
        Ok(())
    }
}

impl FinishableWrite for BufWriter<File> {
    fn finish(mut self: Box<Self>) -> Result<()> {
        self.flush()?;
        Ok(())
    }
}

impl FinishableWrite for GzEncoder<BufWriter<io::Stdout>> {
    fn finish(mut self: Box<Self>) -> Result<()> {
        self.try_finish()?;
        Ok(())
    }
}

impl FinishableWrite for GzEncoder<BufWriter<File>> {
    fn finish(mut self: Box<Self>) -> Result<()> {
        self.try_finish()?;
        Ok(())
    }
}

type FinishableBox = Box<dyn FinishableWrite>;
type StreamOutput = StreamSink<FinishableBox>;
type SingleStreamOutput = SingleOutput<StreamOutput>;
type InterleavedPairedOutput = InterleavedOutput<StreamOutput>;
type SplitPairedOutput = SplitOutput<StreamOutput, StreamOutput>;

/// Opened single-output handle used by the pipeline after a single plan is built.
pub(crate) struct SingleOutputHandle {
    inner: SingleStreamOutput,
}

impl SingleOutputHandle {
    /// Finalize the underlying writer and make the output consumable.
    pub fn finish(self) -> Result<()> {
        self.inner.into_inner().into_inner().finish()
    }
}

impl SingleRecordOutput for SingleOutputHandle {
    fn write_record(&mut self, record: RecordView<'_>) -> Result<()> {
        self.inner.write_record(record)
    }
}

impl UnitOutput for SingleOutputHandle {
    fn write_outcome(
        &mut self,
        outcome: &ExecutionOutcome<'_>,
        stats: &mut ReadStats,
    ) -> Result<()> {
        for record in outcome.emitted() {
            self.write_record(record)?;
            stats.record_emitted(record.sequence().len());
        }
        Ok(())
    }
}

/// Opened interleaved paired-output handle hidden behind [`PairedOutputHandle`].
pub(crate) struct InterleavedOutputHandle {
    inner: InterleavedPairedOutput,
}

impl InterleavedOutputHandle {
    fn finish(self) -> Result<()> {
        self.inner.into_inner().into_inner().finish()
    }
}

/// Opened split paired-output handle hidden behind [`PairedOutputHandle`].
pub(crate) struct SplitOutputHandle {
    inner: SplitPairedOutput,
}

impl SplitOutputHandle {
    fn finish(self) -> Result<()> {
        let (r1, r2) = self.inner.into_parts();
        r1.into_inner().finish()?;
        r2.into_inner().finish()?;
        Ok(())
    }
}

/// Opened paired-output handle used by the pipeline after a paired plan is built.
pub(crate) enum PairedOutputHandle {
    /// Interleaved paired-output sink.
    Interleaved(Box<InterleavedOutputHandle>),
    /// Split paired-output sink.
    Split(Box<SplitOutputHandle>),
}

impl PairedOutputHandle {
    /// Finalize the underlying writer or writers and make the output consumable.
    pub fn finish(self) -> Result<()> {
        match self {
            Self::Interleaved(output) => output.finish(),
            Self::Split(output) => output.finish(),
        }
    }
}

impl UnitOutput for PairedOutputHandle {
    fn write_outcome(
        &mut self,
        outcome: &ExecutionOutcome<'_>,
        stats: &mut ReadStats,
    ) -> Result<()> {
        match outcome.emitted_unit() {
            EmittedUnit::None => {
                if outcome.rejection_count() > 0 {
                    stats.record_pair_rejected();
                }
                Ok(())
            }
            EmittedUnit::Single(record) => match self {
                Self::Interleaved(output) => {
                    output.inner.sink.write_record(record)?;
                    stats.record_emitted(record.sequence().len());
                    stats.record_pair_emitted();
                    Ok(())
                }
                Self::Split(_) => bail!(
                    "split paired output cannot represent a merged single read\n\
                     help: use --interleaved or --out when --merge-pairs is enabled"
                ),
            },
            EmittedUnit::Pair(pair) => {
                self.write_pair(pair.left, pair.right)?;
                stats.record_emitted(pair.left.sequence().len());
                stats.record_emitted(pair.right.sequence().len());
                stats.record_pair_emitted();
                Ok(())
            }
        }
    }
}

impl PairedRecordOutput for PairedOutputHandle {
    fn write_pair(&mut self, r1: RecordView<'_>, r2: RecordView<'_>) -> Result<()> {
        match self {
            Self::Interleaved(output) => output.inner.write_pair(r1, r2),
            Self::Split(output) => output.inner.write_pair(r1, r2),
        }
    }
}

/// Construct a stdout-backed writer with the requested output encoding.
fn finalizable_stdout_writer(encoding: OutputEncoding) -> FinishableBox {
    match encoding {
        OutputEncoding::Plain => Box::new(BufWriter::new(io::stdout())),
        OutputEncoding::Gzip => Box::new(GzEncoder::new(
            BufWriter::new(io::stdout()),
            Compression::default(),
        )),
    }
}

/// Construct a file-backed writer with the requested output encoding.
///
/// # Errors
///
/// Returns an error when the target file cannot be created.
fn finalizable_file_writer(path: &PathBuf, encoding: OutputEncoding) -> Result<FinishableBox> {
    let file = File::create(path).wrap_err_with(|| {
        format!(
            "failed to create output file: {}\n\
             encoding: {encoding:?}\n\
             help: check that the parent directory exists and is writable",
            path.display()
        )
    })?;
    let writer = BufWriter::new(file);

    match encoding {
        OutputEncoding::Plain => Ok(Box::new(writer)),
        OutputEncoding::Gzip => Ok(Box::new(GzEncoder::new(writer, Compression::default()))),
    }
}

#[cfg(test)]
mod tests {
    use std::{
        fs,
        io::{Cursor, Read as _},
        path::PathBuf,
    };

    use color_eyre::Result;
    use flate2::read::GzDecoder;
    use tempfile::tempdir;

    use super::{
        InterleavedOutput, OutputArgs, OutputEncoding, OutputFormat, PairedRecordOutput,
        RecordSink, SingleRecordOutput, StreamSink,
    };
    use crate::record::RecordView;

    #[test]
    fn output_args_resolve_single_stdout_by_default() -> Result<()> {
        let temp = tempdir()?;
        let out = temp.path().join("single.fastq");
        let output = OutputArgs::new(
            false,
            OutputFormat::Fastq,
            Some(OutputEncoding::Plain),
            Some(out.clone()),
            None,
            None,
        );
        let mut handle = output.resolve_single()?;
        handle.write_record(RecordView::new(b"single", b"ACGT", b"IIII"))?;
        handle.finish()?;

        assert_eq!(fs::read_to_string(out)?, "@single\nACGT\n+\nIIII\n");
        Ok(())
    }

    #[test]
    fn output_args_resolve_paired_split_file_targets() -> Result<()> {
        let temp = tempdir()?;
        let r1_path = temp.path().join("r1.fastq");
        let r2_path = temp.path().join("r2.fastq");
        let output = OutputArgs::new(
            false,
            OutputFormat::Fastq,
            Some(OutputEncoding::Plain),
            None,
            Some(r1_path.clone()),
            Some(r2_path.clone()),
        );
        let mut handle = output.resolve_paired()?;
        handle.write_pair(
            RecordView::new(b"r1/1", b"AAAA", b"IIII"),
            RecordView::new(b"r1/2", b"TTTT", b"JJJJ"),
        )?;
        handle.finish()?;

        assert_eq!(fs::read_to_string(r1_path)?, "@r1/1\nAAAA\n+\nIIII\n");
        assert_eq!(fs::read_to_string(r2_path)?, "@r1/2\nTTTT\n+\nJJJJ\n");
        Ok(())
    }

    #[test]
    fn stream_sink_writes_fastq_records() -> Result<()> {
        let mut sink = StreamSink::new(Vec::new(), OutputFormat::Fastq);
        sink.write_record(RecordView::new(b"read1 sample=alpha", b"ACGT", b"IIII"))?;

        assert_eq!(sink.writer, b"@read1 sample=alpha\nACGT\n+\nIIII\n");
        Ok(())
    }

    #[test]
    fn stream_sink_writes_fasta_records() -> Result<()> {
        let mut sink = StreamSink::new(Vec::new(), OutputFormat::Fasta);
        sink.write_record(RecordView::new(b"read1 sample=alpha", b"ACGT", b"IIII"))?;

        assert_eq!(sink.writer, b">read1 sample=alpha\nACGT\n");
        Ok(())
    }

    #[test]
    fn interleaved_output_writes_pairs_into_single_sink() -> Result<()> {
        let sink = StreamSink::new(Cursor::new(Vec::new()), OutputFormat::Fastq);
        let mut output = InterleavedOutput::new(sink);
        output.write_pair(
            RecordView::new(b"r1/1", b"AAAA", b"IIII"),
            RecordView::new(b"r1/2", b"TTTT", b"JJJJ"),
        )?;

        let bytes = output.into_inner().into_inner().into_inner();
        assert_eq!(bytes, b"@r1/1\nAAAA\n+\nIIII\n@r1/2\nTTTT\n+\nJJJJ\n");
        Ok(())
    }

    #[test]
    fn paired_split_plain_output_writes_two_fastq_files() -> Result<()> {
        let temp = tempdir()?;
        let r1_path = temp.path().join("r1.fastq");
        let r2_path = temp.path().join("r2.fastq");

        let output = OutputArgs::new(
            false,
            OutputFormat::Fastq,
            Some(OutputEncoding::Plain),
            None,
            Some(r1_path.clone()),
            Some(r2_path.clone()),
        );
        let mut handle = output.resolve_paired()?;

        handle.write_pair(
            RecordView::new(b"r1/1 sample=alpha", b"AAAA", b"IIII"),
            RecordView::new(b"r1/2 sample=alpha", b"TTTT", b"JJJJ"),
        )?;
        handle.finish()?;

        assert_eq!(
            fs::read_to_string(r1_path)?,
            "@r1/1 sample=alpha\nAAAA\n+\nIIII\n"
        );
        assert_eq!(
            fs::read_to_string(r2_path)?,
            "@r1/2 sample=alpha\nTTTT\n+\nJJJJ\n"
        );
        Ok(())
    }

    #[test]
    fn paired_split_gzip_output_writes_two_gzipped_fastq_files() -> Result<()> {
        let temp = tempdir()?;
        let r1_path = temp.path().join("r1.fastq.gz");
        let r2_path = temp.path().join("r2.fastq.gz");

        let output = OutputArgs::new(
            false,
            OutputFormat::Fastq,
            Some(OutputEncoding::Gzip),
            None,
            Some(r1_path.clone()),
            Some(r2_path.clone()),
        );
        let mut handle = output.resolve_paired()?;

        handle.write_pair(
            RecordView::new(b"r2/1 lane=4", b"CCCC", b"KKKK"),
            RecordView::new(b"r2/2 lane=4", b"GGGG", b"LLLL"),
        )?;
        handle.finish()?;

        assert_eq!(read_gzip_file(&r1_path)?, "@r2/1 lane=4\nCCCC\n+\nKKKK\n");
        assert_eq!(read_gzip_file(&r2_path)?, "@r2/2 lane=4\nGGGG\n+\nLLLL\n");
        Ok(())
    }

    fn read_gzip_file(path: &std::path::Path) -> Result<String> {
        let compressed = fs::read(path)?;
        let mut decoder = GzDecoder::new(compressed.as_slice());
        let mut decompressed = String::new();
        decoder.read_to_string(&mut decompressed)?;
        Ok(decompressed)
    }

    #[test]
    fn single_output_infers_gzip_from_gz_path() -> Result<()> {
        let temp = tempdir()?;
        let out = temp.path().join("single.fastq.gz");
        let output = OutputArgs::new(
            false,
            OutputFormat::Fastq,
            None,
            Some(out.clone()),
            None,
            None,
        );
        let mut handle = output.resolve_single()?;
        handle.write_record(RecordView::new(b"single", b"ACGT", b"IIII"))?;
        handle.finish()?;

        assert_eq!(read_gzip_file(&out)?, "@single\nACGT\n+\nIIII\n");
        Ok(())
    }

    #[test]
    fn paired_output_rejects_inconsistent_inferred_encodings() {
        let output = OutputArgs::new(
            false,
            OutputFormat::Fastq,
            None,
            None,
            Some(PathBuf::from("r1.fastq.gz")),
            Some(PathBuf::from("r2.fastq")),
        );

        let Err(error) = output.resolve_paired() else {
            panic!("mixed paired output suffixes should be rejected");
        };
        assert!(
            error
                .to_string()
                .contains("paired output paths imply different encodings")
        );
    }

    #[test]
    fn paired_output_rejects_split_output_when_merge_pairs_is_enabled() {
        let output = OutputArgs::new(
            false,
            OutputFormat::Fastq,
            Some(OutputEncoding::Plain),
            None,
            Some(PathBuf::from("r1.fastq")),
            Some(PathBuf::from("r2.fastq")),
        )
        .with_merge_pairs(true);

        let Err(error) = output.resolve_paired() else {
            panic!("split paired output cannot represent merged single reads");
        };
        assert!(
            error
                .to_string()
                .contains("--merge-pairs cannot be used with split paired output")
        );
    }
}
