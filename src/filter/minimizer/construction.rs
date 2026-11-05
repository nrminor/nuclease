//! Streaming reference parsing and exact Deacon index construction.

use std::{io::Read, path::Path};

use deacon::{Index, MinimizerVec, Minimizers};
use needletail::{errors::ParseError, parser::Format};

use super::ErrorRecordingReader;
use crate::error::{InternalError, MalformedInputError, Result, UsageError, constructors};

// Reserve bytes, not elements, so both packed widths start at the same size.
const INITIAL_BYTES: usize = 64 * 1024;
const APPEND_CHUNK: usize = 4096;

struct Accumulator<T> {
    values: Vec<T>,
    next_compaction: usize,
}

impl<T: Copy + Ord> Accumulator<T> {
    fn new() -> Self {
        let initial = INITIAL_BYTES / std::mem::size_of::<T>();
        Self {
            values: Vec::with_capacity(initial),
            next_compaction: initial,
        }
    }

    fn append(&mut self, mut incoming: &[T]) {
        while !incoming.is_empty() {
            let count = incoming
                .len()
                .min(APPEND_CHUNK)
                .min(self.next_compaction - self.values.len());
            self.values.extend_from_slice(&incoming[..count]);
            incoming = &incoming[count..];
            if self.values.len() == self.next_compaction {
                self.compact();
                self.next_compaction = self
                    .values
                    .len()
                    .saturating_mul(2)
                    .max(INITIAL_BYTES / std::mem::size_of::<T>());
            }
        }
    }

    fn compact(&mut self) {
        self.values.sort_unstable();
        self.values.dedup();
    }
}

enum PackedAccumulator {
    U64(Accumulator<u64>),
    U128(Accumulator<u128>),
}

impl PackedAccumulator {
    fn new(k: u8) -> Self {
        if k <= 32 {
            Self::U64(Accumulator::new())
        } else {
            Self::U128(Accumulator::new())
        }
    }

    fn append(&mut self, values: &MinimizerVec) {
        match (self, values) {
            (Self::U64(acc), MinimizerVec::U64(values)) => acc.append(values),
            (Self::U128(acc), MinimizerVec::U128(values)) => acc.append(values),
            _ => unreachable!("extractor and accumulator widths are determined by k"),
        }
    }

    fn finish(self) -> MinimizerVec {
        match self {
            Self::U64(mut acc) => {
                acc.compact();
                MinimizerVec::U64(acc.values)
            }
            Self::U128(mut acc) => {
                acc.compact();
                MinimizerVec::U128(acc.values)
            }
        }
    }
}

fn collect_reference(
    reader: &mut (impl Read + Send),
    extractor: &mut Minimizers,
    accumulator: &mut PackedAccumulator,
) -> std::result::Result<(), ParseError> {
    let mut parser = needletail::parse_fastx_reader(reader)?;
    let mut flattened = Vec::new();
    while let Some(record) = parser.next() {
        let record = record?;
        let raw = record.raw_seq();
        let sequence = if record.format() == Format::Fasta
            && raw.iter().any(|&byte| byte == b'\n' || byte == b'\r')
        {
            flattened.clear();
            flattened.reserve(record.num_bases());
            flattened.extend(
                raw.iter()
                    .copied()
                    .filter(|&byte| byte != b'\n' && byte != b'\r'),
            );
            flattened.as_slice()
        } else {
            raw
        };
        accumulator.append(extractor.compute(sequence));
    }
    Ok(())
}

pub(super) fn construct_reference(
    reader: impl Read + Send,
    role: &'static str,
    path: &Path,
    k: u8,
    w: u8,
) -> Result<Index> {
    let mut extractor = Minimizers::new(k, w).map_err(|source| UsageError::IndexParameters {
        k,
        w,
        source: source.into_boxed_dyn_error(),
    })?;
    let mut accumulator = PackedAccumulator::new(k);
    {
        let mut reader = ErrorRecordingReader {
            inner: reader,
            error_kind: None,
        };
        let result = collect_reference(&mut reader, &mut extractor, &mut accumulator);
        result.map_err(|source| {
            constructors::reference_parse_error(role, path, reader.error_kind, source)
        })?;
    }
    // In particular, release the parser, decompressor, flattened sequence, and
    // extractor's largest-record buffers before Deacon allocates its hash set.
    drop(extractor);
    let values = accumulator.finish();
    if values.is_empty() {
        return Err(MalformedInputError::EmptyReferenceIndex {
            path: path.to_path_buf(),
            k,
            w,
        }
        .into());
    }
    Index::from_minimizers(k, w, values).map_err(|source| {
        InternalError::IndexConstruction {
            source: source.into_boxed_dyn_error(),
        }
        .into()
    })
}

#[cfg(test)]
mod tests {
    use std::{
        io::{self, Read},
        path::Path,
        process::ExitCode,
    };

    use deacon::{IndexKind, MinimizerVec};

    use super::{Accumulator, construct_reference};

    fn check_accumulator<T: Copy + Ord + std::fmt::Debug + From<u16>>() {
        let mut acc = Accumulator::<T>::new();
        let repeated = vec![T::from(7); 50_000];
        acc.append(&repeated);
        assert!(acc.values.len() < 10_000);
        let novel: Vec<T> = (0..20_000).map(T::from).collect();
        acc.append(&novel);
        acc.append(&repeated);
        acc.compact();
        assert_eq!(acc.values, novel);
    }

    #[test]
    fn compaction_handles_repeated_and_novel_keys_at_both_widths() {
        check_accumulator::<u64>();
        check_accumulator::<u128>();
    }

    #[test]
    fn boundaries_and_multiline_sequences_are_preserved() {
        let path = Path::new("reference.fa");
        let separate = construct_reference(&b">a\nAA\n>b\nCC\n"[..], "reference", path, 3, 1);
        assert!(matches!(
            separate,
            Err(crate::error::RunError::MalformedInput(_))
        ));
        let joined = construct_reference(&b">a\r\nAA\r\nCC\r\n"[..], "reference", path, 3, 1)
            .expect("multiline reference");
        assert_eq!(joined.kind(), IndexKind::Exact);
        assert!(!joined.is_empty());
        assert_eq!((joined.kmer_length(), joined.window_size()), (3, 1));
    }

    #[test]
    fn malformed_fastq_and_no_minimizers_are_rejected() {
        let path = Path::new("reference.fq");
        for input in [&b"@a\nACGT\n+\n!\n"[..], &b">a\nNNNNNN\n"[..]] {
            let error =
                construct_reference(input, "reference", path, 3, 1).expect_err("invalid reference");
            assert_eq!(error.exit_code(), ExitCode::from(65));
        }
    }

    struct FailingReader {
        data: &'static [u8],
        interrupted: bool,
        fail: bool,
    }

    impl Read for FailingReader {
        fn read(&mut self, output: &mut [u8]) -> io::Result<usize> {
            if self.interrupted {
                self.interrupted = false;
                return Err(io::ErrorKind::Interrupted.into());
            }
            if self.fail && self.data.is_empty() {
                return Err(io::ErrorKind::PermissionDenied.into());
            }
            self.data.read(output)
        }
    }

    #[test]
    fn read_failures_and_retried_interruptions_are_classified() {
        let path = Path::new("reference.fa");
        for data in [&b""[..], &b">a\nACGTACGT\n"[..]] {
            let reader = FailingReader {
                data,
                interrupted: false,
                fail: true,
            };
            let error =
                construct_reference(reader, "reference", path, 3, 1).expect_err("read failure");
            assert_eq!(error.exit_code(), ExitCode::from(74));
        }
        let reader = FailingReader {
            data: b"@a\nACGT\n+\n!\n",
            interrupted: true,
            fail: false,
        };
        let error =
            construct_reference(reader, "reference", path, 3, 1).expect_err("malformed reference");
        assert_eq!(error.exit_code(), ExitCode::from(65));
        let reader = FailingReader {
            data: b">a\nACGTACGT\n",
            interrupted: true,
            fail: false,
        };
        assert!(construct_reference(reader, "reference", path, 3, 1).is_ok());
    }

    #[test]
    fn packed_width_selection_matches_deacon() {
        let mut acc = super::PackedAccumulator::new(33);
        acc.append(&MinimizerVec::U128(vec![1, 1, 2]));
        assert!(matches!(acc.finish(), MinimizerVec::U128(values) if values == [1, 2]));
    }
}
