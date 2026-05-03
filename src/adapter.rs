//! Adapter trimming transforms and curated adapter catalogs.

use wide::u8x16;

use crate::{
    plan::{
        BuildPlan, IntoExecutionStep, ReadTransform, TransformArena, TransformResult, TransformStep,
    },
    record::RecordView,
};

/// Curated adapter catalogs available to trimming transforms.
pub(crate) enum AdapterCatalog {
    IlluminaTruSeq,
}

impl AdapterCatalog {
    /// Construct the curated Illumina `TruSeq` adapter preset.
    pub(crate) const fn illumina_truseq() -> Self {
        Self::IlluminaTruSeq
    }

    /// Return the raw adapter sequences included in this preset.
    fn adapters(&self) -> &'static [&'static [u8]] {
        match self {
            Self::IlluminaTruSeq => &[TRUSEQ_R1, TRUSEQ_R2],
        }
    }

    /// Build the default overlap matcher used for this preset.
    fn default_matcher(&self) -> AdapterMatcher<'_> {
        AdapterMatcher {
            adapters: self.adapters(),
            min_overlap: 8,
            max_mismatch_numerator: 1,
            max_mismatch_denominator: 8,
        }
    }
}

const TRUSEQ_R1: &[u8] = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
const TRUSEQ_R2: &[u8] = b"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT";
/// SIMD chunk width used by the bytewise adapter mismatch counter.
const SIMD_LANES: usize = 16;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct AdapterMatch {
    trim_end: usize,
    overlap: usize,
    mismatches: usize,
}

impl AdapterMatch {
    fn is_better_than(&self, other: &Self) -> bool {
        self.overlap > other.overlap
            || (self.overlap == other.overlap && self.mismatches < other.mismatches)
    }
}

struct AdapterMatcher<'a> {
    adapters: &'a [&'a [u8]],
    min_overlap: usize,
    max_mismatch_numerator: usize,
    max_mismatch_denominator: usize,
}

impl AdapterMatcher<'_> {
    /// Return the exclusive trim end chosen by the best valid adapter overlap, if any.
    fn find_trim_end(&self, read: &[u8]) -> Option<usize> {
        self.best_match(read).map(|candidate| candidate.trim_end)
    }

    /// Return the best valid adapter overlap across all adapters in this matcher.
    fn best_match(&self, read: &[u8]) -> Option<AdapterMatch> {
        let mut best: Option<AdapterMatch> = None;

        for adapter in self.adapters {
            let Some(candidate) = self.best_match_for_adapter(read, adapter) else {
                continue;
            };

            if best
                .as_ref()
                .is_none_or(|current| candidate.is_better_than(current))
            {
                best = Some(candidate);
            }
        }

        best
    }

    /// Return the best valid overlap between one read and one adapter sequence.
    fn best_match_for_adapter(&self, read: &[u8], adapter: &[u8]) -> Option<AdapterMatch> {
        let max_overlap = read.len().min(adapter.len());
        let mut best: Option<AdapterMatch> = None;

        for overlap in self.min_overlap..=max_overlap {
            let read_suffix = &read[read.len() - overlap..];
            let adapter_prefix = &adapter[..overlap];
            let max_mismatches = self.max_mismatches(overlap);

            let Some(mismatches) =
                Self::mismatch_count_bounded(read_suffix, adapter_prefix, max_mismatches)
            else {
                continue;
            };

            let candidate = AdapterMatch {
                trim_end: read.len() - overlap,
                overlap,
                mismatches,
            };

            if best
                .as_ref()
                .is_none_or(|current| candidate.is_better_than(current))
            {
                best = Some(candidate);
            }
        }

        best
    }

    /// Convert the configured mismatch ratio into an absolute mismatch budget for one overlap.
    fn max_mismatches(&self, overlap: usize) -> usize {
        overlap.saturating_mul(self.max_mismatch_numerator) / self.max_mismatch_denominator
    }

    /// Count mismatches up to a fixed budget, dispatching to SIMD when the overlap is wide enough.
    fn mismatch_count_bounded(left: &[u8], right: &[u8], max_mismatches: usize) -> Option<usize> {
        debug_assert_eq!(left.len(), right.len());

        if left.len() >= 16 {
            Self::mismatch_count_bounded_simd(left, right, max_mismatches)
        } else {
            Self::mismatch_count_bounded_scalar(left, right, max_mismatches)
        }
    }

    /// Scalar mismatch counter used for short overlaps and as the SIMD tail fallback.
    fn mismatch_count_bounded_scalar(
        left: &[u8],
        right: &[u8],
        max_mismatches: usize,
    ) -> Option<usize> {
        let mut mismatches = 0_usize;

        for (&lhs, &rhs) in left.iter().zip(right) {
            if lhs != rhs {
                mismatches += 1;
                if mismatches > max_mismatches {
                    return None;
                }
            }
        }

        Some(mismatches)
    }

    /// SIMD-accelerated mismatch counter used for adapter overlap scoring.
    fn mismatch_count_bounded_simd(
        left: &[u8],
        right: &[u8],
        max_mismatches: usize,
    ) -> Option<usize> {
        let mut mismatches = 0_usize;
        let chunks = left.len() / SIMD_LANES;
        let tail_start = chunks * SIMD_LANES;

        for chunk_idx in 0..chunks {
            let start = chunk_idx * SIMD_LANES;
            let lhs = u8x16::from([
                left[start],
                left[start + 1],
                left[start + 2],
                left[start + 3],
                left[start + 4],
                left[start + 5],
                left[start + 6],
                left[start + 7],
                left[start + 8],
                left[start + 9],
                left[start + 10],
                left[start + 11],
                left[start + 12],
                left[start + 13],
                left[start + 14],
                left[start + 15],
            ]);
            let rhs = u8x16::from([
                right[start],
                right[start + 1],
                right[start + 2],
                right[start + 3],
                right[start + 4],
                right[start + 5],
                right[start + 6],
                right[start + 7],
                right[start + 8],
                right[start + 9],
                right[start + 10],
                right[start + 11],
                right[start + 12],
                right[start + 13],
                right[start + 14],
                right[start + 15],
            ]);

            let eq = lhs.simd_eq(rhs).to_array();
            let equal_count = eq.into_iter().filter(|lane| *lane == u8::MAX).count();
            mismatches += SIMD_LANES - equal_count;

            if mismatches > max_mismatches {
                return None;
            }
        }

        let tail_mismatches = Self::mismatch_count_bounded_scalar(
            &left[tail_start..],
            &right[tail_start..],
            max_mismatches.saturating_sub(mismatches),
        )?;

        Some(mismatches + tail_mismatches)
    }
}

/// Read transform that trims known adapter sequence from the end of a read.
pub(crate) struct AdapterTrim {
    catalog: AdapterCatalog,
}

impl AdapterTrim {
    /// Construct a new adapter trimmer backed by the provided curated catalog.
    pub(crate) const fn new(catalog: AdapterCatalog) -> Self {
        Self { catalog }
    }

    fn default_matcher(&self) -> AdapterMatcher<'_> {
        self.catalog.default_matcher()
    }
}

impl ReadTransform for AdapterTrim {
    fn code(&self) -> &'static str {
        "trim_adapters"
    }

    fn apply<'a>(&self, record: RecordView<'a>, arena: &'a TransformArena) -> TransformResult<'a> {
        assert_eq!(
            record.sequence().len(),
            record.quality().len(),
            "adapter trimming requires equal sequence and quality lengths"
        );

        let matcher = self.default_matcher();
        let Some(trim_end) = matcher.find_trim_end(record.sequence()) else {
            return TransformResult {
                record,
                applied: false,
            };
        };

        TransformResult {
            record: record
                .with_sequence_and_quality(
                    arena.alloc_slice_copy(&record.sequence()[..trim_end]),
                    arena.alloc_slice_copy(&record.quality()[..trim_end]),
                )
                .expect("adapter trimming should preserve equal sequence and quality lengths"),
            applied: true,
        }
    }
}

impl IntoExecutionStep for AdapterTrim {
    fn into_execution_step(self) -> Box<dyn crate::plan::ExecutionStep> {
        Box::new(TransformStep(self))
    }
}

/// Fluent extension trait adding the `.trim_adapters(...)` transform combinator to plans.
pub(crate) trait TrimAdaptersTransform: BuildPlan {
    /// Trim adapter overlap using the matching strategy associated with the selected catalog.
    fn trim_adapters(self, catalog: AdapterCatalog) -> Self {
        self.step(AdapterTrim::new(catalog))
    }
}

impl<T> TrimAdaptersTransform for T where T: BuildPlan {}

#[cfg(test)]
mod tests {
    use super::{AdapterCatalog, AdapterMatcher};

    #[test]
    fn overlap_match_prefers_longest_valid_candidate() {
        let matcher = {
            let adapters: &[&[u8]] = &[b"AGATCGGAAGAG"];
            AdapterMatcher {
                adapters,
                min_overlap: 4,
                max_mismatch_numerator: 1,
                max_mismatch_denominator: 8,
            }
        };
        let candidate = matcher
            .best_match(b"TTTTAGATCGGA")
            .expect("overlap match should be found");

        assert_eq!(candidate.overlap, 8);
        assert_eq!(candidate.trim_end, 4);
    }

    #[test]
    fn overlap_match_rejects_short_overlap() {
        let matcher = {
            let adapters: &[&[u8]] = &[b"AGATCGGAAGAG"];
            AdapterMatcher {
                adapters,
                min_overlap: 8,
                max_mismatch_numerator: 1,
                max_mismatch_denominator: 8,
            }
        };
        assert!(matcher.best_match(b"TTTTAGAT").is_none());
    }

    #[test]
    fn overlap_match_accepts_bounded_mismatches() {
        let matcher = {
            let adapters: &[&[u8]] = &[b"AGATCGGAAGAG"];
            AdapterMatcher {
                adapters,
                min_overlap: 8,
                max_mismatch_numerator: 1,
                max_mismatch_denominator: 4,
            }
        };
        let candidate = matcher
            .best_match(b"TTTTAGATCGTA")
            .expect("bounded mismatch overlap should be found");

        assert_eq!(candidate.trim_end, 4);
        assert_eq!(candidate.mismatches, 1);
    }

    #[test]
    fn illumina_catalog_exposes_matcher() {
        let trim_end = AdapterCatalog::illumina_truseq()
            .default_matcher()
            .find_trim_end(b"ACGTAGATCGGAAG");
        assert_eq!(trim_end, Some(4));
    }
}
