//! Quality-based trimming and quality-score binning transforms.

use std::{fmt, str::FromStr};

use crate::{
    plan::{
        BuildPlan, IntoExecutionStep, ReadTransform, TransformArena, TransformResult, TransformStep,
    },
    record::RecordView,
};

/// Supported numbers of quality-score bins for lossy FASTQ quality quantization.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum QualityBinCount {
    /// Collapse quality scores into two bins.
    Two,
    /// Collapse quality scores into three bins.
    Three,
    /// Collapse quality scores into five bins.
    Five,
}

impl fmt::Display for QualityBinCount {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::Two => "2",
            Self::Three => "3",
            Self::Five => "5",
        })
    }
}

impl FromStr for QualityBinCount {
    type Err = String;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        match value {
            "2" => Ok(Self::Two),
            "3" => Ok(Self::Three),
            "5" => Ok(Self::Five),
            _ => Err("quality bin count must be one of 2, 3, or 5".to_owned()),
        }
    }
}

/// Read transform that lossily remaps Phred+33 quality bytes into fewer bins.
pub(crate) struct QualityBin {
    boundaries: &'static [u8],
    representatives: &'static [u8],
}

impl QualityBin {
    /// Construct a quality binner using platform-neutral Phred+33 bins.
    pub(crate) const fn new(count: QualityBinCount) -> Self {
        match count {
            QualityBinCount::Two => Self {
                boundaries: &[33 + 20],
                representatives: &[33 + 10, 33 + 30],
            },
            QualityBinCount::Three => Self {
                boundaries: &[33 + 13, 33 + 27],
                representatives: &[33 + 6, 33 + 20, 33 + 33],
            },
            QualityBinCount::Five => Self {
                boundaries: &[33 + 8, 33 + 16, 33 + 24, 33 + 32],
                representatives: &[33 + 4, 33 + 12, 33 + 20, 33 + 28, 33 + 36],
            },
        }
    }

    fn bin_byte(&self, byte: u8) -> u8 {
        for (index, &boundary) in self.boundaries.iter().enumerate() {
            if byte < boundary {
                return self.representatives[index];
            }
        }

        self.representatives[self.boundaries.len()]
    }
}

impl ReadTransform for QualityBin {
    fn code(&self) -> &'static str {
        "quality_bin"
    }

    fn apply<'a>(&self, record: RecordView<'a>, arena: &'a TransformArena) -> TransformResult<'a> {
        assert_eq!(
            record.sequence().len(),
            record.quality().len(),
            "quality binning requires equal sequence and quality lengths"
        );

        let binned_quality = arena.alloc_slice_copy_mut(record.quality());
        let mut changed = false;

        for byte in binned_quality.iter_mut() {
            let binned = self.bin_byte(*byte);
            changed |= binned != *byte;
            *byte = binned;
        }

        if !changed {
            return TransformResult {
                record,
                applied: false,
            };
        }

        TransformResult {
            record: record
                .with_sequence_and_quality(record.sequence(), binned_quality)
                .expect("quality binning should preserve equal sequence and quality lengths"),
            applied: true,
        }
    }
}

impl IntoExecutionStep for QualityBin {
    fn into_execution_step(self) -> Box<dyn crate::plan::ExecutionStep> {
        Box::new(TransformStep(self))
    }
}

/// Read transform that applies 3' quality trimming using a Phred cutoff.
pub(crate) struct QualityTrim {
    cutoff: u8,
}

impl QualityTrim {
    /// Construct a new quality trimmer with the provided Phred cutoff.
    pub(crate) const fn new(cutoff: u8) -> Self {
        Self { cutoff }
    }

    /// Return the exclusive trim end selected by the 3' quality-trimming algorithm.
    fn quality_trim_end(&self, quality: &[u8]) -> usize {
        let mut best_sum = 0_i32;
        let mut running_sum = 0_i32;
        let mut trim_end = quality.len();

        for (idx, &ascii_q) in quality.iter().enumerate().rev() {
            let q = i32::from(ascii_q.saturating_sub(33));
            let score = i32::from(self.cutoff) - q;
            running_sum += score;

            if running_sum < 0 {
                running_sum = 0;
            } else if running_sum > best_sum {
                best_sum = running_sum;
                trim_end = idx;
            }
        }

        trim_end
    }
}

impl ReadTransform for QualityTrim {
    fn code(&self) -> &'static str {
        "quality_trim"
    }

    fn apply<'a>(&self, record: RecordView<'a>, _arena: &'a TransformArena) -> TransformResult<'a> {
        assert_eq!(
            record.sequence().len(),
            record.quality().len(),
            "quality trimming requires equal sequence and quality lengths"
        );

        let trim_end = self.quality_trim_end(record.quality());

        if trim_end == record.sequence().len() {
            return TransformResult {
                record,
                applied: false,
            };
        }

        TransformResult {
            record: record
                .with_sequence_and_quality(
                    &record.sequence()[..trim_end],
                    &record.quality()[..trim_end],
                )
                .expect("quality trimming should preserve equal sequence and quality lengths"),
            applied: true,
        }
    }
}

impl IntoExecutionStep for QualityTrim {
    fn into_execution_step(self) -> Box<dyn crate::plan::ExecutionStep> {
        Box::new(TransformStep(self))
    }
}

/// Fluent extension trait adding the `.quality_trim(...)` transform combinator to plans.
pub(crate) trait QualityTrimTransform: BuildPlan {
    /// Trim low-quality 3' sequence using the configured Phred cutoff.
    fn quality_trim(self, cutoff: u8) -> Self {
        self.step(QualityTrim::new(cutoff))
    }
}

impl<T> QualityTrimTransform for T where T: BuildPlan {}

/// Fluent extension trait adding the `.quality_bin(...)` transform combinator to plans.
pub(crate) trait QualityBinTransform: BuildPlan {
    /// Lossily remap Phred+33 quality bytes into the configured number of bins.
    fn quality_bin(self, count: QualityBinCount) -> Self {
        self.step(QualityBin::new(count))
    }
}

impl<T> QualityBinTransform for T where T: BuildPlan {}

#[cfg(test)]
mod tests {
    use crate::{
        plan::{ReadTransform, TransformArena},
        record::RecordView,
    };

    use super::{QualityBin, QualityBinCount, QualityTrim};

    fn record(sequence: &'static [u8], quality: &'static [u8]) -> RecordView<'static> {
        RecordView::new(b"read1", sequence, quality)
    }

    #[test]
    fn quality_trim_leaves_high_quality_read_unchanged() {
        let trimmer = QualityTrim::new(20);
        let arena = TransformArena::new();
        let result = trimmer.apply(record(b"ACGT", b"IIII"), &arena);

        assert!(!result.applied);
        assert_eq!(result.record.sequence(), b"ACGT");
    }

    #[test]
    fn quality_trim_removes_low_quality_suffix() {
        let trimmer = QualityTrim::new(20);
        let arena = TransformArena::new();
        let result = trimmer.apply(record(b"ACGTAC", b"IIIII!"), &arena);

        assert!(result.applied);
        assert_eq!(result.record.sequence(), b"ACGTA");
        assert_eq!(result.record.quality(), b"IIIII");
    }

    #[test]
    fn quality_trim_cuts_at_low_quality_tail_start_even_if_tail_recovers() {
        let trimmer = QualityTrim::new(20);
        let arena = TransformArena::new();
        let result = trimmer.apply(record(b"ACGT", b"I!II"), &arena);

        assert!(result.applied);
        assert_eq!(result.record.sequence(), b"A");
        assert_eq!(result.record.quality(), b"I");
    }

    #[test]
    fn quality_trim_can_trim_entire_read_without_rejecting() {
        let trimmer = QualityTrim::new(20);
        let arena = TransformArena::new();
        let result = trimmer.apply(record(b"ACGT", b"!!!!"), &arena);

        assert!(result.applied);
        assert_eq!(result.record.sequence(), b"");
        assert_eq!(result.record.quality(), b"");
    }

    #[test]
    fn quality_bin_two_bins_maps_exact_boundary_to_upper_bin() {
        let binner = QualityBin::new(QualityBinCount::Two);
        let arena = TransformArena::new();
        let result = binner.apply(record(b"ACG", b"!5I"), &arena);

        assert!(result.applied);
        assert_eq!(result.record.sequence(), b"ACG");
        assert_eq!(result.record.quality(), b"+??");
    }

    #[test]
    fn quality_bin_three_bins_uses_platform_neutral_boundaries() {
        let binner = QualityBin::new(QualityBinCount::Three);
        let arena = TransformArena::new();
        let result = binner.apply(record(b"ACG", b"!5I"), &arena);

        assert!(result.applied);
        assert_eq!(result.record.quality(), b"'5B");
    }

    #[test]
    fn quality_bin_five_bins_preserves_already_representative_quality() {
        let binner = QualityBin::new(QualityBinCount::Five);
        let arena = TransformArena::new();
        let result = binner.apply(record(b"ACGTA", b"%-5=E"), &arena);

        assert!(!result.applied);
        assert_eq!(result.record.quality(), b"%-5=E");
    }

    #[test]
    fn quality_bin_count_rejects_unsupported_counts() {
        assert_eq!("2".parse(), Ok(QualityBinCount::Two));
        assert_eq!("3".parse(), Ok(QualityBinCount::Three));
        assert_eq!("5".parse(), Ok(QualityBinCount::Five));
        assert!("4".parse::<QualityBinCount>().is_err());
    }
}
