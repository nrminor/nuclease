//! Quality-based trimming transforms.

use crate::{
    plan::{
        BuildPlan, IntoExecutionStep, ReadTransform, TransformArena, TransformResult, TransformStep,
    },
    record::RecordView,
};

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

#[cfg(test)]
mod tests {
    use crate::{
        plan::{ReadTransform, TransformArena},
        record::RecordView,
    };

    use super::QualityTrim;

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
}
