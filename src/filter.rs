//! Core read filtering operations and their fluent plan extensions.

use bytecount::count;
use std::f64::consts::LN_2;

use crate::{
    plan::{BuildPlan, FilterStep, IntoExecutionStep, ReadFilter, RejectionReason},
    record::RecordView,
};

/// Rejection reason emitted when a record contains too many ambiguous `N` bases.
#[derive(Debug)]
pub(crate) struct TooManyNs {}

impl RejectionReason for TooManyNs {
    fn code(&self) -> &'static str {
        "too_many_ns"
    }
}

/// Zero-copy filter that rejects records containing more than `max_ns` ambiguous bases.
pub(crate) struct MaxNs {
    max_ns: usize,
}

impl MaxNs {
    pub(crate) const fn new(max_ns: usize) -> Self {
        Self { max_ns }
    }
}

impl ReadFilter for MaxNs {
    type Reason = TooManyNs;

    fn evaluate(&self, record: &RecordView<'_>) -> Result<(), Self::Reason> {
        let observed = count(record.sequence(), b'N');
        if observed > self.max_ns {
            Err(TooManyNs {})
        } else {
            Ok(())
        }
    }
}

impl IntoExecutionStep for MaxNs {
    fn into_execution_step(self) -> Box<dyn crate::plan::ExecutionStep> {
        Box::new(FilterStep(self))
    }
}

/// Fluent extension trait adding the `.max_ns(...)` filter combinator to plans.
pub(crate) trait MaxNsFilter: BuildPlan {
    fn max_ns(self, max_ns: usize) -> Self {
        self.step(MaxNs::new(max_ns))
    }
}

impl<T> MaxNsFilter for T where T: BuildPlan {}

/// Rejection reason emitted when a record is shorter than the configured minimum length.
#[derive(Debug)]
pub(crate) struct TooShort {}

impl RejectionReason for TooShort {
    fn code(&self) -> &'static str {
        "too_short"
    }
}

/// Zero-copy filter that rejects records shorter than `min_length`.
pub(crate) struct MinLength {
    min_length: usize,
}

impl MinLength {
    pub(crate) const fn new(min_length: usize) -> Self {
        Self { min_length }
    }
}

impl ReadFilter for MinLength {
    type Reason = TooShort;

    fn evaluate(&self, record: &RecordView<'_>) -> Result<(), Self::Reason> {
        let observed = record.sequence().len();
        if observed < self.min_length {
            Err(TooShort {})
        } else {
            Ok(())
        }
    }
}

impl IntoExecutionStep for MinLength {
    fn into_execution_step(self) -> Box<dyn crate::plan::ExecutionStep> {
        Box::new(FilterStep(self))
    }
}

/// Fluent extension trait adding the `.min_length(...)` filter combinator to plans.
pub(crate) trait MinLengthFilter: BuildPlan {
    fn min_length(self, min_length: usize) -> Self {
        self.step(MinLength::new(min_length))
    }
}

impl<T> MinLengthFilter for T where T: BuildPlan {}

/// Rejection reason emitted when a record's mean quality falls below the threshold.
#[derive(Debug)]
pub(crate) struct LowMeanQuality {}

impl RejectionReason for LowMeanQuality {
    fn code(&self) -> &'static str {
        "low_mean_quality"
    }
}

/// Zero-copy filter that rejects records whose mean Phred quality is too low.
pub(crate) struct MinMeanQuality {
    min_mean_q: f64,
}

impl MinMeanQuality {
    pub(crate) const fn new(min_mean_q: f64) -> Self {
        Self { min_mean_q }
    }
}

impl ReadFilter for MinMeanQuality {
    type Reason = LowMeanQuality;

    fn evaluate(&self, record: &RecordView<'_>) -> Result<(), Self::Reason> {
        let total: u32 = record
            .quality()
            .iter()
            .map(|q| u32::from(q.saturating_sub(33)))
            .sum();
        let observed = f64::from(total)
            / f64::from(
                u32::try_from(record.quality().len()).expect("quality length must fit in u32"),
            );

        if observed < self.min_mean_q {
            Err(LowMeanQuality {})
        } else {
            Ok(())
        }
    }
}

impl IntoExecutionStep for MinMeanQuality {
    fn into_execution_step(self) -> Box<dyn crate::plan::ExecutionStep> {
        Box::new(FilterStep(self))
    }
}

/// Fluent extension trait adding the `.min_mean_q(...)` filter combinator to plans.
pub(crate) trait MinMeanQualityFilter: BuildPlan {
    fn min_mean_q(self, min_mean_q: f64) -> Self {
        self.step(MinMeanQuality::new(min_mean_q))
    }
}

impl<T> MinMeanQualityFilter for T where T: BuildPlan {}

/// Rejection reason emitted when a record's nucleotide composition is too low-entropy.
#[derive(Debug)]
pub(crate) struct LowEntropy {}

impl RejectionReason for LowEntropy {
    fn code(&self) -> &'static str {
        "low_entropy"
    }
}

/// Zero-copy filter that rejects records whose Shannon entropy is below the threshold.
pub(crate) struct MinEntropy {
    min_entropy: f64,
}

impl MinEntropy {
    /// Construct a new entropy filter with the provided minimum Shannon entropy.
    pub(crate) const fn new(min_entropy: f64) -> Self {
        Self { min_entropy }
    }
}

impl ReadFilter for MinEntropy {
    type Reason = LowEntropy;

    fn evaluate(&self, record: &RecordView<'_>) -> Result<(), Self::Reason> {
        let observed = shannon_entropy(record.sequence());
        if observed < self.min_entropy {
            Err(LowEntropy {})
        } else {
            Ok(())
        }
    }
}

impl IntoExecutionStep for MinEntropy {
    fn into_execution_step(self) -> Box<dyn crate::plan::ExecutionStep> {
        Box::new(FilterStep(self))
    }
}

/// Fluent extension trait adding the `.min_entropy(...)` filter combinator to plans.
pub(crate) trait MinEntropyFilter: BuildPlan {
    /// Reject reads whose post-transform Shannon entropy falls below `min_entropy`.
    fn min_entropy(self, min_entropy: f64) -> Self {
        self.step(MinEntropy::new(min_entropy))
    }
}

impl<T> MinEntropyFilter for T where T: BuildPlan {}

/// Compute Shannon entropy over the canonical DNA alphabet `A/C/G/T`, ignoring other symbols.
fn shannon_entropy(sequence: &[u8]) -> f64 {
    let mut counts = [0_u32; 4];

    for &base in sequence {
        match base.to_ascii_uppercase() {
            b'A' => counts[0] += 1,
            b'C' => counts[1] += 1,
            b'G' => counts[2] += 1,
            b'T' => counts[3] += 1,
            _ => {}
        }
    }

    let total = counts.iter().sum::<u32>();
    if total == 0 {
        return 0.0;
    }

    counts
        .into_iter()
        .filter(|count| *count > 0)
        .map(|count| {
            let p = f64::from(count) / f64::from(total);
            -(p * (p.ln() / LN_2))
        })
        .sum()
}

#[cfg(test)]
mod tests {
    use super::shannon_entropy;

    #[test]
    fn entropy_is_zero_for_homopolymer() {
        assert!((shannon_entropy(b"AAAAAA") - 0.0).abs() < f64::EPSILON);
    }

    #[test]
    fn entropy_is_two_for_balanced_bases() {
        let observed = shannon_entropy(b"ACGT");
        assert!((observed - 2.0).abs() < f64::EPSILON);
    }

    #[test]
    fn entropy_ignores_non_acgt_symbols() {
        let observed = shannon_entropy(b"AANNCC");
        assert!((observed - 1.0).abs() < f64::EPSILON);
    }
}
