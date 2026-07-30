//! Paired-read merging backed by `libpairassembly`.

use std::str;

use libpairassembly::{Assembler, CorrectionParams, OverlapParams, PairInput, SeqRecordView};

use crate::{
    error::{IndeterminateInputError, InternalError, MalformedInputError, Result, RunError},
    plan::{
        BuildPlan, IntoExecutionStep, PairTransform, PairTransformResult, PairTransformStep,
        RecordPair, TransformArena,
    },
    record::{InputSource, RecordProvenance, RecordView},
};

/// Pair-aware transform that attempts to merge paired-end reads.
pub(crate) struct MergePairs {
    assembler: Assembler,
}

/// User-facing paired-read merge tuning carried from the CLI into the assembler.
#[derive(Clone, Copy, Debug)]
pub(crate) struct MergePairsConfig {
    pub(crate) min_overlap: usize,
    pub(crate) max_mismatch_rate: f32,
    pub(crate) min_correction_delta_q: u8,
}

impl Default for MergePairsConfig {
    fn default() -> Self {
        let overlap = OverlapParams::default();
        let correction = CorrectionParams::default();
        Self {
            min_overlap: overlap.min_overlap(),
            max_mismatch_rate: overlap.diff_percent_max(),
            min_correction_delta_q: correction.min_base_correction_delta_q,
        }
    }
}

struct AssemblyRecordView<'a> {
    id: &'a str,
    sequence: &'a str,
    quality: &'a str,
}

impl MergePairs {
    /// Construct a merge transform with the provided merge settings.
    pub(crate) fn new(config: MergePairsConfig) -> Result<Self> {
        let overlap_params = OverlapParams::default()
            .with_min_overlap(config.min_overlap)
            .with_diff_percent_max(config.max_mismatch_rate);
        let correction_params = CorrectionParams::default()
            .with_min_base_correction_delta_q(config.min_correction_delta_q);
        let assembler = Assembler::builder()
            .with_overlap_params(overlap_params)
            .with_correction_params(correction_params)
            .build()
            .map_err(|source| InternalError::PairAssembly { source })?;
        Ok(Self { assembler })
    }
}

impl PairTransform for MergePairs {
    fn code(&self) -> &'static str {
        "merge_pairs"
    }

    fn apply_pair<'a>(
        &mut self,
        pair: RecordPair<'a>,
        arena: &'a TransformArena,
    ) -> Result<PairTransformResult<'a>> {
        let input = PairInput::new(
            AssemblyRecordView::try_from_record(pair.left)?,
            AssemblyRecordView::try_from_record(pair.right)?,
        );

        let Some(merged) = self
            .assembler
            .process_pair(&input)
            .map_err(|source| InternalError::PairAssembly { source })?
        else {
            return Ok(PairTransformResult::Pair {
                pair,
                applied: false,
            });
        };

        let provenance = pair.left.provenance().map(|provenance| RecordProvenance {
            source: provenance.source,
            mate: None,
        });
        let mut record = RecordView::new(
            arena.alloc_slice_copy(merged.id().as_bytes()),
            arena.alloc_slice_copy(merged.sequence_bytes()),
            arena.alloc_slice_copy(merged.quality_bytes()),
        );
        if let Some(provenance) = provenance {
            record = record.with_provenance(provenance);
        }

        Ok(PairTransformResult::Single {
            record,
            applied: true,
        })
    }
}

impl IntoExecutionStep for MergePairs {
    fn into_execution_step(self) -> Box<dyn crate::plan::ExecutionStep> {
        Box::new(PairTransformStep(self))
    }
}

/// Fluent extension trait adding the `.merge_pairs(...)` transform combinator to plans.
pub(crate) trait MergePairsTransform: BuildPlan {
    /// Attempt to merge paired reads before downstream per-record filtering.
    fn merge_pairs(self, config: MergePairsConfig) -> Result<Self> {
        Ok(self.step(MergePairs::new(config)?))
    }
}

impl<T> MergePairsTransform for T where T: BuildPlan {}

impl<'a> AssemblyRecordView<'a> {
    fn try_from_record(record: RecordView<'a>) -> Result<Self> {
        Ok(Self {
            id: str::from_utf8(record.pair_key())
                .map_err(|source| invalid_utf8_error(record, "read identifier", source))?,
            sequence: str::from_utf8(record.sequence())
                .map_err(|source| invalid_utf8_error(record, "read sequence", source))?,
            quality: str::from_utf8(record.quality())
                .map_err(|source| invalid_utf8_error(record, "read quality", source))?,
        })
    }
}

fn invalid_utf8_error(
    record: RecordView<'_>,
    field: &'static str,
    source: std::str::Utf8Error,
) -> RunError {
    match record.provenance().map(|provenance| provenance.source) {
        Some(InputSource::Ena { accession }) => IndeterminateInputError::InvalidUtf8 {
            accession: accession.to_owned(),
            field,
            source,
        }
        .into(),
        _ => MalformedInputError::InvalidUtf8 { field, source }.into(),
    }
}

impl SeqRecordView for AssemblyRecordView<'_> {
    fn id(&self) -> &str {
        self.id
    }

    fn seq(&self) -> &str {
        self.sequence
    }

    fn qual(&self) -> &str {
        self.quality
    }
}

#[cfg(test)]
mod tests {
    use crate::plan::{PairTransform, TransformArena};

    use super::*;

    fn pair() -> RecordPair<'static> {
        RecordPair {
            left: RecordView::new(
                b"read-1/1",
                b"ACGTTGCAGTACGATCGTACGGAATTCGCCGATGACTGACCTAGGTCAGTACGATC",
                b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
            ),
            right: RecordView::new(
                b"read-1/2",
                b"GATCGTACTGACCTAGGTCAGTCATCGGCGAATTCCGTACGATCGTACTGCAACGT",
                b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
            ),
        }
    }

    #[test]
    fn merge_pairs_returns_single_merged_record() {
        let mut merge = MergePairs::new(MergePairsConfig::default())
            .expect("default libpairassembly configuration should build");
        let arena = TransformArena::new();

        let result = merge
            .apply_pair(pair(), &arena)
            .expect("fixture pair should be valid libpairassembly input");

        let PairTransformResult::Single { record, applied } = result else {
            panic!("overlapping paired-end fixture should merge into one record");
        };
        assert!(applied);
        assert_eq!(record.header(), b"read-1");
        assert_eq!(record.sequence().len(), record.quality().len());
    }

    #[test]
    fn invalid_utf8_classification_respects_record_origin() {
        let local = RecordView::new(&[0xff], b"ACGT", b"IIII");
        let Err(local_error) = AssemblyRecordView::try_from_record(local) else {
            panic!("invalid local identifier should fail UTF-8 conversion");
        };
        assert!(matches!(
            local_error,
            RunError::MalformedInput(error)
                if matches!(*error, MalformedInputError::InvalidUtf8 { .. })
        ));

        let ena = RecordView::new(&[0xff], b"ACGT", b"IIII").with_provenance(RecordProvenance {
            source: InputSource::Ena {
                accession: "SRR35939766",
            },
            mate: Some(crate::record::MateSide::Left),
        });
        let Err(ena_error) = AssemblyRecordView::try_from_record(ena) else {
            panic!("invalid ENA identifier should fail UTF-8 conversion");
        };
        assert!(matches!(
            ena_error,
            RunError::IndeterminateInput(IndeterminateInputError::InvalidUtf8 { .. })
        ));
    }
}
