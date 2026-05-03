//! Logical authoring and execution substrate for streaming preprocessing plans.

use std::marker::PhantomData;

use bumpalo::Bump;

use color_eyre::eyre::Result;

use crate::record::{ReadStats, RecordView};

/// Typestate marker for a logical preprocessing plan still being authored.
pub(crate) struct Logical;

/// Typestate marker for a compiled execution plan ready to run.
pub(crate) struct Execution;

/// Policy for handling paired reads when exactly one mate survives per-read execution.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub(crate) enum OrphanPolicy {
    /// Drop the entire pair when either mate fails.
    #[default]
    DropPair,
    /// Emit the surviving mate as an orphan.
    EmitOrphan,
}

/// Scratch arena reused across execution of one single record or one pair.
///
/// The same arena is reset after an execution unit has been fully processed and emitted or
/// rejected.
pub(crate) struct TransformArena {
    bump: Bump,
}

impl TransformArena {
    /// Construct an empty reusable transform arena.
    pub fn new() -> Self {
        Self { bump: Bump::new() }
    }

    /// Reset the arena after an execution unit has been fully processed.
    pub fn reset(&mut self) {
        self.bump.reset();
    }

    /// Copy one byte slice into the arena and return a view of the copied bytes.
    pub fn alloc_slice_copy<'a>(&'a self, bytes: &[u8]) -> &'a [u8] {
        self.bump.alloc_slice_copy(bytes)
    }
}

/// Typed rejection reason emitted by a read filter.
pub(crate) trait RejectionReason: std::fmt::Debug + Send + Sync + 'static {
    /// Stable machine-readable reason code suitable for aggregation.
    fn code(&self) -> &'static str;
}

/// Contract for zero-copy read filtering operations.
pub(crate) trait ReadFilter: Send + Sync + 'static {
    /// Typed rejection reason returned when this filter rejects a record.
    type Reason: RejectionReason;

    /// Evaluate the filter on one borrowed record view.
    fn evaluate(&self, record: &RecordView<'_>) -> Result<(), Self::Reason>;
}

/// Contract for read transforms that may rewrite a record using the transform arena.
pub(crate) trait ReadTransform: Send + Sync + 'static {
    /// Stable machine-readable transform code suitable for aggregation.
    fn code(&self) -> &'static str;

    /// Apply the transform and return the resulting record view.
    fn apply<'a>(&self, record: RecordView<'a>, arena: &'a TransformArena) -> TransformResult<'a>;
}

/// Result of applying one transform.
pub(crate) struct TransformResult<'a> {
    /// The resulting record view.
    pub record: RecordView<'a>,
    /// Whether the transform materially changed the record.
    pub applied: bool,
}

/// Outcome of applying one execution step to a record.
pub(crate) enum StepOutcome<'a> {
    /// Continue execution with the returned record view.
    Continue {
        /// Record view after this step.
        record: RecordView<'a>,
        /// Stable transform code when a transform materially changed the record.
        transform_applied: Option<&'static str>,
    },
    /// Reject the record with a boxed rejection reason.
    Reject(Box<dyn RejectionReason>),
}

/// Internal executable step stored by compiled execution plans.
pub(crate) trait ExecutionStep: Send + Sync {
    /// Apply this step to a record.
    fn apply<'a>(&self, record: RecordView<'a>, arena: &'a TransformArena) -> StepOutcome<'a>;
}

/// Internal wrapper that turns a filter operation into a runtime step.
pub(crate) struct FilterStep<F>(pub(crate) F);

impl<F> ExecutionStep for FilterStep<F>
where
    F: ReadFilter,
{
    fn apply<'a>(&self, record: RecordView<'a>, _arena: &'a TransformArena) -> StepOutcome<'a> {
        match self.0.evaluate(&record) {
            Ok(()) => StepOutcome::Continue {
                record,
                transform_applied: None,
            },
            Err(reason) => StepOutcome::Reject(Box::new(reason)),
        }
    }
}

/// Internal wrapper that turns a transform operation into a runtime step.
pub(crate) struct TransformStep<T>(pub(crate) T);

impl<T> ExecutionStep for TransformStep<T>
where
    T: ReadTransform,
{
    fn apply<'a>(&self, record: RecordView<'a>, arena: &'a TransformArena) -> StepOutcome<'a> {
        let result = self.0.apply(record, arena);
        StepOutcome::Continue {
            record: result.record,
            transform_applied: result.applied.then_some(self.0.code()),
        }
    }
}

/// Bridge from authored operations into execution steps.
pub(crate) trait IntoExecutionStep {
    /// Convert the authored operation into a boxed execution step.
    fn into_execution_step(self) -> Box<dyn ExecutionStep>;
}

/// Authored logical plan and compiled execution plan with typestate.
pub(crate) struct Plan<S> {
    steps: Vec<Box<dyn ExecutionStep>>,
    orphan_policy: OrphanPolicy,
    _state: PhantomData<S>,
}

impl Plan<Logical> {
    /// Construct a new empty logical plan.
    pub fn new() -> Self {
        Self {
            steps: Vec::new(),
            orphan_policy: OrphanPolicy::default(),
            _state: PhantomData,
        }
    }
}

/// Trait implemented only by logical plans that can accept new steps and compile.
pub(crate) trait BuildPlan: Sized {
    /// Compiled execution plan type produced when authoring is complete.
    type Execution;

    /// Append one operation to the plan.
    fn step<O>(self, op: O) -> Self
    where
        O: IntoExecutionStep;

    /// Set the orphan policy used by paired execution.
    fn orphan_policy(self, policy: OrphanPolicy) -> Self;

    /// Compile the logical plan into an executable plan.
    fn compile(self) -> Self::Execution;
}

impl BuildPlan for Plan<Logical> {
    type Execution = Plan<Execution>;

    fn step<O>(mut self, op: O) -> Self
    where
        O: IntoExecutionStep,
    {
        self.steps.push(op.into_execution_step());
        self
    }

    fn orphan_policy(mut self, policy: OrphanPolicy) -> Self {
        self.orphan_policy = policy;
        self
    }

    fn compile(self) -> Self::Execution {
        Plan {
            steps: self.steps,
            orphan_policy: self.orphan_policy,
            _state: PhantomData,
        }
    }
}

/// Paired execution unit carrying mate-preserving borrowed record views.
#[derive(Clone, Copy)]
pub(crate) struct RecordPair<'a> {
    /// Left mate.
    pub left: RecordView<'a>,
    /// Right mate.
    pub right: RecordView<'a>,
}

/// Final terminal slot result for one side of an execution unit.
pub(crate) enum ExecutionSlot<'a> {
    /// This side produced a record that should be emitted.
    Emit(RecordView<'a>),
    /// This side was rejected by filtering.
    Reject(Box<dyn RejectionReason>),
    /// This side survived per-read execution but was suppressed by final execution policy.
    Suppress,
}

impl<'a> ExecutionSlot<'a> {
    fn emitted(&self) -> Option<RecordView<'a>> {
        match self {
            Self::Emit(record) => Some(*record),
            Self::Reject(_) | Self::Suppress => None,
        }
    }

    fn is_emit(&self) -> bool {
        matches!(self, Self::Emit(_))
    }

    fn is_reject(&self) -> bool {
        matches!(self, Self::Reject(_))
    }
}

/// Shared final execution result for single-read and paired-read execution.
pub(crate) struct ExecutionOutcome<'a> {
    left: ExecutionSlot<'a>,
    right: Option<ExecutionSlot<'a>>,
}

impl<'a> ExecutionOutcome<'a> {
    /// Construct a single-record outcome.
    pub fn single(slot: ExecutionSlot<'a>) -> Self {
        Self {
            left: slot,
            right: None,
        }
    }

    /// Construct a paired-record outcome.
    pub fn pair(left: ExecutionSlot<'a>, right: ExecutionSlot<'a>) -> Self {
        Self {
            left,
            right: Some(right),
        }
    }

    /// Iterate over records that should be emitted.
    pub fn emitted(&self) -> impl Iterator<Item = RecordView<'a>> + '_ {
        [self.left_emitted(), self.right_emitted()]
            .into_iter()
            .flatten()
    }

    /// Count records that should be emitted.
    pub fn emitted_count(&self) -> usize {
        usize::from(self.left.is_emit())
            + usize::from(self.right.as_ref().is_some_and(ExecutionSlot::is_emit))
    }

    /// Count retained rejection reasons.
    pub fn rejection_count(&self) -> usize {
        usize::from(self.left.is_reject())
            + usize::from(self.right.as_ref().is_some_and(ExecutionSlot::is_reject))
    }

    /// Return the first stable rejection code retained in the final outcome.
    pub fn first_rejection_code(&self) -> Option<&'static str> {
        match (&self.left, &self.right) {
            (ExecutionSlot::Reject(reason), _) | (_, Some(ExecutionSlot::Reject(reason))) => {
                Some(reason.code())
            }
            _ => None,
        }
    }

    /// Return true when execution emitted every slot in the unit.
    pub fn is_fully_emitted(&self) -> bool {
        match &self.right {
            None => self.left.is_emit(),
            Some(right) => self.left.is_emit() && right.is_emit(),
        }
    }

    /// Return true when execution emitted no records.
    pub fn is_fully_rejected(&self) -> bool {
        self.emitted_count() == 0
    }

    /// Return true when paired execution emitted exactly one surviving orphan.
    pub fn is_orphan(&self) -> bool {
        self.right.is_some() && self.emitted_count() == 1
    }

    /// Return the emitted left record when present.
    pub fn left_emitted(&self) -> Option<RecordView<'a>> {
        self.left.emitted()
    }

    /// Return the emitted right record when present.
    pub fn right_emitted(&self) -> Option<RecordView<'a>> {
        self.right.as_ref().and_then(ExecutionSlot::emitted)
    }
}

/// Executable plan surface over a specific execution unit type.
pub(crate) trait Execute<'a, In> {
    /// Execute the compiled plan against one execution unit.
    fn execute(
        &self,
        input: In,
        arena: &'a mut TransformArena,
        stats: &mut ReadStats,
    ) -> ExecutionOutcome<'a>;
}

impl Plan<Execution> {
    /// Execute all compiled steps against one record and return the terminal slot result.
    fn execute_record<'a>(
        &self,
        record: RecordView<'a>,
        arena: &'a TransformArena,
        stats: &mut ReadStats,
    ) -> ExecutionSlot<'a> {
        let mut current = record;

        for step in &self.steps {
            match step.apply(current, arena) {
                StepOutcome::Continue {
                    record,
                    transform_applied,
                } => {
                    current = record;
                    if let Some(code) = transform_applied {
                        stats.record_transform(code);
                    }
                }
                StepOutcome::Reject(reason) => {
                    stats.record_rejected(reason.code());
                    return ExecutionSlot::Reject(reason);
                }
            }
        }

        ExecutionSlot::Emit(current)
    }
}

impl<'a> Execute<'a, RecordView<'a>> for Plan<Execution> {
    fn execute(
        &self,
        input: RecordView<'a>,
        arena: &'a mut TransformArena,
        stats: &mut ReadStats,
    ) -> ExecutionOutcome<'a> {
        ExecutionOutcome::single(self.execute_record(input, arena, stats))
    }
}

impl<'a> Execute<'a, RecordPair<'a>> for Plan<Execution> {
    fn execute(
        &self,
        input: RecordPair<'a>,
        arena: &'a mut TransformArena,
        stats: &mut ReadStats,
    ) -> ExecutionOutcome<'a> {
        let left = self.execute_record(input.left, arena, stats);
        let right = self.execute_record(input.right, arena, stats);

        match self.orphan_policy {
            OrphanPolicy::DropPair => match (left, right) {
                (ExecutionSlot::Emit(left_record), ExecutionSlot::Emit(right_record)) => {
                    ExecutionOutcome::pair(
                        ExecutionSlot::Emit(left_record),
                        ExecutionSlot::Emit(right_record),
                    )
                }
                (ExecutionSlot::Emit(_), right) => {
                    ExecutionOutcome::pair(ExecutionSlot::Suppress, right)
                }
                (left, ExecutionSlot::Emit(_)) => {
                    ExecutionOutcome::pair(left, ExecutionSlot::Suppress)
                }
                (left, right) => ExecutionOutcome::pair(left, right),
            },
            OrphanPolicy::EmitOrphan => ExecutionOutcome::pair(left, right),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{
        BuildPlan, Execute, ExecutionOutcome, FilterStep, IntoExecutionStep, Logical, OrphanPolicy,
        Plan, ReadFilter, ReadTransform, RecordPair, RejectionReason, TransformArena,
        TransformResult, TransformStep,
    };
    use crate::record::{ReadStats, RecordView};

    #[derive(Debug)]
    struct TooShort;

    impl RejectionReason for TooShort {
        fn code(&self) -> &'static str {
            "too_short"
        }
    }

    struct MinLength {
        min_length: usize,
    }

    impl MinLength {
        fn new(min_length: usize) -> Self {
            Self { min_length }
        }
    }

    impl ReadFilter for MinLength {
        type Reason = TooShort;

        fn evaluate(&self, record: &RecordView<'_>) -> Result<(), Self::Reason> {
            if record.sequence().len() < self.min_length {
                Err(TooShort)
            } else {
                Ok(())
            }
        }
    }

    impl IntoExecutionStep for MinLength {
        fn into_execution_step(self) -> Box<dyn super::ExecutionStep> {
            Box::new(FilterStep(self))
        }
    }

    struct TrimPrefix {
        amount: usize,
    }

    impl TrimPrefix {
        fn new(amount: usize) -> Self {
            Self { amount }
        }
    }

    impl ReadTransform for TrimPrefix {
        fn code(&self) -> &'static str {
            "trim_prefix"
        }

        fn apply<'a>(
            &self,
            record: RecordView<'a>,
            arena: &'a TransformArena,
        ) -> TransformResult<'a> {
            TransformResult {
                record: record
                    .with_sequence_and_quality(
                        arena.alloc_slice_copy(&record.sequence()[self.amount..]),
                        arena.alloc_slice_copy(&record.quality()[self.amount..]),
                    )
                    .expect("trim prefix should preserve equal sequence and quality lengths"),
                applied: true,
            }
        }
    }

    impl IntoExecutionStep for TrimPrefix {
        fn into_execution_step(self) -> Box<dyn super::ExecutionStep> {
            Box::new(TransformStep(self))
        }
    }

    fn record(sequence: &'static [u8]) -> RecordView<'static> {
        let quality = match sequence.len() {
            6 => b"IIIIII".as_slice(),
            _ => b"IIII".as_slice(),
        };
        RecordView::new(b"read1", sequence, quality)
    }

    #[test]
    fn logical_plan_accumulates_steps_and_compiles() {
        let _plan = Plan::<Logical>::new()
            .step(MinLength::new(4))
            .step(TrimPrefix::new(1))
            .compile();
    }

    #[test]
    fn single_execution_rejects_record_when_filter_fails() {
        let plan = Plan::<Logical>::new().step(MinLength::new(6)).compile();
        let mut arena = TransformArena::new();
        let mut stats = ReadStats::default();

        let outcome = plan.execute(record(b"ACGT"), &mut arena, &mut stats);

        assert!(outcome.is_fully_rejected());
        assert_eq!(outcome.rejection_count(), 1);
        assert_eq!(stats.rejection_counts.get("too_short"), Some(&1));
    }

    #[test]
    fn single_execution_applies_transform_in_order() {
        let plan = Plan::<Logical>::new().step(TrimPrefix::new(1)).compile();
        let mut arena = TransformArena::new();
        let mut stats = ReadStats::default();

        let outcome = plan.execute(record(b"ACGT"), &mut arena, &mut stats);

        assert!(outcome.is_fully_emitted());
        assert_eq!(
            outcome
                .left_emitted()
                .expect("record should emit")
                .sequence(),
            b"CGT"
        );
    }

    #[test]
    fn paired_execution_drops_orphan_by_default() {
        let plan = Plan::<Logical>::new().step(MinLength::new(6)).compile();
        let mut arena = TransformArena::new();
        let mut stats = ReadStats::default();

        let outcome = plan.execute(
            RecordPair {
                left: record(b"ACGTAC"),
                right: record(b"ACGT"),
            },
            &mut arena,
            &mut stats,
        );

        assert!(outcome.is_fully_rejected());
        assert!(!outcome.is_orphan());
        assert_eq!(outcome.emitted_count(), 0);
        assert_eq!(outcome.rejection_count(), 1);
    }

    #[test]
    fn paired_execution_can_emit_orphan_when_policy_allows() {
        let plan = Plan::<Logical>::new()
            .step(MinLength::new(6))
            .orphan_policy(OrphanPolicy::EmitOrphan)
            .compile();
        let mut arena = TransformArena::new();
        let mut stats = ReadStats::default();

        let outcome = plan.execute(
            RecordPair {
                left: record(b"ACGTAC"),
                right: record(b"ACGT"),
            },
            &mut arena,
            &mut stats,
        );

        assert!(outcome.is_orphan());
        assert_eq!(outcome.emitted_count(), 1);
        assert_eq!(outcome.rejection_count(), 1);
    }

    #[test]
    fn emitted_iterator_preserves_slot_order() {
        let outcome = ExecutionOutcome::pair(
            super::ExecutionSlot::Emit(record(b"AAAAAA")),
            super::ExecutionSlot::Emit(record(b"CCCCCC")),
        );

        let emitted = outcome
            .emitted()
            .map(|record| record.sequence())
            .collect::<Vec<_>>();
        assert_eq!(emitted, vec![b"AAAAAA".as_slice(), b"CCCCCC".as_slice()]);
    }
}
