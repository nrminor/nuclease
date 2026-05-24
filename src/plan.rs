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
    #[allow(
        dead_code,
        reason = "orphan emission is supported by the plan model even though the current CLI keeps paired output conservative"
    )]
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

    /// Copy one byte slice into the arena and return a mutable view of the copied bytes.
    pub fn alloc_slice_copy_mut<'a>(&'a self, bytes: &[u8]) -> &'a mut [u8] {
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

/// Contract for read transforms that may rewrite one record using the transform arena.
pub(crate) trait ReadTransform: Send + Sync + 'static {
    /// Stable machine-readable transform code suitable for aggregation.
    fn code(&self) -> &'static str;

    /// Apply the transform and return the resulting record view.
    fn apply<'a>(&self, record: RecordView<'a>, arena: &'a TransformArena) -> TransformResult<'a>;
}

/// Active execution unit carried between plan steps.
#[derive(Clone, Copy)]
pub(crate) enum ActiveUnit<'a> {
    /// One single-end record, or one paired record that has become a single record.
    Single(RecordView<'a>),
    /// One paired-end execution unit.
    Pair(RecordPair<'a>),
}

/// Emitted unit produced by an execution plan.
#[derive(Clone, Copy)]
pub(crate) enum EmittedUnit<'a> {
    /// Emit no records.
    None,
    /// Emit one record.
    Single(RecordView<'a>),
    /// Emit a paired record group.
    Pair(RecordPair<'a>),
}

/// Result of a pair-aware transform.
pub(crate) enum PairTransformResult<'a> {
    /// Continue as a pair.
    Pair {
        /// Pair after transformation.
        pair: RecordPair<'a>,
        /// Whether the transform materially changed the unit.
        applied: bool,
    },
    /// Continue as one record.
    Single {
        /// Record after transformation.
        record: RecordView<'a>,
        /// Whether the transform materially changed the unit.
        applied: bool,
    },
    /// Drop the entire unit with a stable rejection reason.
    #[allow(
        dead_code,
        reason = "pair-aware extension points may reject whole units even though merge-pairs currently keeps unmerged pairs"
    )]
    Drop {
        /// Stable rejection reason code.
        reason: &'static str,
    },
}

/// Contract for transforms that operate on a paired unit and may change output cardinality.
pub(crate) trait PairTransform: Send {
    /// Stable machine-readable transform code suitable for aggregation.
    fn code(&self) -> &'static str;

    /// Apply this transform to a paired unit.
    fn apply_pair<'a>(
        &mut self,
        pair: RecordPair<'a>,
        arena: &'a TransformArena,
    ) -> Result<PairTransformResult<'a>>;

    /// Apply this transform to a single unit.
    fn apply_single<'a>(
        &mut self,
        record: RecordView<'a>,
        _arena: &'a TransformArena,
    ) -> Result<PairTransformResult<'a>> {
        Ok(PairTransformResult::Single {
            record,
            applied: false,
        })
    }
}

/// Result of applying one transform.
pub(crate) struct TransformResult<'a> {
    /// The resulting record view.
    pub record: RecordView<'a>,
    /// Whether the transform materially changed the record.
    pub applied: bool,
}

/// Outcome of applying one execution step to an active unit.
pub(crate) enum StepOutcome<'a> {
    /// Continue execution with the returned unit.
    Continue(ActiveUnit<'a>),
    /// Stop execution and return the terminal outcome.
    Stop(ExecutionOutcome<'a>),
}

/// Mutable execution context shared by runtime steps for one execution unit.
pub(crate) struct StepContext<'a, 'stats> {
    arena: &'a TransformArena,
    stats: &'stats mut ReadStats,
    orphan_policy: OrphanPolicy,
    rejection_count: usize,
}

impl StepContext<'_, '_> {
    fn record_rejection(&mut self, code: &'static str) {
        self.stats.record_rejected(code);
        self.rejection_count += 1;
    }

    fn record_transform(&mut self, code: &'static str) {
        self.stats.record_transform(code);
    }

    fn outcome<'a>(&self, emitted: EmittedUnit<'a>) -> ExecutionOutcome<'a> {
        ExecutionOutcome {
            emitted,
            rejection_count: self.rejection_count,
        }
    }
}

/// Internal executable step stored by compiled execution plans.
pub(crate) trait ExecutionStep: Send {
    /// Apply this step to an active execution unit.
    fn apply<'a>(
        &mut self,
        unit: ActiveUnit<'a>,
        context: &mut StepContext<'a, '_>,
    ) -> Result<StepOutcome<'a>>;
}

/// Internal wrapper that turns a filter operation into a runtime step.
pub(crate) struct FilterStep<F>(pub(crate) F);

impl<F> ExecutionStep for FilterStep<F>
where
    F: ReadFilter,
{
    fn apply<'a>(
        &mut self,
        unit: ActiveUnit<'a>,
        context: &mut StepContext<'a, '_>,
    ) -> Result<StepOutcome<'a>> {
        Ok(match unit {
            ActiveUnit::Single(record) => match self.0.evaluate(&record) {
                Ok(()) => StepOutcome::Continue(ActiveUnit::Single(record)),
                Err(reason) => {
                    context.record_rejection(reason.code());
                    StepOutcome::Stop(context.outcome(EmittedUnit::None))
                }
            },
            ActiveUnit::Pair(pair) => apply_filter_to_pair(&self.0, pair, context),
        })
    }
}

fn apply_filter_to_pair<'a, F>(
    filter: &F,
    pair: RecordPair<'a>,
    context: &mut StepContext<'a, '_>,
) -> StepOutcome<'a>
where
    F: ReadFilter,
{
    let left = filter.evaluate(&pair.left).map_err(|reason| reason.code());
    let right = filter.evaluate(&pair.right).map_err(|reason| reason.code());

    match (left, right) {
        (Ok(()), Ok(())) => StepOutcome::Continue(ActiveUnit::Pair(pair)),
        (Err(left_reason), Ok(())) => {
            context.record_rejection(left_reason);
            match context.orphan_policy {
                OrphanPolicy::DropPair => StepOutcome::Stop(context.outcome(EmittedUnit::None)),
                OrphanPolicy::EmitOrphan => StepOutcome::Continue(ActiveUnit::Single(pair.right)),
            }
        }
        (Ok(()), Err(right_reason)) => {
            context.record_rejection(right_reason);
            match context.orphan_policy {
                OrphanPolicy::DropPair => StepOutcome::Stop(context.outcome(EmittedUnit::None)),
                OrphanPolicy::EmitOrphan => StepOutcome::Continue(ActiveUnit::Single(pair.left)),
            }
        }
        (Err(left_reason), Err(right_reason)) => {
            context.record_rejection(left_reason);
            context.record_rejection(right_reason);
            StepOutcome::Stop(context.outcome(EmittedUnit::None))
        }
    }
}

/// Internal wrapper that turns a transform operation into a runtime step.
pub(crate) struct TransformStep<T>(pub(crate) T);

impl<T> ExecutionStep for TransformStep<T>
where
    T: ReadTransform,
{
    fn apply<'a>(
        &mut self,
        unit: ActiveUnit<'a>,
        context: &mut StepContext<'a, '_>,
    ) -> Result<StepOutcome<'a>> {
        let code = self.0.code();
        let mut map_one = |record| {
            let result = self.0.apply(record, context.arena);
            if result.applied {
                context.record_transform(code);
            }
            result.record
        };

        Ok(StepOutcome::Continue(match unit {
            ActiveUnit::Single(record) => ActiveUnit::Single(map_one(record)),
            ActiveUnit::Pair(pair) => ActiveUnit::Pair(RecordPair {
                left: map_one(pair.left),
                right: map_one(pair.right),
            }),
        }))
    }
}

/// Internal wrapper that turns a pair-aware transform into a runtime step.
pub(crate) struct PairTransformStep<T>(pub(crate) T);

impl<T> ExecutionStep for PairTransformStep<T>
where
    T: PairTransform,
{
    fn apply<'a>(
        &mut self,
        unit: ActiveUnit<'a>,
        context: &mut StepContext<'a, '_>,
    ) -> Result<StepOutcome<'a>> {
        let result = match unit {
            ActiveUnit::Single(record) => self.0.apply_single(record, context.arena)?,
            ActiveUnit::Pair(pair) => self.0.apply_pair(pair, context.arena)?,
        };

        Ok(match result {
            PairTransformResult::Pair { pair, applied } => {
                if applied {
                    context.record_transform(self.0.code());
                }
                StepOutcome::Continue(ActiveUnit::Pair(pair))
            }
            PairTransformResult::Single { record, applied } => {
                if applied {
                    context.record_transform(self.0.code());
                }
                StepOutcome::Continue(ActiveUnit::Single(record))
            }
            PairTransformResult::Drop { reason } => {
                context.record_rejection(reason);
                StepOutcome::Stop(context.outcome(EmittedUnit::None))
            }
        })
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

/// Shared final execution result for single-read and paired-read execution.
pub(crate) struct ExecutionOutcome<'a> {
    emitted: EmittedUnit<'a>,
    rejection_count: usize,
}

impl<'a> ExecutionOutcome<'a> {
    /// Return the emitted unit.
    pub fn emitted_unit(&self) -> EmittedUnit<'a> {
        self.emitted
    }

    /// Iterate over records that should be emitted.
    pub fn emitted(&self) -> impl Iterator<Item = RecordView<'a>> {
        let records = match self.emitted {
            EmittedUnit::None => [None, None],
            EmittedUnit::Single(record) => [Some(record), None],
            EmittedUnit::Pair(pair) => [Some(pair.left), Some(pair.right)],
        };
        records.into_iter().flatten()
    }

    /// Count retained rejection reasons.
    pub fn rejection_count(&self) -> usize {
        self.rejection_count
    }
}

/// Executable plan surface over a specific execution unit type.
pub(crate) trait Execute<'a, In> {
    /// Execute the compiled plan against one execution unit.
    fn execute(
        &mut self,
        input: In,
        arena: &'a mut TransformArena,
        stats: &mut ReadStats,
    ) -> Result<ExecutionOutcome<'a>>;
}

impl Plan<Execution> {
    /// Execute all compiled steps against one active unit.
    fn execute_unit<'a>(
        &mut self,
        mut unit: ActiveUnit<'a>,
        arena: &'a mut TransformArena,
        stats: &mut ReadStats,
    ) -> Result<ExecutionOutcome<'a>> {
        let mut context = StepContext {
            arena,
            stats,
            orphan_policy: self.orphan_policy,
            rejection_count: 0,
        };

        for step in &mut self.steps {
            match step.apply(unit, &mut context)? {
                StepOutcome::Continue(next) => unit = next,
                StepOutcome::Stop(outcome) => return Ok(outcome),
            }
        }

        Ok(context.outcome(match unit {
            ActiveUnit::Single(record) => EmittedUnit::Single(record),
            ActiveUnit::Pair(pair) => EmittedUnit::Pair(pair),
        }))
    }
}

impl<'a> Execute<'a, RecordView<'a>> for Plan<Execution> {
    fn execute(
        &mut self,
        input: RecordView<'a>,
        arena: &'a mut TransformArena,
        stats: &mut ReadStats,
    ) -> Result<ExecutionOutcome<'a>> {
        self.execute_unit(ActiveUnit::Single(input), arena, stats)
    }
}

impl<'a> Execute<'a, RecordPair<'a>> for Plan<Execution> {
    fn execute(
        &mut self,
        input: RecordPair<'a>,
        arena: &'a mut TransformArena,
        stats: &mut ReadStats,
    ) -> Result<ExecutionOutcome<'a>> {
        self.execute_unit(ActiveUnit::Pair(input), arena, stats)
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
            _arena: &'a TransformArena,
        ) -> TransformResult<'a> {
            TransformResult {
                record: record
                    .with_sequence_and_quality(
                        &record.sequence()[self.amount..],
                        &record.quality()[self.amount..],
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
        let mut plan = Plan::<Logical>::new().step(MinLength::new(6)).compile();
        let mut arena = TransformArena::new();
        let mut stats = ReadStats::default();

        let outcome = plan
            .execute(record(b"ACGT"), &mut arena, &mut stats)
            .expect("single-record filter failure should produce a rejected outcome");

        assert_eq!(outcome.emitted().count(), 0);
        assert_eq!(outcome.rejection_count(), 1);
        assert_eq!(stats.rejection_counts.get("too_short"), Some(&1));
    }

    #[test]
    fn single_execution_applies_transform_in_order() {
        let mut plan = Plan::<Logical>::new().step(TrimPrefix::new(1)).compile();
        let mut arena = TransformArena::new();
        let mut stats = ReadStats::default();

        let outcome = plan
            .execute(record(b"ACGT"), &mut arena, &mut stats)
            .expect("single-record transform should produce an emitted outcome");

        assert_eq!(outcome.emitted().count(), 1);
        assert_eq!(
            outcome
                .emitted()
                .next()
                .expect("record should emit")
                .sequence(),
            b"CGT"
        );
    }

    #[test]
    fn paired_execution_drops_orphan_by_default() {
        let mut plan = Plan::<Logical>::new().step(MinLength::new(6)).compile();
        let mut arena = TransformArena::new();
        let mut stats = ReadStats::default();

        let outcome = plan
            .execute(
                RecordPair {
                    left: record(b"ACGTAC"),
                    right: record(b"ACGT"),
                },
                &mut arena,
                &mut stats,
            )
            .expect("paired filter failure should produce a rejected outcome");

        assert_eq!(outcome.emitted().count(), 0);
        assert_eq!(outcome.rejection_count(), 1);
    }

    #[test]
    fn paired_execution_can_emit_orphan_when_policy_allows() {
        let mut plan = Plan::<Logical>::new()
            .step(MinLength::new(6))
            .orphan_policy(OrphanPolicy::EmitOrphan)
            .compile();
        let mut arena = TransformArena::new();
        let mut stats = ReadStats::default();

        let outcome = plan
            .execute(
                RecordPair {
                    left: record(b"ACGTAC"),
                    right: record(b"ACGT"),
                },
                &mut arena,
                &mut stats,
            )
            .expect("paired orphan policy should produce an emitted orphan outcome");

        assert_eq!(outcome.emitted().count(), 1);
        assert_eq!(outcome.rejection_count(), 1);
    }

    #[test]
    fn emitted_iterator_preserves_slot_order() {
        let outcome = ExecutionOutcome {
            emitted: super::EmittedUnit::Pair(RecordPair {
                left: record(b"AAAAAA"),
                right: record(b"CCCCCC"),
            }),
            rejection_count: 0,
        };

        let emitted = outcome
            .emitted()
            .map(|record| record.sequence())
            .collect::<Vec<_>>();
        assert_eq!(emitted, vec![b"AAAAAA".as_slice(), b"CCCCCC".as_slice()]);
    }
}
