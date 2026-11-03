//! Logical authoring and execution substrate for streaming preprocessing plans.

use std::marker::PhantomData;

use bumpalo::Bump;

use crate::{
    error::{InternalError, Result},
    observer::RunObserver,
    record::{RecordPair, RecordView},
};

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
pub(crate) trait ReadFilter: Send + 'static {
    /// Typed rejection reason returned when this filter rejects a record.
    type Reason: RejectionReason;

    /// Evaluate the filter on one borrowed record view.
    fn evaluate(&mut self, record: &RecordView<'_>) -> std::result::Result<(), Self::Reason>;

    /// Evaluate a pair, independently by default.
    fn evaluate_pair(&mut self, pair: &RecordPair<'_>) -> PairRelation<Self::Reason> {
        PairRelation::PerMate {
            left: self.evaluate(&pair.left()),
            right: self.evaluate(&pair.right()),
        }
    }
}

/// A pair decision, either per mate or atomic for the whole pair.
pub(crate) enum PairRelation<R> {
    /// Independent decisions for each mate.
    PerMate {
        left: std::result::Result<(), R>,
        right: std::result::Result<(), R>,
    },
    /// One decision covering both mates.
    #[allow(dead_code, reason = "whole-pair filters are an extension point")]
    WholePair(std::result::Result<(), R>),
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

impl<'record> From<RecordPair<'record>> for ActiveUnit<'record> {
    fn from(pair: RecordPair<'record>) -> Self {
        Self::Pair(pair)
    }
}

impl<'record> ActiveUnit<'record> {
    fn map_records(self, mut map: impl FnMut(RecordView<'record>) -> RecordView<'record>) -> Self {
        match self {
            Self::Single(record) => Self::Single(map(record)),
            Self::Pair(pair) => {
                let left = map(pair.left());
                let right = map(pair.right());
                Self::Pair(pair.with_records(left, right))
            }
        }
    }
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
pub(crate) struct StepContext<'a, 'observer> {
    arena: &'a TransformArena,
    observer: &'observer mut RunObserver,
    orphan_policy: OrphanPolicy,
    rejection_count: usize,
}

impl StepContext<'_, '_> {
    fn record_rejection(&mut self, code: &'static str) {
        self.record_rejection_for(code, 1);
    }

    fn record_rejection_for(&mut self, code: &'static str, records: usize) {
        for _ in 0..records {
            self.observer.record_rejected(code);
        }
        self.rejection_count += 1;
    }

    fn record_transform(&mut self, code: &'static str) {
        self.observer.record_transform(code);
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

/// Source of active read units for whole-read-set transforms.
#[allow(
    dead_code,
    reason = "read-set execution is scaffolded before concrete read-set transforms"
)]
pub(crate) trait ActiveUnitSource {
    /// Return the next active unit, or `None` when the source is exhausted.
    fn next_unit(&mut self) -> Result<Option<ActiveUnit<'_>>>;
}

/// Sink for active read units emitted by whole-read-set transforms.
#[allow(
    dead_code,
    reason = "read-set execution is scaffolded before concrete read-set transforms"
)]
pub(crate) trait ActiveUnitSink {
    /// Emit one active unit downstream.
    fn emit_unit(&mut self, unit: ActiveUnit<'_>) -> Result<()>;
}

/// Contract for transforms that operate over a streamable read set.
#[allow(
    dead_code,
    reason = "read-set execution is scaffolded before concrete read-set transforms"
)]
pub(crate) trait ReadSetTransform: Send + 'static {
    /// Stable machine-readable transform code suitable for aggregation.
    fn code(&self) -> &'static str;

    /// Consume active units from `input` and emit transformed units to `output`.
    ///
    /// Units borrowed from `input` are short-lived. Implementations that need to retain or reorder
    /// records must copy, spill, or otherwise store them using operation-owned storage.
    fn transform_read_set(
        &mut self,
        input: &mut dyn ActiveUnitSource,
        output: &mut dyn ActiveUnitSink,
    ) -> Result<()>;
}

enum LogicalPlanItem {
    PerUnit(Box<dyn PrepareExecutionStep>),
    #[allow(
        dead_code,
        reason = "read-set logical items are currently exercised by tests only"
    )]
    ReadSet(Box<dyn PrepareExecutionStep<dyn ReadSetTransform>>),
}

enum PhysicalPlanItem {
    PerUnitChain(PerUnitChain),
    ReadSet(Box<dyn ReadSetTransform>),
}

struct PerUnitChain {
    steps: Vec<Box<dyn ExecutionStep>>,
}

impl PerUnitChain {
    fn new(steps: Vec<Box<dyn ExecutionStep>>) -> Self {
        Self { steps }
    }

    fn push(&mut self, step: Box<dyn ExecutionStep>) {
        self.steps.push(step);
    }

    fn execute_unit<'a>(
        &mut self,
        mut unit: ActiveUnit<'a>,
        context: &mut StepContext<'a, '_>,
    ) -> Result<StepOutcome<'a>> {
        for step in &mut self.steps {
            match step.apply(unit, context)? {
                StepOutcome::Continue(next) => unit = next,
                StepOutcome::Stop(outcome) => return Ok(StepOutcome::Stop(outcome)),
            }
        }

        Ok(StepOutcome::Continue(unit))
    }
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
            ActiveUnit::Pair(pair) => apply_filter_to_pair(&mut self.0, pair, context),
        })
    }
}

fn apply_filter_to_pair<'a, F>(
    filter: &mut F,
    pair: RecordPair<'a>,
    context: &mut StepContext<'a, '_>,
) -> StepOutcome<'a>
where
    F: ReadFilter,
{
    let (left, right) = match filter.evaluate_pair(&pair) {
        PairRelation::WholePair(Ok(())) => return StepOutcome::Continue(ActiveUnit::Pair(pair)),
        PairRelation::WholePair(Err(reason)) => {
            context.record_rejection_for(reason.code(), 2);
            return StepOutcome::Stop(context.outcome(EmittedUnit::None));
        }
        PairRelation::PerMate { left, right } => {
            (left.map_err(|r| r.code()), right.map_err(|r| r.code()))
        }
    };

    match (left, right) {
        (Ok(()), Ok(())) => StepOutcome::Continue(ActiveUnit::Pair(pair)),
        (Err(left_reason), Ok(())) => {
            context.record_rejection(left_reason);
            match context.orphan_policy {
                OrphanPolicy::DropPair => {
                    context.observer.reads_suppressed += 1;
                    StepOutcome::Stop(context.outcome(EmittedUnit::None))
                }
                OrphanPolicy::EmitOrphan => StepOutcome::Continue(ActiveUnit::Single(pair.right())),
            }
        }
        (Ok(()), Err(right_reason)) => {
            context.record_rejection(right_reason);
            match context.orphan_policy {
                OrphanPolicy::DropPair => {
                    context.observer.reads_suppressed += 1;
                    StepOutcome::Stop(context.outcome(EmittedUnit::None))
                }
                OrphanPolicy::EmitOrphan => StepOutcome::Continue(ActiveUnit::Single(pair.left())),
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

        Ok(StepOutcome::Continue(unit.map_records(&mut map_one)))
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
        let cardinality = match unit {
            ActiveUnit::Single(_) => 1,
            ActiveUnit::Pair(_) => 2,
        };
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
                context.record_rejection_for(reason, cardinality);
                StepOutcome::Stop(context.outcome(EmittedUnit::None))
            }
        })
    }
}

/// Bridge from reusable operation definitions to owned execution state.
///
/// Each preparation creates fresh mutable state. Immutable resources may be shared, but the
/// returned worker must not borrow this definition. Preparation must leave the definition reusable
/// on both success and failure. Workers need only be `Send`, not `Sync`.
pub(crate) trait PrepareExecutionStep<T: ?Sized = dyn ExecutionStep>:
    Send + Sync + 'static
{
    /// Prepare one independent execution instance, releasing partial resources on failure.
    fn prepare(&self) -> Result<Box<T>>;
}

/// Authored logical plan and compiled execution plan with typestate.
pub(crate) struct Plan<S> {
    logical_items: Vec<LogicalPlanItem>,
    physical_items: Vec<PhysicalPlanItem>,
    orphan_policy: OrphanPolicy,
    _state: PhantomData<S>,
}

impl Plan<Logical> {
    /// Construct a new empty logical plan.
    pub fn new() -> Self {
        Self {
            logical_items: Vec::new(),
            physical_items: Vec::new(),
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
        O: PrepareExecutionStep;

    /// Append one whole-read-set transform to the plan.
    #[allow(
        dead_code,
        reason = "read-set transforms are currently scaffolded for future features"
    )]
    fn read_set_transform<T>(self, transform: T) -> Self
    where
        T: PrepareExecutionStep<dyn ReadSetTransform>;

    /// Set the orphan policy used by paired execution.
    fn orphan_policy(self, policy: OrphanPolicy) -> Self;

    /// Prepare an owned executable plan without consuming the logical plan.
    fn compile(&self) -> Result<Self::Execution>;
}

impl BuildPlan for Plan<Logical> {
    type Execution = Plan<Execution>;

    fn step<O>(mut self, op: O) -> Self
    where
        O: PrepareExecutionStep,
    {
        self.logical_items
            .push(LogicalPlanItem::PerUnit(Box::new(op)));
        self
    }

    fn read_set_transform<T>(mut self, transform: T) -> Self
    where
        T: PrepareExecutionStep<dyn ReadSetTransform>,
    {
        self.logical_items
            .push(LogicalPlanItem::ReadSet(Box::new(transform)));
        self
    }

    fn orphan_policy(mut self, policy: OrphanPolicy) -> Self {
        self.orphan_policy = policy;
        self
    }

    fn compile(&self) -> Result<Self::Execution> {
        let stage_count = self
            .logical_items
            .chunk_by(|left, right| {
                matches!(left, LogicalPlanItem::PerUnit(_))
                    && matches!(right, LogicalPlanItem::PerUnit(_))
            })
            .count();
        let mut physical_items = Vec::with_capacity(stage_count);

        for (index, item) in self.logical_items.iter().enumerate() {
            match item {
                LogicalPlanItem::PerUnit(definition) => {
                    let step = definition.prepare()?;
                    if let Some(PhysicalPlanItem::PerUnitChain(chain)) = physical_items.last_mut() {
                        chain.push(step);
                    } else {
                        let chain_len = self.logical_items[index..]
                            .iter()
                            .take_while(|item| matches!(item, LogicalPlanItem::PerUnit(_)))
                            .count();
                        let mut steps = Vec::with_capacity(chain_len);
                        steps.push(step);
                        physical_items
                            .push(PhysicalPlanItem::PerUnitChain(PerUnitChain::new(steps)));
                    }
                }
                LogicalPlanItem::ReadSet(transform) => {
                    physical_items.push(PhysicalPlanItem::ReadSet(transform.prepare()?));
                }
            }
        }

        Ok(Plan {
            logical_items: Vec::new(),
            physical_items,
            orphan_policy: self.orphan_policy,
            _state: PhantomData,
        })
    }
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
            EmittedUnit::Pair(pair) => [Some(pair.left()), Some(pair.right())],
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
        observer: &mut RunObserver,
    ) -> Result<ExecutionOutcome<'a>>;
}

impl Plan<Execution> {
    /// Execute all compiled steps against one active unit.
    fn execute_unit<'a>(
        &mut self,
        mut unit: ActiveUnit<'a>,
        arena: &'a mut TransformArena,
        observer: &mut RunObserver,
    ) -> Result<ExecutionOutcome<'a>> {
        let mut context = StepContext {
            arena,
            observer,
            orphan_policy: self.orphan_policy,
            rejection_count: 0,
        };

        for item in &mut self.physical_items {
            match item {
                PhysicalPlanItem::PerUnitChain(chain) => {
                    match chain.execute_unit(unit, &mut context)? {
                        StepOutcome::Continue(next) => unit = next,
                        StepOutcome::Stop(outcome) => return Ok(outcome),
                    }
                }
                PhysicalPlanItem::ReadSet(transform) => {
                    return Err(InternalError::PlanInvariant {
                        detail: format!(
                            "read-set transform `{}` requires whole-read-set execution",
                            transform.code()
                        ),
                    }
                    .into());
                }
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
        observer: &mut RunObserver,
    ) -> Result<ExecutionOutcome<'a>> {
        self.execute_unit(ActiveUnit::Single(input), arena, observer)
    }
}

impl<'a> Execute<'a, RecordPair<'a>> for Plan<Execution> {
    fn execute(
        &mut self,
        input: RecordPair<'a>,
        arena: &'a mut TransformArena,
        observer: &mut RunObserver,
    ) -> Result<ExecutionOutcome<'a>> {
        self.execute_unit(input.into(), arena, observer)
    }
}

#[cfg(test)]
mod tests {
    use super::{
        ActiveUnitSink, ActiveUnitSource, BuildPlan, Execute, ExecutionOutcome, FilterStep,
        Logical, OrphanPolicy, PairRelation, PhysicalPlanItem, Plan, PrepareExecutionStep,
        ReadFilter, ReadSetTransform, ReadTransform, RecordPair, RejectionReason, TransformArena,
        TransformResult, TransformStep,
    };
    use crate::{
        error::{InternalError, RunError},
        observer::RunObserver,
        record::RecordView,
    };
    use color_eyre::eyre::Result;
    use std::{
        cell::Cell,
        sync::{
            Arc,
            atomic::{AtomicUsize, Ordering},
        },
    };

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

        fn evaluate(&mut self, record: &RecordView<'_>) -> Result<(), Self::Reason> {
            if record.sequence().len() < self.min_length {
                Err(TooShort)
            } else {
                Ok(())
            }
        }
    }

    impl PrepareExecutionStep for MinLength {
        fn prepare(&self) -> crate::error::Result<Box<dyn super::ExecutionStep>> {
            Ok(Box::new(FilterStep(Self::new(self.min_length))))
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

    impl PrepareExecutionStep for TrimPrefix {
        fn prepare(&self) -> crate::error::Result<Box<dyn super::ExecutionStep>> {
            Ok(Box::new(TransformStep(Self::new(self.amount))))
        }
    }

    struct FakeReadSetTransform;

    impl PrepareExecutionStep<dyn ReadSetTransform> for FakeReadSetTransform {
        fn prepare(&self) -> crate::error::Result<Box<dyn ReadSetTransform>> {
            Ok(Box::new(Self))
        }
    }

    impl ReadSetTransform for FakeReadSetTransform {
        fn code(&self) -> &'static str {
            "fake_read_set"
        }

        fn transform_read_set(
            &mut self,
            _input: &mut dyn ActiveUnitSource,
            _output: &mut dyn ActiveUnitSink,
        ) -> crate::error::Result<()> {
            Ok(())
        }
    }

    fn record(sequence: &'static [u8]) -> RecordView<'static> {
        let quality = match sequence.len() {
            6 => b"IIIIII".as_slice(),
            _ => b"IIII".as_slice(),
        };
        RecordView::new(b"read1", sequence, quality)
    }

    fn pair(left: &'static [u8], right: &'static [u8]) -> RecordPair<'static> {
        RecordPair::try_new(record(left).into(), record(right).into())
            .expect("test records should form a valid pair")
    }

    struct TrackedWorker(Arc<AtomicUsize>);
    impl Drop for TrackedWorker {
        fn drop(&mut self) {
            self.0.fetch_add(1, Ordering::SeqCst);
        }
    }
    impl super::ExecutionStep for TrackedWorker {
        fn apply<'a>(
            &mut self,
            unit: super::ActiveUnit<'a>,
            _context: &mut super::StepContext<'a, '_>,
        ) -> crate::error::Result<super::StepOutcome<'a>> {
            Ok(super::StepOutcome::Continue(unit))
        }
    }
    impl ReadSetTransform for TrackedWorker {
        fn code(&self) -> &'static str {
            "tracked"
        }
        fn transform_read_set(
            &mut self,
            _input: &mut dyn ActiveUnitSource,
            _output: &mut dyn ActiveUnitSink,
        ) -> crate::error::Result<()> {
            Ok(())
        }
    }
    struct TrackedDefinition {
        drops: Arc<AtomicUsize>,
        fail: bool,
    }
    impl TrackedDefinition {
        fn worker(&self) -> crate::error::Result<TrackedWorker> {
            let worker = TrackedWorker(Arc::clone(&self.drops));
            if self.fail {
                return Err(InternalError::PlanInvariant {
                    detail: "prepare failure".to_owned(),
                }
                .into());
            }
            Ok(worker)
        }
    }
    impl PrepareExecutionStep for TrackedDefinition {
        fn prepare(&self) -> crate::error::Result<Box<dyn super::ExecutionStep>> {
            Ok(Box::new(self.worker()?))
        }
    }
    impl PrepareExecutionStep<dyn ReadSetTransform> for TrackedDefinition {
        fn prepare(&self) -> crate::error::Result<Box<dyn ReadSetTransform>> {
            Ok(Box::new(self.worker()?))
        }
    }

    #[test]
    fn repeated_preparation_failures_clean_up_both_stage_kinds() {
        for fail_read_set in [false, true] {
            let drops = Arc::new(AtomicUsize::new(0));
            let definition = |fail| TrackedDefinition {
                drops: Arc::clone(&drops),
                fail,
            };
            let plan = if fail_read_set {
                Plan::<Logical>::new()
                    .step(definition(false))
                    .read_set_transform(definition(true))
                    .step(definition(false))
            } else {
                Plan::<Logical>::new()
                    .read_set_transform(definition(false))
                    .step(definition(true))
                    .read_set_transform(definition(false))
            };
            for attempt in 1..=2 {
                assert!(
                    matches!(plan.compile(), Err(RunError::Internal(InternalError::PlanInvariant { detail })) if detail == "prepare failure")
                );
                assert_eq!(drops.load(Ordering::SeqCst), attempt * 2);
            }
        }
    }

    #[test]
    fn logical_plan_accumulates_steps_and_compiles() {
        let _plan = Plan::<Logical>::new()
            .step(MinLength::new(4))
            .step(TrimPrefix::new(1))
            .compile()
            .expect("logical operations should prepare");
    }

    #[test]
    fn compiled_executions_outlive_the_reusable_logical_plan() {
        let recipe = Plan::<Logical>::new()
            .step(TrimPrefix::new(1))
            .step(MinLength::new(3));
        let mut first = recipe.compile().expect("first execution should prepare");
        let mut second = recipe.compile().expect("second execution should prepare");
        drop(recipe);

        for execution in [&mut first, &mut second] {
            let mut arena = TransformArena::new();
            let mut observer = RunObserver::new("test".to_owned());
            let outcome = execution
                .execute(record(b"ACGT"), &mut arena, &mut observer)
                .expect("prepared execution should own its operations");
            let sequences = outcome
                .emitted()
                .map(|record| record.sequence())
                .collect::<Vec<_>>();
            assert_eq!(sequences, [b"CGT".as_slice()]);
            assert_eq!(observer.reads_rejected, 0);
            assert_eq!(observer.transform_counts["trim_prefix"], 1);
        }
    }

    #[test]
    fn compile_groups_adjacent_per_unit_steps_around_read_set_transforms() {
        let plan = Plan::<Logical>::new()
            .step(MinLength::new(4))
            .step(TrimPrefix::new(1))
            .read_set_transform(FakeReadSetTransform)
            .step(MinLength::new(3))
            .compile()
            .expect("all stages should prepare");

        assert_eq!(plan.physical_items.len(), 3);
        assert!(matches!(
            plan.physical_items[0],
            PhysicalPlanItem::PerUnitChain(_)
        ));
        assert!(matches!(
            plan.physical_items[1],
            PhysicalPlanItem::ReadSet(_)
        ));
        assert!(matches!(
            plan.physical_items[2],
            PhysicalPlanItem::PerUnitChain(_)
        ));
    }

    #[test]
    fn streaming_only_plan_compiles_to_one_per_unit_chain() {
        let plan = Plan::<Logical>::new()
            .step(MinLength::new(4))
            .step(TrimPrefix::new(1))
            .compile()
            .expect("streaming operations should prepare");

        assert_eq!(plan.physical_items.len(), 1);
        assert!(matches!(
            plan.physical_items[0],
            PhysicalPlanItem::PerUnitChain(_)
        ));
    }

    #[test]
    fn per_unit_execution_rejects_read_set_transforms() {
        let mut plan = Plan::<Logical>::new()
            .read_set_transform(FakeReadSetTransform)
            .compile()
            .expect("read-set operation should prepare");
        let mut arena = TransformArena::new();
        let mut observer = RunObserver::new("test".to_owned());

        let Err(error) = plan.execute(record(b"ACGT"), &mut arena, &mut observer) else {
            panic!("per-unit execution should reject read-set transforms");
        };

        assert!(matches!(
            error,
            RunError::Internal(InternalError::PlanInvariant { detail })
                if detail.contains("fake_read_set")
        ));
    }

    #[test]
    fn single_execution_rejects_record_when_filter_fails() {
        let mut plan = Plan::<Logical>::new()
            .step(MinLength::new(6))
            .compile()
            .expect("filter should prepare");
        let mut arena = TransformArena::new();
        let mut observer = RunObserver::new("test".to_owned());

        let outcome = plan
            .execute(record(b"ACGT"), &mut arena, &mut observer)
            .expect("single-record filter failure should produce a rejected outcome");

        assert_eq!(outcome.emitted().count(), 0);
        assert_eq!(outcome.rejection_count(), 1);
        assert_eq!(observer.rejection_counts.get("too_short"), Some(&1));
    }

    #[test]
    fn single_execution_applies_transform_in_order() {
        let mut plan = Plan::<Logical>::new()
            .step(TrimPrefix::new(1))
            .compile()
            .expect("transform should prepare");
        let mut arena = TransformArena::new();
        let mut observer = RunObserver::new("test".to_owned());

        let outcome = plan
            .execute(record(b"ACGT"), &mut arena, &mut observer)
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
        let mut plan = Plan::<Logical>::new()
            .step(MinLength::new(6))
            .compile()
            .expect("filter should prepare");
        let mut arena = TransformArena::new();
        let mut observer = RunObserver::new("test".to_owned());

        let outcome = plan
            .execute(pair(b"ACGTAC", b"ACGT"), &mut arena, &mut observer)
            .expect("paired filter failure should produce a rejected outcome");

        assert_eq!(outcome.emitted().count(), 0);
        assert_eq!(outcome.rejection_count(), 1);
        assert_eq!(observer.reads_rejected, 1);
        assert_eq!(observer.reads_suppressed, 1);
    }

    #[test]
    fn paired_execution_can_emit_orphan_when_policy_allows() {
        let mut plan = Plan::<Logical>::new()
            .step(MinLength::new(6))
            .orphan_policy(OrphanPolicy::EmitOrphan)
            .compile()
            .expect("filter should prepare");
        let mut arena = TransformArena::new();
        let mut observer = RunObserver::new("test".to_owned());

        let outcome = plan
            .execute(pair(b"ACGTAC", b"ACGT"), &mut arena, &mut observer)
            .expect("paired orphan policy should produce an emitted orphan outcome");

        assert_eq!(outcome.emitted().count(), 1);
        assert_eq!(outcome.rejection_count(), 1);
        assert_eq!(observer.reads_suppressed, 0);
    }

    #[test]
    fn default_pair_filter_evaluates_both_mates_under_both_policies() {
        for policy in [OrphanPolicy::DropPair, OrphanPolicy::EmitOrphan] {
            for (left, right) in [
                (b"ACGT".as_slice(), b"ACGT".as_slice()),
                (b"ACGT".as_slice(), b"ACGTAC".as_slice()),
                (b"ACGTAC".as_slice(), b"ACGT".as_slice()),
                (b"ACGTAC".as_slice(), b"ACGTAC".as_slice()),
            ] {
                let mut plan = Plan::<Logical>::new()
                    .step(MinLength::new(6))
                    .orphan_policy(policy)
                    .compile()
                    .expect("prepare");
                let mut observer = RunObserver::new("test".to_owned());
                let mut arena = TransformArena::new();
                let outcome = plan
                    .execute(pair(left, right), &mut arena, &mut observer)
                    .expect("execute");
                let failures = usize::from(left.len() < 6) + usize::from(right.len() < 6);
                let emitted = if failures == 0 {
                    2
                } else {
                    usize::from(failures == 1 && policy == OrphanPolicy::EmitOrphan)
                };
                let records: Vec<_> = outcome.emitted().map(|r| r.sequence()).collect();
                assert_eq!(records.len(), emitted);
                if emitted == 1 {
                    assert_eq!(records[0], if left.len() == 6 { left } else { right });
                }
                assert_eq!(outcome.rejection_count(), failures);
                assert_eq!(observer.reads_rejected, failures as u64);
                assert_eq!(
                    observer
                        .rejection_counts
                        .get("too_short")
                        .copied()
                        .unwrap_or(0),
                    failures as u64
                );
                assert_eq!(
                    observer.reads_suppressed,
                    u64::from(failures == 1 && policy == OrphanPolicy::DropPair)
                );
            }
        }
    }

    struct CollapsePair;
    impl super::PairTransform for CollapsePair {
        fn code(&self) -> &'static str {
            "collapse"
        }
        fn apply_pair<'a>(
            &mut self,
            pair: RecordPair<'a>,
            _arena: &'a TransformArena,
        ) -> crate::error::Result<super::PairTransformResult<'a>> {
            Ok(super::PairTransformResult::Single {
                record: pair.left(),
                applied: true,
            })
        }
    }
    impl PrepareExecutionStep for CollapsePair {
        fn prepare(&self) -> crate::error::Result<Box<dyn super::ExecutionStep>> {
            Ok(Box::new(super::PairTransformStep(Self)))
        }
    }
    struct DropUnit;
    impl super::PairTransform for DropUnit {
        fn code(&self) -> &'static str {
            "drop"
        }
        fn apply_pair<'a>(
            &mut self,
            _pair: RecordPair<'a>,
            _arena: &'a TransformArena,
        ) -> crate::error::Result<super::PairTransformResult<'a>> {
            Ok(super::PairTransformResult::Drop { reason: "dropped" })
        }
        fn apply_single<'a>(
            &mut self,
            _record: RecordView<'a>,
            _arena: &'a TransformArena,
        ) -> crate::error::Result<super::PairTransformResult<'a>> {
            Ok(super::PairTransformResult::Drop { reason: "dropped" })
        }
    }
    impl PrepareExecutionStep for DropUnit {
        fn prepare(&self) -> crate::error::Result<Box<dyn super::ExecutionStep>> {
            Ok(Box::new(super::PairTransformStep(Self)))
        }
    }

    #[test]
    fn merged_single_rejection_counts_only_current_record() {
        let mut plan = Plan::<Logical>::new()
            .step(CollapsePair)
            .step(MinLength::new(6))
            .step(TrimPrefix::new(1))
            .compile()
            .expect("prepare");
        let mut arena = TransformArena::new();
        let mut observer = RunObserver::new("test".to_owned());
        let outcome = plan
            .execute(pair(b"ACGT", b"ACGTAC"), &mut arena, &mut observer)
            .expect("execute");
        assert_eq!(outcome.emitted().count(), 0);
        assert_eq!(outcome.rejection_count(), 1);
        assert_eq!(observer.reads_rejected, 1);
        assert_eq!(observer.rejection_counts["too_short"], 1);
        assert_eq!(observer.reads_suppressed, 0);
        assert!(!observer.transform_counts.contains_key("trim_prefix"));
    }

    #[test]
    fn pair_transform_drop_counts_current_unit_cardinality() {
        for single in [false, true] {
            let recipe = Plan::<Logical>::new()
                .step(DropUnit)
                .step(TrimPrefix::new(1));
            let mut plan = recipe.compile().expect("prepare");
            let mut arena = TransformArena::new();
            let mut observer = RunObserver::new("test".to_owned());
            let outcome = if single {
                plan.execute(record(b"ACGT"), &mut arena, &mut observer)
            } else {
                plan.execute(pair(b"ACGT", b"ACGT"), &mut arena, &mut observer)
            }
            .expect("execute");
            assert_eq!(outcome.emitted().count(), 0);
            assert_eq!(outcome.rejection_count(), 1);
            let expected = if single { 1 } else { 2 };
            assert_eq!(observer.reads_rejected, expected);
            assert_eq!(observer.rejection_counts["dropped"], expected);
            assert_eq!(observer.reads_suppressed, 0);
            assert!(observer.transform_counts.is_empty());
        }
    }

    struct AtomicFilter(bool);

    impl ReadFilter for AtomicFilter {
        type Reason = TooShort;

        fn evaluate(&mut self, _record: &RecordView<'_>) -> Result<(), Self::Reason> {
            panic!("pair hook must be used")
        }

        fn evaluate_pair(&mut self, _pair: &RecordPair<'_>) -> PairRelation<Self::Reason> {
            PairRelation::WholePair(if self.0 { Ok(()) } else { Err(TooShort) })
        }
    }

    impl PrepareExecutionStep for AtomicFilter {
        fn prepare(&self) -> crate::error::Result<Box<dyn super::ExecutionStep>> {
            Ok(Box::new(FilterStep(Self(self.0))))
        }
    }

    #[test]
    fn whole_pair_rejection_is_atomic_under_both_orphan_policies() {
        for policy in [OrphanPolicy::DropPair, OrphanPolicy::EmitOrphan] {
            let mut plan = Plan::<Logical>::new()
                .step(AtomicFilter(false))
                .step(TrimPrefix::new(1))
                .orphan_policy(policy)
                .compile()
                .expect("prepare");
            let mut arena = TransformArena::new();
            let mut observer = RunObserver::new("test".to_owned());
            let outcome = plan
                .execute(pair(b"ACGTAC", b"ACGT"), &mut arena, &mut observer)
                .expect("execute");
            assert_eq!(outcome.emitted().count(), 0);
            assert_eq!(outcome.rejection_count(), 1);
            assert_eq!(observer.reads_rejected, 2);
            assert_eq!(observer.rejection_counts["too_short"], 2);
            assert_eq!(observer.reads_suppressed, 0);
            assert!(observer.transform_counts.is_empty());
        }
    }

    #[test]
    fn whole_pair_acceptance_preserves_both_mates() {
        for policy in [OrphanPolicy::DropPair, OrphanPolicy::EmitOrphan] {
            let mut plan = Plan::<Logical>::new()
                .step(AtomicFilter(true))
                .orphan_policy(policy)
                .compile()
                .expect("prepare");
            let mut arena = TransformArena::new();
            let mut observer = RunObserver::new("test".to_owned());
            let outcome = plan
                .execute(pair(b"ACGTAC", b"ACGT"), &mut arena, &mut observer)
                .expect("execute");
            assert_eq!(
                outcome.emitted().map(|r| r.sequence()).collect::<Vec<_>>(),
                [b"ACGTAC".as_slice(), b"ACGT".as_slice()]
            );
            assert_eq!(outcome.rejection_count(), 0);
            assert_eq!(observer.reads_rejected, 0);
            assert_eq!(observer.reads_suppressed, 0);
        }
    }

    struct StatefulFilter(Cell<usize>);

    impl ReadFilter for StatefulFilter {
        type Reason = TooShort;
        fn evaluate(&mut self, _record: &RecordView<'_>) -> Result<(), Self::Reason> {
            let count = self.0.get();
            self.0.set(count + 1);
            if count == 0 { Ok(()) } else { Err(TooShort) }
        }
    }

    struct StatefulDefinition;
    impl PrepareExecutionStep for StatefulDefinition {
        fn prepare(&self) -> crate::error::Result<Box<dyn super::ExecutionStep>> {
            Ok(Box::new(FilterStep(StatefulFilter(Cell::new(0)))))
        }
    }

    #[test]
    fn prepared_non_sync_workers_have_independent_state_and_accept_borrowed_records() {
        let recipe = Plan::<Logical>::new().step(StatefulDefinition);
        let mut first = recipe.compile().expect("prepare");
        let mut second = recipe.compile().expect("prepare");
        let sequence = Vec::from(&b"ACGT"[..]);
        let quality = Vec::from(&b"IIII"[..]);
        let record = RecordView::new(b"read", &sequence, &quality);
        for plan in [&mut first, &mut second] {
            let mut observer = RunObserver::new("test".to_owned());
            let mut arena = TransformArena::new();
            assert_eq!(
                plan.execute(record, &mut arena, &mut observer)
                    .expect("execute")
                    .emitted()
                    .count(),
                1
            );
            let mut arena = TransformArena::new();
            assert_eq!(
                plan.execute(record, &mut arena, &mut observer)
                    .expect("execute")
                    .emitted()
                    .count(),
                0
            );
        }
    }

    #[test]
    fn emitted_iterator_preserves_slot_order() {
        let outcome = ExecutionOutcome {
            emitted: super::EmittedUnit::Pair(pair(b"AAAAAA", b"CCCCCC")),
            rejection_count: 0,
        };

        let emitted = outcome
            .emitted()
            .map(|record| record.sequence())
            .collect::<Vec<_>>();
        assert_eq!(emitted, vec![b"AAAAAA".as_slice(), b"CCCCCC".as_slice()]);
    }
}
