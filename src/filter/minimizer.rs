//! Deacon-backed selection with reusable immutable indexes and per-plan workers.

use std::{
    fs::File,
    io::{self, BufReader, Read},
    path::Path,
    sync::Arc,
};

use deacon::{FilterKernel, FilterParams, Index};

use crate::{
    error::{InternalError, Result, UnavailableInputError, UsageError, constructors},
    plan::{
        ExecutionStep, FilterStep, PairRelation, PrepareExecutionStep, ReadFilter, RejectionReason,
    },
    record::{RecordPair, RecordView},
};

/// How enrichment and depletion decisions are combined.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum SelectionPolicy {
    /// Require enrichment and reject depletion.
    #[default]
    DepletionVeto,
    /// Keep enrichment hits or reads without depletion hits.
    EnrichmentRescue,
}

#[derive(Clone)]
struct SelectionIndex {
    index: Arc<Index>,
    min_hits: usize,
}

/// Immutable selection definition shared by independently prepared workers.
#[derive(Clone)]
pub(crate) struct MinimizerFilter {
    enrich: Option<SelectionIndex>,
    deplete: Option<SelectionIndex>,
    policy: SelectionPolicy,
}

impl MinimizerFilter {
    /// Configure index paths, thresholds, and policy through named builder methods.
    pub(crate) fn builder<'a>() -> MinimizerFilterBuilder<'a> {
        MinimizerFilterBuilder::default()
    }

    /// Load explicitly configured indexes into a reusable filter definition.
    ///
    /// Returns `None` when both paths are absent. Thresholds for supplied indexes must be positive.
    /// Prefer [`Self::builder`] when assembling configuration from individual settings.
    pub(crate) fn load(
        enrich_path: Option<&Path>,
        enrich_min_hits: usize,
        deplete_path: Option<&Path>,
        deplete_min_hits: usize,
        policy: SelectionPolicy,
    ) -> Result<Option<Self>> {
        for (path, hits, role) in [
            (enrich_path, enrich_min_hits, "enrichment"),
            (deplete_path, deplete_min_hits, "depletion"),
        ] {
            if path.is_some() && hits == 0 {
                return Err(UsageError::ZeroSelectionThreshold { role }.into());
            }
        }
        let enrich = enrich_path
            .map(|path| load_index(path, "enrichment", enrich_min_hits))
            .transpose()?;
        let deplete = deplete_path
            .map(|path| load_index(path, "depletion", deplete_min_hits))
            .transpose()?;
        if enrich.is_none() && deplete.is_none() {
            return Ok(None);
        }
        Ok(Some(Self {
            enrich,
            deplete,
            policy,
        }))
    }
}

/// Borrowed index paths and optional overrides for a minimizer filter.
#[derive(Default)]
#[must_use]
pub(crate) struct MinimizerFilterBuilder<'a> {
    enrich_index: Option<&'a Path>,
    enrich_min_hits: Option<usize>,
    deplete_index: Option<&'a Path>,
    deplete_min_hits: Option<usize>,
    policy: Option<SelectionPolicy>,
}

impl<'a> MinimizerFilterBuilder<'a> {
    /// Set the enrichment index path.
    pub(crate) fn enrich_index(mut self, path: impl Into<Option<&'a Path>>) -> Self {
        self.enrich_index = path.into();
        self
    }

    /// Set the enrichment hit threshold; `None` uses the default of two hits.
    pub(crate) fn enrich_min_hits(mut self, hits: impl Into<Option<usize>>) -> Self {
        self.enrich_min_hits = hits.into();
        self
    }

    /// Set the depletion index path.
    pub(crate) fn deplete_index(mut self, path: impl Into<Option<&'a Path>>) -> Self {
        self.deplete_index = path.into();
        self
    }

    /// Set the depletion hit threshold; `None` uses the default of two hits.
    pub(crate) fn deplete_min_hits(mut self, hits: impl Into<Option<usize>>) -> Self {
        self.deplete_min_hits = hits.into();
        self
    }

    /// Set the combination policy; `None` uses depletion veto.
    pub(crate) fn policy(mut self, policy: impl Into<Option<SelectionPolicy>>) -> Self {
        self.policy = policy.into();
        self
    }

    /// Load the configured indexes into an owned filter definition.
    pub(crate) fn load(self) -> Result<Option<MinimizerFilter>> {
        MinimizerFilter::load(
            self.enrich_index,
            self.enrich_min_hits.unwrap_or(2),
            self.deplete_index,
            self.deplete_min_hits.unwrap_or(2),
            self.policy.unwrap_or_default(),
        )
    }
}

fn load_index(path: &Path, role: &'static str, min_hits: usize) -> Result<SelectionIndex> {
    let file = File::open(path).map_err(|source| UnavailableInputError::SelectionIndex {
        role,
        path: path.to_path_buf(),
        source,
    })?;
    let mut reader = BufReader::with_capacity(
        1 << 20,
        ErrorRecordingReader {
            inner: file,
            error_kind: None,
        },
    );
    let index = deacon::load_index(&mut reader).map_err(|source| {
        constructors::index_load_error(
            role,
            path,
            reader.get_ref().error_kind,
            source.into_boxed_dyn_error(),
        )
    })?;
    Ok(SelectionIndex {
        index: Arc::new(index),
        min_hits,
    })
}

/// Wraps a reader so index loading can distinguish read failures from malformed data.
///
/// With the pinned Deacon version, Bincode keeps some underlying I/O errors out
/// of the standard error chain. Walking `Error::source()` can therefore miss a
/// real read failure and misreport it as a corrupt index. Inspecting Bincode's
/// concrete error type would require a direct dependency just for this check.
///
/// Instead, this wrapper records the error kind at read time and leaves the read
/// result unchanged. [`constructors::index_load_error`] uses that record to
/// classify a failed load while retaining Deacon's original error context.
///
/// A successful read clears the record. Otherwise, an interrupted read that the
/// decoder retried successfully could cause a later malformed-index error to be
/// reported as an I/O failure.
struct ErrorRecordingReader<R> {
    inner: R,
    error_kind: Option<io::ErrorKind>,
}

impl<R: Read> Read for ErrorRecordingReader<R> {
    fn read(&mut self, buffer: &mut [u8]) -> io::Result<usize> {
        let result = self.inner.read(buffer);
        // Clear failures recovered by the decoder, including retried interrupted reads.
        self.error_kind = result.as_ref().err().map(io::Error::kind);
        result
    }
}

fn prepare_index(index: &SelectionIndex) -> Result<FilterKernel> {
    FilterKernel::new(
        Arc::clone(&index.index),
        FilterParams {
            deplete: false,
            abs_threshold: index.min_hits,
            rel_threshold: 0.0,
            prefix_length: 0,
        },
    )
    .map_err(|source| {
        InternalError::SelectionWorker {
            source: source.into_boxed_dyn_error(),
        }
        .into()
    })
}

impl PrepareExecutionStep for MinimizerFilter {
    fn prepare(&self) -> Result<Box<dyn ExecutionStep>> {
        Ok(Box::new(FilterStep(MinimizerWorker {
            enrich: self.enrich.as_ref().map(prepare_index).transpose()?,
            deplete: self.deplete.as_ref().map(prepare_index).transpose()?,
            policy: self.policy,
        })))
    }
}

struct MinimizerWorker {
    enrich: Option<FilterKernel>,
    deplete: Option<FilterKernel>,
    policy: SelectionPolicy,
}

#[derive(Debug)]
struct NotSelected;
impl RejectionReason for NotSelected {
    fn code(&self) -> &'static str {
        "not_selected"
    }
}

impl MinimizerWorker {
    fn keep(&mut self, mut matches: impl FnMut(&mut FilterKernel) -> bool) -> bool {
        match (&mut self.enrich, &mut self.deplete) {
            (Some(enrich), Some(deplete)) => match self.policy {
                SelectionPolicy::DepletionVeto => matches(enrich) && !matches(deplete),
                SelectionPolicy::EnrichmentRescue => matches(enrich) || !matches(deplete),
            },
            (Some(enrich), None) => matches(enrich),
            (None, Some(deplete)) => !matches(deplete),
            (None, None) => true,
        }
    }
}

impl ReadFilter for MinimizerWorker {
    type Reason = NotSelected;
    fn evaluate(&mut self, record: &RecordView<'_>) -> std::result::Result<(), Self::Reason> {
        self.keep(|kernel| kernel.classify_read(record.sequence()).is_match)
            .then_some(())
            .ok_or(NotSelected)
    }
    fn evaluate_pair(&mut self, pair: &RecordPair<'_>) -> PairRelation<Self::Reason> {
        PairRelation::WholePair(
            self.keep(|kernel| {
                kernel
                    .classify_pair(pair.left().sequence(), pair.right().sequence())
                    .is_match
            })
            .then_some(())
            .ok_or(NotSelected),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::{ErrorRecordingReader, MinimizerFilter, SelectionPolicy};
    use crate::{
        error::{Result, RunError, UsageError, constructors},
        observer::RunObserver,
        plan::{BuildPlan, Execute, Logical, OrphanPolicy, Plan, TransformArena},
        record::{RecordPair, RecordView},
    };
    use deacon::{Index, MinimizerVec};
    use std::{
        io::{self, BufReader, Read},
        path::{Path, PathBuf},
        process::ExitCode,
    };

    struct ErrorInjectingReader<'a> {
        before: &'a [u8],
        error: Option<io::ErrorKind>,
        after: &'a [u8],
    }

    impl Read for ErrorInjectingReader<'_> {
        fn read(&mut self, buffer: &mut [u8]) -> io::Result<usize> {
            if buffer.is_empty() || !self.before.is_empty() {
                return self.before.read(buffer);
            }
            if let Some(kind) = self.error.take() {
                return Err(io::Error::from(kind));
            }
            self.after.read(buffer)
        }
    }

    fn decode_index(reader: impl Read) -> Result<Index> {
        let mut reader = BufReader::new(ErrorRecordingReader {
            inner: reader,
            error_kind: None,
        });
        deacon::load_index(&mut reader).map_err(|source| {
            constructors::index_load_error(
                "enrichment",
                Path::new("enrich.idx"),
                reader.get_ref().error_kind,
                source.into_boxed_dyn_error(),
            )
        })
    }

    fn index_bytes() -> Vec<u8> {
        let index = Index::from_minimizers(3, 1, MinimizerVec::U64(vec![0, 21])).expect("index");
        let mut bytes = Vec::new();
        deacon::write_index(&index, &mut bytes).expect("serialize index");
        bytes
    }

    #[test]
    fn index_load_classifies_read_failures_at_every_offset() {
        let bytes = index_bytes();
        for (kind, expected_status) in [
            (io::ErrorKind::PermissionDenied, 74),
            (io::ErrorKind::Other, 74),
            (io::ErrorKind::UnexpectedEof, 65),
        ] {
            for offset in 0..bytes.len() {
                let error = decode_index(ErrorInjectingReader {
                    before: &bytes[..offset],
                    error: Some(kind),
                    after: &[],
                })
                .expect_err("read failure");
                assert_eq!(
                    error.exit_code(),
                    ExitCode::from(expected_status),
                    "{kind:?} at {offset}",
                );
            }
        }
    }

    #[test]
    fn recovered_interruptions_preserve_index_load_outcomes() {
        let valid = index_bytes();
        for (bytes, succeeds) in [
            (valid.as_slice(), true),
            (b"bad".as_slice(), false),
            (&valid[..2], false),
        ] {
            for offset in 1..bytes.len() {
                let result = decode_index(ErrorInjectingReader {
                    before: &bytes[..offset],
                    error: Some(io::ErrorKind::Interrupted),
                    after: &bytes[offset..],
                });
                if succeeds {
                    assert!(result.is_ok(), "interrupted at {offset}");
                } else {
                    let error = result.expect_err("malformed index");
                    assert_eq!(error.exit_code(), ExitCode::from(65), "at {offset}");
                }
            }
        }

        // Deacon's initial fill_buf does not retry an interrupted read.
        let error = decode_index(ErrorInjectingReader {
            before: &[],
            error: Some(io::ErrorKind::Interrupted),
            after: &valid,
        })
        .expect_err("unrecovered interruption");
        assert_eq!(error.exit_code(), ExitCode::from(74));
    }

    // Exact indexes constructed with the pinned Deacon library API, not an independent producer.
    fn index(dir: &Path, name: &str, values: Vec<u64>) -> PathBuf {
        let path = dir.join(name);
        let index = Index::from_minimizers(3, 1, MinimizerVec::U64(values)).expect("index");
        deacon::write_index_to_path(&index, &path).expect("write index");
        path
    }

    #[test]
    fn builder_defaults_and_raw_constructor_validate_before_loading() {
        let dir = tempfile::tempdir().expect("tempdir");
        let path = index(dir.path(), "index", vec![0, 21]);
        let built = MinimizerFilter::builder()
            .enrich_index(path.as_path())
            .load()
            .expect("builder")
            .expect("filter");
        let raw = MinimizerFilter::load(Some(&path), 2, None, 2, SelectionPolicy::DepletionVeto)
            .expect("raw")
            .expect("filter");
        for filter in [built, raw] {
            let mut plan = Plan::<Logical>::new()
                .step(filter)
                .compile()
                .expect("compile");
            let record = RecordView::new(b"read", b"AAA", b"III");
            let mut arena = TransformArena::new();
            let mut observer = RunObserver::new("test".to_owned());
            assert_eq!(
                plan.execute(record, &mut arena, &mut observer)
                    .expect("execute")
                    .emitted()
                    .count(),
                0
            );
        }
        let missing = dir.path().join("missing");
        for result in [
            MinimizerFilter::load(
                Some(&missing),
                1,
                Some(&missing),
                0,
                SelectionPolicy::default(),
            ),
            MinimizerFilter::builder()
                .enrich_index(missing.as_path())
                .enrich_min_hits(0)
                .load(),
        ] {
            assert!(matches!(
                result,
                Err(RunError::Usage(UsageError::ZeroSelectionThreshold { .. }))
            ));
        }
        assert!(
            MinimizerFilter::load(None, 0, None, 0, SelectionPolicy::default())
                .expect("empty")
                .is_none()
        );
    }

    #[test]
    fn selection_truth_table_and_independent_compilations() {
        let dir = tempfile::tempdir().expect("tempdir");
        let enrich = index(dir.path(), "enrich", vec![0]);
        let deplete = index(dir.path(), "deplete", vec![21]);
        for policy in [
            SelectionPolicy::DepletionVeto,
            SelectionPolicy::EnrichmentRescue,
        ] {
            for (e_index, d_index) in [(false, false), (true, false), (false, true), (true, true)] {
                let selection = MinimizerFilter::builder()
                    .enrich_index(e_index.then_some(enrich.as_path()))
                    .deplete_index(d_index.then_some(deplete.as_path()))
                    .enrich_min_hits(1)
                    .deplete_min_hits(1)
                    .policy(policy)
                    .load()
                    .expect("load");
                let recipe = if let Some(selection) = selection {
                    Plan::<Logical>::new().step(selection)
                } else {
                    Plan::<Logical>::new()
                };
                let mut first = recipe.compile().expect("first worker");
                let mut second = recipe.compile().expect("second worker");
                drop(recipe);
                for plan in [&mut first, &mut second] {
                    for (seq, e, d) in [
                        (b"NNN".as_slice(), false, false),
                        (b"AAA", true, false),
                        (b"CCC", false, true),
                        (b"AAANCCC", true, true),
                    ] {
                        let keep = match (e_index, d_index, policy) {
                            (false, false, _) => true,
                            (true, false, _) => e,
                            (false, true, _) => !d,
                            (true, true, SelectionPolicy::DepletionVeto) => e && !d,
                            (true, true, SelectionPolicy::EnrichmentRescue) => e || !d,
                        };
                        let quality = vec![b'I'; seq.len()];
                        let record = RecordView::new(b"read", seq, &quality);
                        let mut arena = TransformArena::new();
                        let mut observer = RunObserver::new("test".to_owned());
                        let outcome = plan
                            .execute(record, &mut arena, &mut observer)
                            .expect("execute");
                        assert_eq!(
                            outcome.emitted().count(),
                            usize::from(keep),
                            "{policy:?} {e_index} {d_index} {seq:?}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn pairs_pool_distinct_evidence_and_are_atomic() {
        let dir = tempfile::tempdir().expect("tempdir");
        let path = index(dir.path(), "enrich", vec![0, 21]);
        for policy in [OrphanPolicy::DropPair, OrphanPolicy::EmitOrphan] {
            let recipe = Plan::<Logical>::new()
                .step(
                    MinimizerFilter::builder()
                        .enrich_index(path.as_path())
                        .load()
                        .expect("load")
                        .expect("filter"),
                )
                .orphan_policy(policy);
            let mut plan = recipe.compile().expect("worker");
            for (left, right, keep) in [
                (b"AAA".as_slice(), b"CCC".as_slice(), true),
                (b"CCC", b"AAA", true),
                (b"AAA", b"AAA", false),
                (b"AAA", b"NNN", false),
            ] {
                let left = RecordView::new(b"read", left, b"III");
                let right = RecordView::new(b"read", right, b"III");
                let pair = RecordPair::try_new(left.into(), right.into()).expect("pair");
                let mut arena = TransformArena::new();
                let mut observer = RunObserver::new("test".to_owned());
                let outcome = plan
                    .execute(pair, &mut arena, &mut observer)
                    .expect("execute");
                assert_eq!(outcome.emitted().count(), if keep { 2 } else { 0 });
                assert_eq!(outcome.rejection_count(), usize::from(!keep));
                assert_eq!(observer.reads_rejected, if keep { 0 } else { 2 });
                assert_eq!(observer.reads_suppressed, 0);
            }
        }
    }
}
