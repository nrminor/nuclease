//! Generated fragment histories against an independent, two-class evidence model.

use std::path::{Path, PathBuf};

use deacon::{Index, MinimizerVec};
use proptest::{prelude::*, test_runner::FileFailurePersistence};

use crate::{
    filter::minimizer::{MinimizerFilter, SelectionPolicy},
    observer::RunObserver,
    plan::{BuildPlan, Execute, Execution, Logical, OrphanPolicy, Plan, TransformArena},
    record::{RecordPair, RecordView},
};

#[derive(Clone, Debug)]
struct Settings {
    enrich: Option<u8>,
    deplete: Option<u8>,
    enrich_threshold: usize,
    deplete_threshold: usize,
    policy: SelectionPolicy,
    orphan_policy: OrphanPolicy,
}

#[derive(Clone, Debug)]
enum Unit {
    Single(Vec<u8>),
    Pair(Vec<u8>, Vec<u8>),
}

impl Unit {
    fn tokens(&self) -> (&[u8], Option<&[u8]>) {
        match self {
            Self::Single(left) => (left, None),
            Self::Pair(left, right) => (left, Some(right)),
        }
    }
}

fn settings() -> impl Strategy<Value = Settings> {
    (
        proptest::option::of(0_u8..=3),
        proptest::option::of(0_u8..=3),
        1_usize..=3,
        1_usize..=3,
        any::<bool>(),
        any::<bool>(),
    )
        .prop_map(
            |(enrich, deplete, enrich_threshold, deplete_threshold, rescue, orphan)| Settings {
                enrich,
                deplete,
                enrich_threshold,
                deplete_threshold,
                policy: if rescue {
                    SelectionPolicy::EnrichmentRescue
                } else {
                    SelectionPolicy::DepletionVeto
                },
                orphan_policy: if orphan {
                    OrphanPolicy::EmitOrphan
                } else {
                    OrphanPolicy::DropPair
                },
            },
        )
}

fn unit() -> impl Strategy<Value = Unit> {
    let tokens = proptest::collection::vec(0_u8..4, 0..=8);
    prop_oneof![
        tokens.clone().prop_map(Unit::Single),
        (tokens.clone(), tokens).prop_map(|(left, right)| Unit::Pair(left, right)),
    ]
}

// Token identity is the oracle: AAA/TTT contribute bit 0, CCC/GGG bit 1.
// No Deacon classifier, minimizer extraction, or production policy helper is used.
fn evidence(tokens: &[u8]) -> u8 {
    tokens
        .iter()
        .fold(0, |mask, token| mask | (1 << (token / 2)))
}

fn expected_keep(settings: &Settings, unit: &Unit) -> bool {
    let (left, right) = unit.tokens();
    let query = evidence(left) | right.map_or(0, evidence);
    let matches = |index: u8, threshold: usize| {
        [1, 2]
            .into_iter()
            .filter(|bit| query & index & bit != 0)
            .count()
            >= threshold
    };
    let enrich = settings
        .enrich
        .is_some_and(|index| matches(index, settings.enrich_threshold));
    let deplete = settings
        .deplete
        .is_some_and(|index| matches(index, settings.deplete_threshold));
    match (settings.enrich, settings.deplete) {
        (None, None) => true,
        (Some(_), None) => enrich,
        (None, Some(_)) => !deplete,
        (Some(_), Some(_)) => {
            // Rows are (E,D) = 00, 01, 10, 11; evaluate both predicates eagerly.
            let row = 2 * usize::from(enrich) + usize::from(deplete);
            match settings.policy {
                SelectionPolicy::DepletionVeto => [false, false, true, false][row],
                SelectionPolicy::EnrichmentRescue => [true, false, true, true][row],
            }
        }
    }
}

fn render(tokens: &[u8]) -> Vec<u8> {
    let mut sequence = vec![b'N'];
    for token in tokens {
        sequence.extend_from_slice(match token {
            0 => b"AAA",
            1 => b"TTT",
            2 => b"CCC",
            3 => b"GGG",
            _ => unreachable!("token strategies only generate 0..4"),
        });
        // Ambiguity breaks windows, so adjacent tokens cannot create extra 3-mers.
        sequence.push(b'N');
    }
    sequence
}

fn write_index(dir: &Path, name: &str, mask: u8) -> PathBuf {
    let values = [(1, 0), (2, 21)]
        .into_iter()
        .filter_map(|(bit, value)| (mask & bit != 0).then_some(value))
        .collect();
    let index = Index::from_minimizers(3, 1, MinimizerVec::U64(values)).expect("exact index");
    let path = dir.join(name);
    deacon::write_index_to_path(&index, &path).expect("serialize exact index");
    path
}

fn check_unit(
    plan: &mut Plan<Execution>,
    observer: &mut RunObserver,
    settings: &Settings,
    unit: &Unit,
) {
    let (left, right) = unit.tokens();
    let left_sequence = render(left);
    let left_quality = vec![b'I'; left_sequence.len()];
    let right_sequence = right.map(render);
    let right_quality = right_sequence.as_ref().map(|seq| vec![b'J'; seq.len()]);
    let left = RecordView::new(b"fragment/1", &left_sequence, &left_quality);
    let right = right_sequence
        .as_ref()
        .zip(right_quality.as_ref())
        .map(|(seq, quality)| RecordView::new(b"fragment/2", seq, quality));
    let keep = expected_keep(settings, unit);
    let cardinality = if right.is_some() { 2_u64 } else { 1 };
    let before = (observer.reads_rejected, observer.reads_suppressed);
    let mut arena = TransformArena::new();
    let outcome = if let Some(right) = right {
        let pair = RecordPair::try_new(left.into(), right.into()).expect("valid generated pair");
        plan.execute(pair, &mut arena, observer)
    } else {
        plan.execute(left, &mut arena, observer)
    }
    .expect("execute generated unit");
    let snapshot = |record: RecordView<'_>| {
        (
            record.header().to_vec(),
            record.sequence().to_vec(),
            record.quality().to_vec(),
        )
    };
    let expected = if keep {
        std::iter::once(left)
            .chain(right)
            .map(snapshot)
            .collect::<Vec<_>>()
    } else {
        Vec::new()
    };
    assert_eq!(
        outcome.emitted().map(snapshot).collect::<Vec<_>>(),
        expected,
        "{settings:?}: {unit:?}",
    );
    assert_eq!(outcome.rejection_count(), usize::from(!keep), "{unit:?}");
    assert_eq!(
        observer.reads_rejected - before.0,
        if keep { 0 } else { cardinality },
        "{settings:?}: {unit:?}",
    );
    assert_eq!(observer.reads_suppressed - before.1, 0, "{unit:?}");
}

fn check_history(settings: &Settings, history: &[Unit]) {
    let dir = tempfile::tempdir().expect("index directory");
    let enrich = settings
        .enrich
        .map(|mask| write_index(dir.path(), "enrich", mask));
    let deplete = settings
        .deplete
        .map(|mask| write_index(dir.path(), "deplete", mask));
    let selection = MinimizerFilter::builder()
        .enrich_index(enrich.as_deref())
        .deplete_index(deplete.as_deref())
        .enrich_min_hits(settings.enrich_threshold)
        .deplete_min_hits(settings.deplete_threshold)
        .policy(settings.policy)
        .load()
        .expect("load generated indexes");
    let recipe = Plan::<Logical>::new().orphan_policy(settings.orphan_policy);
    let recipe = match selection {
        Some(selection) => recipe.step(selection),
        None => recipe,
    };
    let mut first = recipe.compile().expect("first independent worker");
    let mut second = recipe.compile().expect("second independent worker");
    drop(recipe);
    let mut first_observer = RunObserver::new("forward".to_owned());
    let mut second_observer = RunObserver::new("reverse".to_owned());
    for (forward, reverse) in history.iter().zip(history.iter().rev()) {
        check_unit(&mut first, &mut first_observer, settings, forward);
        check_unit(&mut second, &mut second_observer, settings, reverse);
    }
}

proptest! {
    #![proptest_config(ProptestConfig {
        cases: 128,
        max_shrink_iters: 2048,
        failure_persistence: Some(Box::new(FileFailurePersistence::SourceParallel(
            "proptest-regressions",
        ))),
        ..ProptestConfig::default()
    })]

    #[test]
    fn generated_histories_preserve_fragment_selection_and_worker_isolation(
        settings in settings(),
        history in proptest::collection::vec(unit(), 1..=12),
    ) {
        check_history(&settings, &history);
    }
}

#[test]
fn unequal_thresholds_and_duplicate_evidence_have_literal_witnesses() {
    let history = [
        Unit::Single(vec![2]),
        Unit::Pair(vec![0], vec![2]),
        Unit::Pair(vec![0], vec![0]),
        Unit::Pair(vec![0], vec![1]),
        Unit::Single(vec![0, 0, 1]),
        Unit::Single(vec![]),
        Unit::Pair(vec![2], vec![3]),
        Unit::Pair(vec![2], vec![0]),
    ];
    for (policy, expected) in [
        (SelectionPolicy::DepletionVeto, [false; 8]),
        (
            SelectionPolicy::EnrichmentRescue,
            [false, true, true, true, true, true, false, true],
        ),
    ] {
        for orphan_policy in [OrphanPolicy::DropPair, OrphanPolicy::EmitOrphan] {
            let settings = Settings {
                enrich: Some(3),
                deplete: Some(2),
                enrich_threshold: 2,
                deplete_threshold: 1,
                policy,
                orphan_policy,
            };
            for (unit, keep) in history.iter().zip(expected) {
                assert_eq!(expected_keep(&settings, unit), keep, "{unit:?}");
            }
            check_history(&settings, &history);
        }
    }
}

#[test]
fn pooled_evidence_and_empty_indexes_do_not_depend_on_history() {
    let history = [
        Unit::Pair(vec![0], vec![2]),
        Unit::Single(vec![0, 0, 1]),
        Unit::Pair(vec![0], vec![1]),
        Unit::Single(vec![]),
        Unit::Single(vec![0, 2]),
    ];
    // These literal outcomes distinguish absent roles, empty indexes and distinct hits.
    for (enrich, deplete, expected) in [
        (Some(3), None, [true, false, false, false, true]),
        (None, Some(3), [false, true, true, true, false]),
        (Some(0), None, [false; 5]),
        (None, Some(0), [true; 5]),
        (None, None, [true; 5]),
    ] {
        let settings = Settings {
            enrich,
            deplete,
            enrich_threshold: 2,
            deplete_threshold: 2,
            policy: SelectionPolicy::DepletionVeto,
            orphan_policy: OrphanPolicy::EmitOrphan,
        };
        for (unit, keep) in history.iter().zip(expected) {
            assert_eq!(expected_keep(&settings, unit), keep, "{unit:?}");
        }
        check_history(&settings, &history);
    }
}
