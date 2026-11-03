//! Offline interoperability tests using independently produced Deacon indexes.
//!
//! Producer: official Deacon 0.19.1 aarch64-apple-darwin release executable,
//! <https://github.com/bede/deacon/releases/download/0.19.1/deacon-aarch64-apple-darwin>
//! SHA-256: fc61da27a16c15fa4434649940039850c7fdd46d4e0880026d644502ea3430b1
//! From the repository root (with `deacon` denoting that executable):
//! ```text
//! deacon index build tests/fixtures/deacon-0.19.1/reference.fa -k 7 -w 3 -o tests/fixtures/deacon-0.19.1/exact-k7-w3.idx -q
//! deacon index freeze tests/fixtures/deacon-0.19.1/exact-k7-w3.idx -b 32 -o tests/fixtures/deacon-0.19.1/fuse32-k7-w3.pidx -q
//! deacon index build tests/fixtures/deacon-0.19.1/reference.fa -k 33 -w 3 -o tests/fixtures/deacon-0.19.1/exact-k33-w3.idx -q
//! ```
//! SHA-256, in the same order as the commands above:
//! 52f307574710e058fd8a33b694055d3248afd4c98eaf10bed78ad5af580af13a;
//! f23a67b86ec877167c09fb3b88f6c88a67d1b595634e83ba512b56bca89bda68;
//! 364c3d0b953808aa405731cd218fac8117ace82aac7a8e0866a0298caea35e51.
//! Reference SHA-256: d97f216bb28607385e44a9de97e16fd1829db57e5d88d3d52ad875bfb9e5399f.
//!
//! The literal examples were checked with `deacon filter --abs-threshold N
//! --rel-threshold 0 --ordered -t 1 -q --debug debug.tsv INDEX reads.fastq
//! -o output.fastq`, for N=1 and 2 on exact/fuse, and N=1 on wide exact.
//! Index creation and filtering occurred in a separate release executable; the
//! producer and nuclease still share Deacon's algorithm, so this is not an
//! independent algorithmic oracle.
//! Ideas inspired by <https://github.com/bede/deacon/blob/0.19.1/tests/index_tests.rs#L907-L955>
//! and <https://github.com/bede/deacon/blob/0.19.1/tests/filter_tests.rs#L1175-L1248>;
//! the sequences and test code here are original.

use std::{fs, path::PathBuf, process::Command};

const FORWARD: &str = "GCTACGATTCGGTACCAAGTCGATGCTAGGACCTTACGTCAGATCCGATTCAGGCTAACGTCCATGATCGTACGGAATCTGACCTTAGCGATCA";
const REVERSE: &str = "TGATCGCTAAGGTCAGATTCCGTACGATCATGGACGTTAGCCTGAATCGGATCTGACGTAAGGTCCTAGCATCGACTTGGTACCGAATCGTAGC";
const LOWER: &str = "gctacgattcggtaccaagtcgatgctaggaccttacgtcagatccgattcaggctaacgtccatgatcgtacggaatctgaccttagcgatca";
const REPEAT: &str = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
const MISS: &str = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC";
const AMBIGUOUS: &str = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
const MIXED: &str = "NNNNNGCTACGATTCGGTACCAAGTCGATGCTAGGNNNNN";
const MULTIPLE: &str = "AAAAAAAAAAAAAAAAAAAANGCTACGATTCGGTACCAAGTCGATGCTAGG";

fn check_selection(index: &str, threshold: usize, reads: &[(&str, &str, bool)]) {
    let dir = tempfile::tempdir().expect("temporary directory");
    let input_path = dir.path().join("reads.fastq");
    let mut input = String::new();
    let mut expected = String::new();
    for &(name, sequence, selected) in reads {
        let record = format!("@{name}\n{sequence}\n+\n{}\n", "I".repeat(sequence.len()));
        input.push_str(&record);
        if selected {
            expected.push_str(&record);
        }
    }
    fs::write(&input_path, input).expect("write reads");
    let index_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/deacon-0.19.1")
        .join(index);
    let output = Command::new(env!("CARGO_BIN_EXE_nuclease"))
        .arg("--in")
        .arg(input_path)
        .arg("--enrich-index")
        .arg(index_path)
        .arg("--enrich-min-hits")
        .arg(threshold.to_string())
        .args([
            "--min-length",
            "0",
            "--min-mean-q",
            "0",
            "--max-ns",
            "1000",
            "--min-entropy",
            "0",
            "--trim-min-q",
            "0",
            "-qqq",
        ])
        .output()
        .expect("run nuclease");
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        output.stdout,
        expected.as_bytes(),
        "index={index}, threshold={threshold}"
    );
}

#[test]
fn exact_producer_boundaries_strands_and_ambiguity() {
    check_selection(
        "exact-k7-w3.idx",
        1,
        &[
            ("short6", "AAAAAA", false),
            ("k7", "AAAAAAA", false),
            ("window_minus1_8", "AAAAAAAA", false),
            ("window9", "AAAAAAAAA", true),
            ("repeat", REPEAT, true),
            ("miss", MISS, false),
            ("forward", FORWARD, true),
            ("ambiguous", AMBIGUOUS, false),
            ("reverse", REVERSE, true),
            ("lower", LOWER, true),
            ("mixed", MIXED, true),
        ],
    );
}

#[test]
fn exact_and_fuse_have_the_observed_fixture_decisions() {
    // These fixed reads happen to agree, not a claim that fuse has no false positives.
    for index in ["exact-k7-w3.idx", "fuse32-k7-w3.pidx"] {
        check_selection(
            index,
            1,
            &[
                ("repeat", REPEAT, true),
                ("multiple", MULTIPLE, true),
                ("miss", MISS, false),
                ("mixed", MIXED, true),
            ],
        );
        check_selection(
            index,
            2,
            &[
                ("repeat", REPEAT, false),
                ("multiple", MULTIPLE, true),
                ("miss", MISS, false),
                ("mixed", MIXED, true),
            ],
        );
    }
}

#[test]
fn wide_exact_producer_supports_u128_keys() {
    check_selection(
        "exact-k33-w3.idx",
        1,
        &[
            ("wide32", "GCTACGATTCGGTACCAAGTCGATGCTAGGAC", false),
            ("wide33", "GCTACGATTCGGTACCAAGTCGATGCTAGGACC", false),
            ("wide34", "GCTACGATTCGGTACCAAGTCGATGCTAGGACCT", false),
            ("wide35", "GCTACGATTCGGTACCAAGTCGATGCTAGGACCTT", true),
            ("miss", MISS, false),
            ("forward", FORWARD, true),
            ("reverse", REVERSE, true),
            ("lower", LOWER, true),
            ("ambiguous", AMBIGUOUS, false),
        ],
    );
}
