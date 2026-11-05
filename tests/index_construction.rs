//! Reference construction, producer compatibility, and saved-index lifecycle regressions.

use std::{
    fs,
    io::Write,
    path::{Path, PathBuf},
    process::{Command, Output},
};

use flate2::{Compression, write::GzEncoder};
use tempfile::{TempDir, tempdir};

const HIT: &str = "GCTACGATTCGGTACCAAGTCGATGCTAGGACCTTACGTCAGATCCGATTCAGGCTAACGTCCATGATCGTACGGAATCTGACCTTAGCGATCA";
const MISS: &str = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC";

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/deacon-0.19.1")
        .join(name)
}

fn record(name: &str, sequence: &str) -> String {
    format!("@{name}\n{sequence}\n+\n{}\n", "I".repeat(sequence.len()))
}

fn setup() -> (TempDir, PathBuf) {
    let dir = tempdir().expect("temporary directory");
    let reads = dir.path().join("reads.fastq");
    fs::write(&reads, record("hit", HIT) + &record("miss", MISS)).expect("write reads");
    (dir, reads)
}

fn command(reads: &Path) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_nuclease"));
    command.arg("--in").arg(reads).args([
        "--min-length",
        "0",
        "--trim-min-q",
        "0",
        "--min-mean-q",
        "0",
        "--max-ns",
        "1000",
        "--min-entropy",
        "0",
    ]);
    command
}

fn parameters(command: &mut Command, k: u8, w: u8) -> &mut Command {
    command
        .arg("--index-kmer-length")
        .arg(k.to_string())
        .arg("--index-window-size")
        .arg(w.to_string())
}

fn run(command: &mut Command) -> Output {
    command.output().expect("run nuclease")
}

fn success(output: &Output) {
    assert!(
        output.status.success(),
        "status={:?}; stderr={}",
        output.status.code(),
        String::from_utf8_lossy(&output.stderr)
    );
}

fn failure(output: &Output) {
    assert!(!output.status.success(), "unexpected success: {output:?}");
    assert!(output.stdout.is_empty(), "failure emitted reads");
}

fn canonical_index(path: &Path) -> Vec<u8> {
    let index = deacon::load_index_from_path(path).expect("load saved exact index");
    let mut bytes = Vec::new();
    // Deacon 0.19.1 write_index sorts exact keys by (hash bucket, value),
    // not numeric order or hash-set iteration order. Normalize both sides.
    deacon::write_index(&index, &mut bytes).expect("canonical serialization");
    bytes
}

fn default_save(reference: &Path, k: u8, w: u8) -> PathBuf {
    let mut name = reference.as_os_str().to_os_string();
    name.push(format!(".k{k}.w{w}.idx"));
    PathBuf::from(name)
}

#[test]
fn constructed_exact_indexes_match_independent_producer_and_filtering_in_both_widths() {
    for k in [7, 33] {
        let (dir, reads) = setup();
        let saved = dir.path().join("constructed.idx");
        let expected_index = fixture(&format!("exact-k{k}-w3.idx"));
        let mut build = command(&reads);
        build
            .arg("--enrich-index")
            .arg(fixture("reference.fa"))
            .arg("--save-enrich-index")
            .arg(&saved);
        let output = run(parameters(&mut build, k, 3));
        success(&output);
        assert_eq!(output.stdout, record("hit", HIT).as_bytes());
        assert_eq!(canonical_index(&saved), canonical_index(&expected_index));
        let loaded = run(command(&reads).arg("--enrich-index").arg(expected_index));
        success(&loaded);
        assert_eq!(output.stdout, loaded.stdout);
    }
}

#[test]
fn fasta_and_fastq_references_are_detected_by_content_not_filename() {
    let (dir, reads) = setup();
    for (name, content) in [
        ("fasta.idx", format!(">reference\n{HIT}\n")),
        ("fastq.fa", record("reference", HIT)),
    ] {
        let reference = dir.path().join(name);
        fs::write(&reference, content).expect("reference");
        for (flag, expected) in [
            ("--enrich-index", record("hit", HIT)),
            ("--deplete-index", record("miss", MISS)),
        ] {
            let output = run(command(&reads).arg(flag).arg(&reference));
            success(&output);
            assert_eq!(output.stdout, expected.as_bytes(), "{name} {flag}");
        }
        assert!(!default_save(&reference, 31, 15).exists());
    }
}

#[test]
fn gzip_reference_without_compressed_extension_builds_the_same_index() {
    let (dir, reads) = setup();
    let reference = dir.path().join("compressed.idx");
    let saved = dir.path().join("saved.idx");
    let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
    encoder
        .write_all(&fs::read(fixture("reference.fa")).expect("fixture"))
        .expect("compress reference");
    fs::write(&reference, encoder.finish().expect("finish gzip")).expect("reference");
    let mut build = command(&reads);
    build
        .arg("--enrich-index")
        .arg(&reference)
        .arg("--save-enrich-index")
        .arg(&saved);
    let output = run(parameters(&mut build, 7, 3));
    success(&output);
    assert_eq!(output.stdout, record("hit", HIT).as_bytes());
    assert_eq!(
        canonical_index(&saved),
        canonical_index(&fixture("exact-k7-w3.idx"))
    );
}

#[test]
fn bare_save_appends_parameters_to_the_full_reference_filename() {
    for (source_flag, save_flag) in [
        ("--enrich-index", "--save-enrich-index"),
        ("--deplete-index", "--save-deplete-index"),
    ] {
        let (dir, reads) = setup();
        let reference = dir.path().join("reference.fa.gz");
        fs::copy(fixture("reference.fa"), &reference).expect("reference");
        let output = run(command(&reads)
            .arg(source_flag)
            .arg(&reference)
            .arg(save_flag));
        success(&output);
        let saved = default_save(&reference, 31, 15);
        assert!(saved.is_file());
        let loaded = run(command(&reads).arg(source_flag).arg(saved));
        success(&loaded);
        assert_eq!(output.stdout, loaded.stdout);
    }
}

#[test]
fn bare_save_uses_explicit_parameters_and_explicit_save_does_not_create_default() {
    let (dir, reads) = setup();
    let reference = dir.path().join("reference.fa");
    fs::copy(fixture("reference.fa"), &reference).expect("reference");
    let mut build = command(&reads);
    build
        .arg("--enrich-index")
        .arg(&reference)
        .arg("--save-enrich-index");
    success(&run(parameters(&mut build, 7, 3)));
    assert_eq!(
        canonical_index(&default_save(&reference, 7, 3)),
        canonical_index(&fixture("exact-k7-w3.idx"))
    );
    let other = dir.path().join("other.fa");
    fs::copy(&reference, &other).expect("second reference");
    let saved = dir.path().join("chosen.idx");
    success(&run(command(&reads)
        .arg("--deplete-index")
        .arg(&other)
        .arg("--save-deplete-index")
        .arg(&saved)));
    assert!(saved.is_file());
    assert!(!default_save(&other, 31, 15).exists());
}

#[test]
fn save_requires_a_reference_not_a_prebuilt_index_or_an_absent_role() {
    let (dir, reads) = setup();
    for (source_flag, save_flag) in [
        ("--enrich-index", "--save-enrich-index"),
        ("--deplete-index", "--save-deplete-index"),
    ] {
        for source in [None, Some(fixture("exact-k7-w3.idx"))] {
            let saved = dir.path().join("must-not-exist.idx");
            let mut cli = command(&reads);
            if let Some(source) = source {
                cli.arg(source_flag).arg(source);
            }
            failure(&run(cli.arg(save_flag).arg(&saved)));
            assert!(!saved.exists());
        }
    }
}

#[test]
fn existing_save_target_is_never_overwritten() {
    let (dir, reads) = setup();
    let saved = dir.path().join("existing.idx");
    fs::write(&saved, b"keep me").expect("existing target");
    failure(&run(command(&reads)
        .arg("--enrich-index")
        .arg(fixture("reference.fa"))
        .arg("--save-enrich-index")
        .arg(&saved)));
    assert_eq!(fs::read(saved).expect("preserved target"), b"keep me");
}

#[test]
fn both_save_plans_are_validated_before_either_reference_is_built() {
    let (dir, reads) = setup();
    let saved = dir.path().join("first.idx");
    let existing = dir.path().join("occupied.idx");
    let second = dir.path().join("other.fa");
    fs::copy(fixture("reference.fa"), &second).expect("second reference");
    fs::write(&existing, b"preserved").expect("occupied target");
    failure(&run(command(&reads)
        .arg("--enrich-index")
        .arg(fixture("reference.fa"))
        .arg("--save-enrich-index")
        .arg(&saved)
        .arg("--deplete-index")
        .arg(second)
        .arg("--save-deplete-index")
        .arg(&existing)));
    assert!(
        !saved.exists(),
        "first index built before second save validation"
    );
    assert_eq!(fs::read(existing).expect("target"), b"preserved");
}

#[test]
fn missing_second_source_is_validated_before_first_save() {
    let (dir, reads) = setup();
    let saved = dir.path().join("first.idx");
    failure(&run(command(&reads)
        .arg("--enrich-index")
        .arg(fixture("reference.fa"))
        .arg("--save-enrich-index")
        .arg(&saved)
        .arg("--deplete-index")
        .arg(dir.path().join("missing.fa"))));
    assert!(!saved.exists());
}

#[test]
fn save_targets_cannot_collide_with_each_other_or_run_outputs() {
    let (dir, reads) = setup();
    for output_flag in [
        "--out",
        "--out1",
        "--out2",
        "--summary",
        "--invalid-input-report",
    ] {
        let target = dir.path().join("collision");
        let other_output = dir.path().join("other-output.fastq");
        let mut cli = command(&reads);
        cli.arg("--enrich-index")
            .arg(fixture("reference.fa"))
            .arg("--save-enrich-index")
            .arg(&target)
            .arg(output_flag)
            .arg(&target);
        // Supply a valid split-output plan: a missing mate option would make
        // this pass at argument parsing without exercising collision checks.
        if output_flag == "--out1" || output_flag == "--out2" {
            let mate_flag = if output_flag == "--out1" {
                "--out2"
            } else {
                "--out1"
            };
            cli.arg("--paired").arg(mate_flag).arg(&other_output);
        }
        failure(&run(&mut cli));
        assert!(
            !target.exists(),
            "created collision target for {output_flag}"
        );
        assert!(!other_output.exists());
    }
    let reference = dir.path().join("other.fa");
    fs::copy(fixture("reference.fa"), &reference).expect("reference");
    let target = dir.path().join("shared.idx");
    failure(&run(command(&reads)
        .arg("--enrich-index")
        .arg(fixture("reference.fa"))
        .arg("--save-enrich-index")
        .arg(&target)
        .arg("--deplete-index")
        .arg(reference)
        .arg("--save-deplete-index")
        .arg(&target)));
    assert!(!target.exists());
}

#[test]
fn identical_source_paths_are_rejected_but_distinct_files_with_identical_content_are_allowed() {
    let (dir, reads) = setup();
    let reference = fixture("reference.fa");
    failure(&run(command(&reads)
        .arg("--enrich-index")
        .arg(&reference)
        .arg("--deplete-index")
        .arg(&reference)));
    let separate = dir.path().join("separate.fa");
    fs::copy(&reference, &separate).expect("independent copy");
    success(&run(command(&reads)
        .arg("--enrich-index")
        .arg(reference)
        .arg("--deplete-index")
        .arg(separate)));
}

#[cfg(unix)]
#[test]
fn symlink_and_hardlink_source_aliases_are_rejected_before_saving() {
    let (dir, reads) = setup();
    let reference = dir.path().join("reference.fa");
    fs::copy(fixture("reference.fa"), &reference).expect("reference");
    let symbolic = dir.path().join("symbolic.fa");
    let hard = dir.path().join("hard.fa");
    std::os::unix::fs::symlink(&reference, &symbolic).expect("symlink");
    fs::hard_link(&reference, &hard).expect("hard link");
    for alias in [symbolic, hard] {
        let saved = dir.path().join("must-not-exist.idx");
        failure(&run(command(&reads)
            .arg("--enrich-index")
            .arg(&reference)
            .arg("--save-enrich-index")
            .arg(&saved)
            .arg("--deplete-index")
            .arg(alias)));
        assert!(!saved.exists());
    }
}

#[cfg(unix)]
#[test]
fn dangling_save_symlink_is_not_followed_or_replaced() {
    let (dir, reads) = setup();
    let target = dir.path().join("missing-target.idx");
    let link = dir.path().join("save.idx");
    std::os::unix::fs::symlink(&target, &link).expect("dangling symlink");
    failure(&run(command(&reads)
        .arg("--enrich-index")
        .arg(fixture("reference.fa"))
        .arg("--save-enrich-index")
        .arg(&link)));
    assert_eq!(fs::read_link(link).expect("symlink preserved"), target);
    assert!(!target.exists());
}

#[test]
fn saved_enrichment_survives_later_depletion_construction_failure() {
    let (dir, reads) = setup();
    let empty = dir.path().join("too-short.fa");
    let saved = dir.path().join("enrich.idx");
    let failed_save = dir.path().join("deplete.idx");
    fs::write(&empty, b">short\nACG\n").expect("short reference");
    let mut cli = command(&reads);
    cli.arg("--enrich-index")
        .arg(fixture("reference.fa"))
        .arg("--save-enrich-index")
        .arg(&saved)
        .arg("--deplete-index")
        .arg(empty)
        .arg("--save-deplete-index")
        .arg(&failed_save);
    failure(&run(parameters(&mut cli, 7, 3)));
    assert_eq!(
        canonical_index(&saved),
        canonical_index(&fixture("exact-k7-w3.idx"))
    );
    assert!(!failed_save.exists());
}

#[test]
fn saved_index_survives_later_read_ingress_failure_and_reads_remain_fastq_only() {
    let (dir, _) = setup();
    let reads = dir.path().join("reads.fa");
    let saved = dir.path().join("saved.idx");
    fs::write(&reads, format!(">read\n{HIT}\n")).expect("FASTA reads");
    let mut cli = command(&reads);
    cli.arg("--enrich-index")
        .arg(fixture("reference.fa"))
        .arg("--save-enrich-index")
        .arg(&saved);
    failure(&run(parameters(&mut cli, 7, 3)));
    assert_eq!(
        canonical_index(&saved),
        canonical_index(&fixture("exact-k7-w3.idx"))
    );
}

#[test]
fn no_minimizers_reports_reference_path_and_parameters_without_saving() {
    let (dir, reads) = setup();
    for content in [">short\nACG\n", ">ambiguous\nNNNNNNNNNNNNNNNNNNNN\n"] {
        let reference = dir.path().join("empty-minimizers.fa");
        let saved = dir.path().join("never.idx");
        fs::write(&reference, content).expect("reference");
        let mut cli = command(&reads);
        cli.arg("--enrich-index")
            .arg(&reference)
            .arg("--save-enrich-index")
            .arg(&saved);
        let output = run(parameters(&mut cli, 7, 3));
        failure(&output);
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(stderr.contains("empty-minimizers.fa"), "{stderr}");
        assert!(stderr.contains('7') && stderr.contains('3'), "{stderr}");
        assert!(!saved.exists());
    }
}

#[test]
fn malformed_fastq_reference_is_not_retried_as_a_prebuilt_index() {
    let (dir, reads) = setup();
    let reference = dir.path().join("malformed.idx");
    let saved = dir.path().join("never.idx");
    fs::write(&reference, format!("@broken\n{HIT}\n+\nI\n")).expect("reference");
    let output = run(command(&reads)
        .arg("--enrich-index")
        .arg(&reference)
        .arg("--save-enrich-index")
        .arg(&saved));
    failure(&output);
    let stderr = String::from_utf8_lossy(&output.stderr).to_lowercase();
    assert!(stderr.contains("reference"), "{stderr}");
    assert!(!saved.exists());
}

#[test]
fn explicit_unused_parameters_warn_without_overriding_prebuilt_parameters() {
    let (_dir, reads) = setup();
    fs::write(&reads, record("short-hit", "AAAAAAAAA")).expect("boundary read");
    let index = fixture("exact-k7-w3.idx");
    let implicit = run(command(&reads)
        .arg("--enrich-index")
        .arg(&index)
        .args(["--enrich-min-hits", "1"]));
    success(&implicit);
    assert_eq!(implicit.stdout, record("short-hit", "AAAAAAAAA").as_bytes());
    let stderr = String::from_utf8_lossy(&implicit.stderr).to_lowercase();
    assert!(
        !stderr.contains("unused") && !stderr.contains("ignored"),
        "{stderr}"
    );
    for args in [
        vec!["--index-kmer-length", "33"],
        vec!["--index-window-size", "15"],
        vec!["--index-kmer-length", "33", "--index-window-size", "15"],
    ] {
        let output = run(command(&reads)
            .arg("--enrich-index")
            .arg(&index)
            .args(["--enrich-min-hits", "1"])
            .args(args));
        success(&output);
        assert_eq!(output.stdout, implicit.stdout);
        let stderr = String::from_utf8_lossy(&output.stderr).to_lowercase();
        assert!(
            stderr.contains("warn"),
            "missing unused-setting warning: {stderr}"
        );
    }
}

#[test]
fn prebuilt_indexes_keep_their_own_parameters_when_mixing_widths_and_reference_inputs() {
    let (dir, reads) = setup();
    for (enrich, deplete) in [
        (fixture("exact-k7-w3.idx"), fixture("exact-k33-w3.idx")),
        (fixture("exact-k33-w3.idx"), fixture("exact-k7-w3.idx")),
    ] {
        let output = run(command(&reads)
            .arg("--enrich-index")
            .arg(enrich)
            .arg("--deplete-index")
            .arg(deplete)
            .args(["--selection-policy", "depletion-veto"]));
        success(&output);
        assert!(
            output.stdout.is_empty(),
            "both fixture indexes must recognize the hit"
        );
    }
    fs::write(
        &reads,
        record("short-hit", "AAAAAAAAA") + &record("miss", MISS),
    )
    .expect("boundary reads");
    let other = dir.path().join("miss.fa");
    fs::write(&other, format!(">miss\n{MISS}\n")).expect("reference");
    for k in [7, 33] {
        let mut cli = command(&reads);
        cli.arg("--enrich-index")
            .arg(fixture("exact-k7-w3.idx"))
            .args(["--enrich-min-hits", "1"])
            .arg("--deplete-index")
            .arg(&other)
            .args(["--selection-policy", "depletion-veto"]);
        let output = run(parameters(&mut cli, k, 3));
        success(&output);
        assert_eq!(output.stdout, record("short-hit", "AAAAAAAAA").as_bytes());
    }
}

#[test]
fn prebuilt_index_is_detected_by_content_even_with_a_reference_extension() {
    let (dir, reads) = setup();
    let renamed = dir.path().join("prebuilt.fa");
    fs::copy(fixture("exact-k7-w3.idx"), &renamed).expect("renamed prebuilt");
    let output = run(command(&reads).arg("--enrich-index").arg(&renamed));
    success(&output);
    assert_eq!(output.stdout, record("hit", HIT).as_bytes());
    assert!(!default_save(&renamed, 31, 15).exists());
}

#[test]
fn shared_index_settings_reject_values_outside_u8_at_argument_parsing() {
    let (_dir, reads) = setup();
    for setting in ["--index-kmer-length", "--index-window-size"] {
        for invalid in ["256", "-1", "not-a-number"] {
            let output = run(command(&reads)
                .arg("--enrich-index")
                .arg(fixture("reference.fa"))
                .arg(format!("{setting}={invalid}")));
            assert_eq!(
                output.status.code(),
                Some(2),
                "{setting}={invalid}: {output:?}"
            );
            assert!(output.stdout.is_empty());
        }
    }
}

#[test]
fn other_compression_formats_match_plain_reference_and_reject_truncation() {
    // `printf '>ref\nACGTACGTTGCAACGT\n' | {bzip2,xz,zstd} -c`.
    // Embedded bytes keep the tests independent of system codec executables.
    let compressed: &[&[u8]] = &[
        &[
            0x42, 0x5a, 0x68, 0x39, 0x31, 0x41, 0x59, 0x26, 0x53, 0x59, 0x18, 0xfa, 0xec, 0x50,
            0x00, 0x00, 0x01, 0x4f, 0x80, 0x00, 0x10, 0x00, 0x01, 0x28, 0x80, 0x04, 0x00, 0x03,
            0x00, 0x10, 0x00, 0x20, 0x00, 0x22, 0x0c, 0x23, 0x42, 0x0c, 0x98, 0x87, 0xc1, 0xb0,
            0xb4, 0xd8, 0x56, 0x68, 0x64, 0x93, 0xc5, 0xdc, 0x91, 0x4e, 0x14, 0x24, 0x06, 0x3e,
            0xbb, 0x14, 0x00,
        ],
        &[
            0xfd, 0x37, 0x7a, 0x58, 0x5a, 0x00, 0x00, 0x04, 0xe6, 0xd6, 0xb4, 0x46, 0x04, 0xc0,
            0x1a, 0x16, 0x21, 0x01, 0x16, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
            0x1c, 0xaa, 0x3b, 0x74, 0x01, 0x00, 0x15, 0x3e, 0x72, 0x65, 0x66, 0x0a, 0x41, 0x43,
            0x47, 0x54, 0x41, 0x43, 0x47, 0x54, 0x54, 0x47, 0x43, 0x41, 0x41, 0x43, 0x47, 0x54,
            0x0a, 0x00, 0x00, 0x00, 0x61, 0xea, 0x11, 0xe9, 0x49, 0x7a, 0x84, 0x3f, 0x00, 0x01,
            0x36, 0x16, 0x0f, 0x91, 0x4e, 0x5d, 0x1f, 0xb6, 0xf3, 0x7d, 0x01, 0x00, 0x00, 0x00,
            0x00, 0x04, 0x59, 0x5a,
        ],
        &[
            0x28, 0xb5, 0x2f, 0xfd, 0x04, 0x58, 0xb1, 0x00, 0x00, 0x3e, 0x72, 0x65, 0x66, 0x0a,
            0x41, 0x43, 0x47, 0x54, 0x41, 0x43, 0x47, 0x54, 0x54, 0x47, 0x43, 0x41, 0x41, 0x43,
            0x47, 0x54, 0x0a, 0xc3, 0x15, 0x80, 0x5e,
        ],
    ];
    let (dir, reads) = setup();
    let reference = dir.path().join("reference.idx");
    let plain = dir.path().join("plain.idx");
    fs::write(&reference, b">ref\nACGTACGTTGCAACGT\n").expect("reference");
    success(&run(parameters(
        command(&reads)
            .arg("--enrich-index")
            .arg(&reference)
            .arg("--save-enrich-index")
            .arg(&plain),
        3,
        1,
    )));
    for bytes in compressed {
        let saved = dir.path().join("compressed.idx");
        fs::write(&reference, bytes).expect("compressed reference");
        success(&run(parameters(
            command(&reads)
                .arg("--enrich-index")
                .arg(&reference)
                .arg("--save-enrich-index")
                .arg(&saved),
            3,
            1,
        )));
        assert_eq!(canonical_index(&saved), canonical_index(&plain));
        fs::remove_file(&saved).expect("remove saved index");
        fs::write(&reference, &bytes[..bytes.len() / 2]).expect("truncated reference");
        let output = run(parameters(
            command(&reads)
                .arg("--enrich-index")
                .arg(&reference)
                .arg("--save-enrich-index")
                .arg(&saved),
            3,
            1,
        ));
        assert_eq!(output.status.code(), Some(65), "{output:?}");
        assert!(!saved.exists());
        assert!(!fs::read_dir(dir.path()).expect("entries").any(|entry| {
            entry
                .expect("entry")
                .file_name()
                .to_string_lossy()
                .starts_with(".nuclease-")
        }));
    }
}

#[cfg(unix)]
#[test]
fn dangling_output_links_and_symlinked_output_parents_collide_before_saving() {
    let (dir, reads) = setup();
    let destination = dir.path().join("saved.idx");
    let alias = dir.path().join("output-link");
    std::os::unix::fs::symlink("saved.idx", &alias).expect("dangling output link");
    let parent_alias = dir.path().join("parent-link");
    std::os::unix::fs::symlink(dir.path(), &parent_alias).expect("parent link");
    for output in [alias, parent_alias.join("saved.idx")] {
        failure(&run(command(&reads)
            .arg("--enrich-index")
            .arg(fixture("reference.fa"))
            .arg("--save-enrich-index")
            .arg(&destination)
            .arg("--out")
            .arg(output)));
        assert!(!destination.exists());
    }
}

#[test]
fn destination_spelling_uses_filesystem_rules_before_construction() {
    let (dir, reads) = setup();
    let probe = dir.path().join("case-check");
    fs::write(&probe, b"").expect("case probe");
    let insensitive = dir.path().join("CASE-CHECK").exists();
    fs::remove_file(probe).expect("remove probe");
    let reference = dir.path().join("other.fa");
    fs::copy(fixture("reference.fa"), &reference).expect("second reference");
    for dual_save in [false, true] {
        let destination = dir.path().join("saved.idx");
        let alias = dir.path().join("SAVED.idx");
        let mut cli = command(&reads);
        cli.arg("--enrich-index")
            .arg(fixture("reference.fa"))
            .arg("--save-enrich-index")
            .arg(&destination);
        if dual_save {
            cli.arg("--deplete-index")
                .arg(&reference)
                .arg("--save-deplete-index")
                .arg(&alias);
        } else {
            cli.arg("--out").arg(&alias);
        }
        let output = run(&mut cli);
        if insensitive {
            failure(&output);
            assert!(
                !destination.exists(),
                "collision was detected after publication"
            );
        } else {
            success(&output);
            fs::remove_file(&destination).expect("saved index");
            fs::remove_file(&alias).expect("distinct output");
        }
    }
}

#[test]
fn invalid_deacon_parameters_fail_before_saving_or_opening_outputs() {
    let (dir, reads) = setup();
    let saved = dir.path().join("saved.idx");
    let output = dir.path().join("reads.out");
    for (k, w) in [(0, 1), (1, 0), (2, 1), (3, 2), (63, 1), (61, 37)] {
        let result = run(parameters(
            command(&reads)
                .arg("--enrich-index")
                .arg(fixture("reference.fa"))
                .arg("--save-enrich-index")
                .arg(&saved)
                .arg("--out")
                .arg(&output),
            k,
            w,
        ));
        assert_eq!(result.status.code(), Some(2), "{result:?}");
        assert!(!saved.exists());
        assert!(!output.exists());
    }
}
