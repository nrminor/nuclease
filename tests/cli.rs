//! CLI contract tests for `nuclease`.

use std::{fs, process::Command};

use serde_json::Value;
use tempfile::tempdir;

fn nuclease() -> Command {
    Command::new(env!("CARGO_BIN_EXE_nuclease"))
}

#[test]
fn local_single_fastq_streams_cleaned_reads_and_writes_summary() {
    let temp = tempdir().expect("tempdir should be created");
    let input = temp.path().join("reads.fastq");
    let summary = temp.path().join("summary.json");
    fs::write(&input, b"@read1\nACGT\n+\nIIII\n@read2\nTGCA\n+\nJJJJ\n")
        .expect("fixture FASTQ should be writable");

    let output = nuclease()
        .args([
            "--in1",
            input.to_str().expect("fixture path should be UTF-8"),
            "--min-length",
            "4",
            "--trim-min-q",
            "0",
            "--min-mean-q",
            "0",
            "--summary",
            summary.to_str().expect("summary path should be UTF-8"),
            "-qqq",
        ])
        .output()
        .expect("nuclease should run");

    assert!(
        output.status.success(),
        "nuclease failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        output.stdout,
        b"@read1\nACGT\n+\nIIII\n@read2\nTGCA\n+\nJJJJ\n"
    );

    let summary_json = fs::read_to_string(summary).expect("summary should be readable");
    let summary: Value = serde_json::from_str(&summary_json).expect("summary should be JSON");
    assert_eq!(summary["reads_seen"], 2);
    assert_eq!(summary["reads_emitted"], 2);
    assert_eq!(summary["reads_rejected"], 0);
    assert_eq!(summary["invalid_reads"], 0);
}

#[test]
fn adapter_preset_defaults_to_illumina_truseq() {
    let temp = tempdir().expect("tempdir should be created");
    let input = temp.path().join("reads.fastq");
    fs::write(&input, b"@read1\nACGTAGATCGGAAG\n+\nIIIIIIIIIIIIII\n")
        .expect("fixture FASTQ should be writable");

    let output = nuclease()
        .args([
            "--in1",
            input.to_str().expect("fixture path should be UTF-8"),
            "--min-length",
            "1",
            "--trim-min-q",
            "0",
            "--min-mean-q",
            "0",
            "-qqq",
        ])
        .output()
        .expect("nuclease should run");

    assert!(
        output.status.success(),
        "nuclease failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(output.stdout, b"@read1\nACGT\n+\nIIII\n");
}

#[test]
fn adapter_preset_none_skips_adapter_trimming() {
    let temp = tempdir().expect("tempdir should be created");
    let input = temp.path().join("reads.fastq");
    let fixture = b"@read1\nACGTAGATCGGAAG\n+\nIIIIIIIIIIIIII\n";
    fs::write(&input, fixture).expect("fixture FASTQ should be writable");

    let output = nuclease()
        .args([
            "--in1",
            input.to_str().expect("fixture path should be UTF-8"),
            "--adapter-preset",
            "none",
            "--min-length",
            "1",
            "--trim-min-q",
            "0",
            "--min-mean-q",
            "0",
            "-qqq",
        ])
        .output()
        .expect("nuclease should run");

    assert!(
        output.status.success(),
        "nuclease failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(output.stdout, fixture);
}

#[test]
fn warn_drop_invalid_fastq_policy_does_not_recover_parser_error() {
    let temp = tempdir().expect("tempdir should be created");
    let input = temp.path().join("reads.fastq");
    fs::write(
        &input,
        b"@good1\nACGT\n+\nIIII\n@bad\nAAAA\n+\nI\n@good2\nTGCA\n+\nJJJJ\n",
    )
    .expect("fixture FASTQ should be writable");

    let output = nuclease()
        .args([
            "--in1",
            input.to_str().expect("fixture path should be UTF-8"),
            "--min-length",
            "4",
            "--trim-min-q",
            "0",
            "--min-mean-q",
            "0",
            "--invalid-fastq-policy",
            "warn-drop",
            "-qqq",
        ])
        .output()
        .expect("nuclease should run");

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        !output.status.success(),
        "parser-level FASTQ errors should remain fatal under warn-drop"
    );
    assert!(
        stderr.contains("FASTQ parser rejected malformed input"),
        "stderr did not include controlled parser diagnostic: {stderr}"
    );
    assert!(
        stderr.contains("invalid_fastq_policy=warn_drop"),
        "stderr did not include active invalid FASTQ policy: {stderr}"
    );
    assert!(
        stderr.contains("parser_error_kind=UnequalLengths"),
        "stderr did not include needletail error kind: {stderr}"
    );
    assert!(
        !stderr.contains("The application panicked"),
        "parser-level error should not surface as a panic: {stderr}"
    );
}

#[test]
fn invalid_fastq_report_writes_fatal_parser_error_as_jsonl() {
    let temp = tempdir().expect("tempdir should be created");
    let input = temp.path().join("reads.fastq");
    let report = temp.path().join("invalid-fastq.jsonl");
    fs::write(&input, b"@bad1\nAAAA\n+\nI\n@good\nTGCA\n+\nJJJJ\n")
        .expect("fixture FASTQ should be writable");

    let output = nuclease()
        .args([
            "--in1",
            input.to_str().expect("fixture path should be UTF-8"),
            "--min-length",
            "4",
            "--trim-min-q",
            "0",
            "--min-mean-q",
            "0",
            "--invalid-fastq-policy",
            "silent-drop",
            "--invalid-fastq-report",
            report.to_str().expect("report path should be UTF-8"),
            "-qqq",
        ])
        .output()
        .expect("nuclease should run");

    assert!(
        !output.status.success(),
        "parser-level FASTQ errors should be fatal: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let report = fs::read_to_string(report).expect("invalid FASTQ report should be readable");
    let events = report
        .lines()
        .map(|line| serde_json::from_str::<Value>(line).expect("event should be JSON"))
        .collect::<Vec<_>>();
    assert_eq!(events.len(), 1);
    assert_eq!(events[0]["kind"], "fastq_parse_error");
    assert_eq!(events[0]["mate"], "single");
    assert_eq!(events[0]["policy"], "silent_drop");
    assert_eq!(events[0]["recoverable"], false);
    assert_eq!(events[0]["fatal"], true);
    assert_eq!(events[0]["parser_error_kind"], "UnequalLengths");
    assert_eq!(events[0]["parser_error_line"], 1);
    assert!(
        events[0]["parser_error_message"]
            .as_str()
            .is_some_and(|message| message.contains("quality length is 1")),
        "event did not preserve parser error message: {:?}",
        events[0]
    );
}

#[test]
fn malformed_fastq_does_not_surface_raw_parser_slice_panic() {
    let temp = tempdir().expect("tempdir should be created");
    let input = temp.path().join("malformed.fastq");
    fs::write(
        &input,
        concat!(
            "@padding\n",
            "A\n",
            "+\n",
            "I\n",
            "@bad\n",
            "ACGT\n",
            "+\n",
            "\n",
            "\n",
            "ACGTACGT\n",
            "+\n",
            "!!!!!!!!\n",
            "@after\n",
            "ACGT\n",
            "+\n",
            "IIII\n",
        ),
    )
    .expect("malformed FASTQ fixture should be writable");

    let output = nuclease()
        .args([
            "--in1",
            input.to_str().expect("fixture path should be UTF-8"),
            "--min-length",
            "1",
            "--trim-min-q",
            "0",
            "--min-mean-q",
            "0",
            "--invalid-fastq-policy",
            "warn-drop",
            "-qqq",
        ])
        .output()
        .expect("nuclease should run");

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        !stderr.contains("slice index starts"),
        "raw parser slice panic leaked to stderr: {stderr}"
    );
    assert!(
        !stderr.contains("The application panicked"),
        "raw parser panic banner leaked to stderr: {stderr}"
    );
    assert!(
        output.status.success() || stderr.contains("FASTQ parser rejected malformed input"),
        "malformed input should either warn/drop successfully or fail with controlled parser diagnostics: {stderr}"
    );
}

#[test]
fn paired_fastq_streams_interleaved_reads_and_writes_summary() {
    let temp = tempdir().expect("tempdir should be created");
    let input1 = temp.path().join("reads_1.fastq");
    let input2 = temp.path().join("reads_2.fastq");
    let summary = temp.path().join("summary.json");
    fs::write(
        &input1,
        b"@read1/1\nAAAA\n+\nIIII\n@read2/1\nCCCC\n+\nJJJJ\n",
    )
    .expect("read 1 fixture should be writable");
    fs::write(
        &input2,
        b"@read1/2\nTTTT\n+\nKKKK\n@read2/2\nGGGG\n+\nLLLL\n",
    )
    .expect("read 2 fixture should be writable");

    let output = nuclease()
        .args([
            "--in1",
            input1.to_str().expect("read 1 path should be UTF-8"),
            "--in2",
            input2.to_str().expect("read 2 path should be UTF-8"),
            "--interleaved",
            "--min-length",
            "4",
            "--trim-min-q",
            "0",
            "--min-mean-q",
            "0",
            "--summary",
            summary.to_str().expect("summary path should be UTF-8"),
            "-qqq",
        ])
        .output()
        .expect("nuclease should run");

    assert!(
        output.status.success(),
        "nuclease failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        output.stdout,
        b"@read1/1\nAAAA\n+\nIIII\n@read1/2\nTTTT\n+\nKKKK\n@read2/1\nCCCC\n+\nJJJJ\n@read2/2\nGGGG\n+\nLLLL\n"
    );

    let summary_json = fs::read_to_string(summary).expect("summary should be readable");
    let summary: Value = serde_json::from_str(&summary_json).expect("summary should be JSON");
    assert_eq!(summary["reads_seen"], 4);
    assert_eq!(summary["reads_emitted"], 4);
    assert_eq!(summary["pairs_seen"], 2);
    assert_eq!(summary["pairs_emitted"], 2);
}

#[test]
fn paired_fastq_merge_pairs_emits_merged_record() {
    let temp = tempdir().expect("tempdir should be created");
    let input1 = temp.path().join("reads_1.fastq");
    let input2 = temp.path().join("reads_2.fastq");
    let summary = temp.path().join("summary.json");
    fs::write(
        &input1,
        b"@read-1/1\nACGTTGCAGTACGATCGTACGGAATTCGCCGATGACTGACCTAGGTCAGTACGATC\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n",
    )
    .expect("read 1 fixture should be writable");
    fs::write(
        &input2,
        b"@read-1/2\nGATCGTACTGACCTAGGTCAGTCATCGGCGAATTCCGTACGATCGTACTGCAACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n",
    )
    .expect("read 2 fixture should be writable");

    let output = nuclease()
        .args([
            "--in1",
            input1.to_str().expect("read 1 path should be UTF-8"),
            "--in2",
            input2.to_str().expect("read 2 path should be UTF-8"),
            "--merge-pairs",
            "--interleaved",
            "--adapter-preset",
            "none",
            "--min-length",
            "1",
            "--trim-min-q",
            "0",
            "--min-mean-q",
            "0",
            "--summary",
            summary.to_str().expect("summary path should be UTF-8"),
            "-qqq",
        ])
        .output()
        .expect("nuclease should run");

    assert!(
        output.status.success(),
        "nuclease failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stdout = String::from_utf8(output.stdout).expect("FASTQ output should be UTF-8");
    assert!(
        stdout.starts_with("@read-1\n"),
        "merged read should use normalized pair id as header: {stdout}"
    );
    assert_eq!(
        stdout.matches('\n').count(),
        4,
        "merged pair should emit one FASTQ record: {stdout}"
    );

    let summary_json = fs::read_to_string(summary).expect("summary should be readable");
    let summary: Value = serde_json::from_str(&summary_json).expect("summary should be JSON");
    assert_eq!(summary["reads_seen"], 2);
    assert_eq!(summary["reads_emitted"], 1);
    assert_eq!(summary["pairs_seen"], 1);
    assert_eq!(summary["pairs_emitted"], 1);
}

#[test]
fn paired_fastq_merge_pairs_keeps_unmerged_pair() {
    let temp = tempdir().expect("tempdir should be created");
    let input1 = temp.path().join("reads_1.fastq");
    let input2 = temp.path().join("reads_2.fastq");
    let summary = temp.path().join("summary.json");
    fs::write(&input1, b"@read1/1\nAAAAAAAAAAAA\n+\nIIIIIIIIIIII\n")
        .expect("read 1 fixture should be writable");
    fs::write(&input2, b"@read1/2\nCCCCCCCCCCCC\n+\nJJJJJJJJJJJJ\n")
        .expect("read 2 fixture should be writable");

    let output = nuclease()
        .args([
            "--in1",
            input1.to_str().expect("read 1 path should be UTF-8"),
            "--in2",
            input2.to_str().expect("read 2 path should be UTF-8"),
            "--merge-pairs",
            "--interleaved",
            "--adapter-preset",
            "none",
            "--min-length",
            "1",
            "--trim-min-q",
            "0",
            "--min-mean-q",
            "0",
            "--summary",
            summary.to_str().expect("summary path should be UTF-8"),
            "-qqq",
        ])
        .output()
        .expect("nuclease should run");

    assert!(
        output.status.success(),
        "nuclease failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        output.stdout,
        b"@read1/1\nAAAAAAAAAAAA\n+\nIIIIIIIIIIII\n@read1/2\nCCCCCCCCCCCC\n+\nJJJJJJJJJJJJ\n"
    );

    let summary_json = fs::read_to_string(summary).expect("summary should be readable");
    let summary: Value = serde_json::from_str(&summary_json).expect("summary should be JSON");
    assert_eq!(summary["reads_seen"], 2);
    assert_eq!(summary["reads_emitted"], 2);
    assert_eq!(summary["pairs_seen"], 1);
    assert_eq!(summary["pairs_emitted"], 1);
}

#[test]
fn paired_fastq_merge_min_overlap_can_reject_shorter_overlap() {
    let temp = tempdir().expect("tempdir should be created");
    let input1 = temp.path().join("reads_1.fastq");
    let input2 = temp.path().join("reads_2.fastq");
    fs::write(
        &input1,
        b"@read-1/1\nACGTTGCAGTACGATCGTACGGAATTCGCCGATGACTGACCTAGGTCAGTACGATC\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n",
    )
    .expect("read 1 fixture should be writable");
    fs::write(
        &input2,
        b"@read-1/2\nGATCGTACTGACCTAGGTCAGTCATCGGCGAATTCCGTACGATCGTACTGCAACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n",
    )
    .expect("read 2 fixture should be writable");

    let output = nuclease()
        .args([
            "--in1",
            input1.to_str().expect("read 1 path should be UTF-8"),
            "--in2",
            input2.to_str().expect("read 2 path should be UTF-8"),
            "--merge-pairs",
            "--merge-min-overlap",
            "80",
            "--interleaved",
            "--adapter-preset",
            "none",
            "--min-length",
            "1",
            "--trim-min-q",
            "0",
            "--min-mean-q",
            "0",
            "-qqq",
        ])
        .output()
        .expect("nuclease should run");

    assert!(
        output.status.success(),
        "nuclease failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        bytecount::count(&output.stdout, b'@'),
        2,
        "high merge-min-overlap should preserve the original pair: {}",
        String::from_utf8_lossy(&output.stdout)
    );
}

#[test]
fn merge_pairs_rejects_single_end_input() {
    let temp = tempdir().expect("tempdir should be created");
    let input = temp.path().join("reads.fastq");
    fs::write(&input, b"@read1\nACGT\n+\nIIII\n").expect("fixture FASTQ should be writable");

    let output = nuclease()
        .args([
            "--in1",
            input.to_str().expect("input path should be UTF-8"),
            "--merge-pairs",
            "-qqq",
        ])
        .output()
        .expect("nuclease should run");

    assert!(
        !output.status.success(),
        "single-end merge-pairs input should fail"
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("--merge-pairs requires paired-end input"),
        "stderr did not explain paired-input requirement: {stderr}"
    );
}

#[test]
fn merge_pairs_rejects_split_paired_output() {
    let temp = tempdir().expect("tempdir should be created");
    let input1 = temp.path().join("reads_1.fastq");
    let input2 = temp.path().join("reads_2.fastq");
    let out1 = temp.path().join("out_1.fastq");
    let out2 = temp.path().join("out_2.fastq");
    fs::write(&input1, b"@read1/1\nAAAA\n+\nIIII\n").expect("read 1 fixture should be writable");
    fs::write(&input2, b"@read1/2\nTTTT\n+\nIIII\n").expect("read 2 fixture should be writable");

    let output = nuclease()
        .args([
            "--in1",
            input1.to_str().expect("read 1 path should be UTF-8"),
            "--in2",
            input2.to_str().expect("read 2 path should be UTF-8"),
            "--merge-pairs",
            "--out1",
            out1.to_str().expect("out1 path should be UTF-8"),
            "--out2",
            out2.to_str().expect("out2 path should be UTF-8"),
            "-qqq",
        ])
        .output()
        .expect("nuclease should run");

    assert!(
        !output.status.success(),
        "split output should fail when merge-pairs is enabled"
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("--merge-pairs cannot be used with split paired output"),
        "stderr did not explain merge-pairs output constraint: {stderr}"
    );
}

#[test]
fn paired_fastq_count_mismatch_reports_source_and_progress() {
    let temp = tempdir().expect("tempdir should be created");
    let input1 = temp.path().join("reads_1.fastq");
    let input2 = temp.path().join("reads_2.fastq");
    fs::write(
        &input1,
        b"@read1/1\nAAAA\n+\nIIII\n@read2/1\nCCCC\n+\nJJJJ\n",
    )
    .expect("read 1 fixture should be writable");
    fs::write(&input2, b"@read1/2\nTTTT\n+\nKKKK\n").expect("read 2 fixture should be writable");

    let output = nuclease()
        .args([
            "--in1",
            input1.to_str().expect("read 1 path should be UTF-8"),
            "--in2",
            input2.to_str().expect("read 2 path should be UTF-8"),
            "--interleaved",
            "--min-length",
            "4",
            "--trim-min-q",
            "0",
            "--min-mean-q",
            "0",
            "-qqq",
        ])
        .output()
        .expect("nuclease should run");

    assert!(!output.status.success(), "mismatched inputs should fail");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("paired FASTQ inputs have different record counts"),
        "stderr did not explain count mismatch: {stderr}"
    );
    assert!(
        stderr.contains("complete_pairs_seen: 1"),
        "stderr did not include completed pair count: {stderr}"
    );
    assert!(
        stderr.contains("local-paired:"),
        "stderr did not include input source label: {stderr}"
    );
}

#[test]
fn paired_fastq_mate_id_mismatch_errors_by_default() {
    let temp = tempdir().expect("tempdir should be created");
    let input1 = temp.path().join("reads_1.fastq");
    let input2 = temp.path().join("reads_2.fastq");
    fs::write(&input1, b"@read1/1\nAAAA\n+\nIIII\n").expect("read 1 fixture should be writable");
    fs::write(&input2, b"@other/2\nTTTT\n+\nKKKK\n").expect("read 2 fixture should be writable");

    let output = nuclease()
        .args([
            "--in1",
            input1.to_str().expect("read 1 path should be UTF-8"),
            "--in2",
            input2.to_str().expect("read 2 path should be UTF-8"),
            "--interleaved",
            "--min-length",
            "4",
            "--trim-min-q",
            "0",
            "--min-mean-q",
            "0",
            "-qqq",
        ])
        .output()
        .expect("nuclease should run");

    assert!(!output.status.success(), "mismatched mate IDs should fail");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("paired FASTQ headers do not agree"),
        "stderr did not explain mate mismatch: {stderr}"
    );
    assert!(
        stderr.contains("read1/1"),
        "stderr missed left header: {stderr}"
    );
    assert!(
        stderr.contains("other/2"),
        "stderr missed right header: {stderr}"
    );
}

#[test]
fn paired_fastq_mate_id_mismatch_warn_drop_continues_with_later_pairs() {
    let temp = tempdir().expect("tempdir should be created");
    let input1 = temp.path().join("reads_1.fastq");
    let input2 = temp.path().join("reads_2.fastq");
    let summary = temp.path().join("summary.json");
    fs::write(&input1, b"@bad1/1\nAAAA\n+\nIIII\n@good/1\nCCCC\n+\nJJJJ\n")
        .expect("read 1 fixture should be writable");
    fs::write(&input2, b"@bad2/2\nTTTT\n+\nKKKK\n@good/2\nGGGG\n+\nLLLL\n")
        .expect("read 2 fixture should be writable");

    let output = nuclease()
        .args([
            "--in1",
            input1.to_str().expect("read 1 path should be UTF-8"),
            "--in2",
            input2.to_str().expect("read 2 path should be UTF-8"),
            "--interleaved",
            "--min-length",
            "4",
            "--trim-min-q",
            "0",
            "--min-mean-q",
            "0",
            "--invalid-fastq-policy",
            "warn-drop",
            "--summary",
            summary.to_str().expect("summary path should be UTF-8"),
        ])
        .output()
        .expect("nuclease should run");

    assert!(
        output.status.success(),
        "warn-drop mate mismatch should continue: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        output.stdout,
        b"@good/1\nCCCC\n+\nJJJJ\n@good/2\nGGGG\n+\nLLLL\n"
    );

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("dropping invalid FASTQ pair"),
        "stderr did not warn about invalid pair: {stderr}"
    );

    let summary_json = fs::read_to_string(summary).expect("summary should be readable");
    let summary: Value = serde_json::from_str(&summary_json).expect("summary should be JSON");
    assert_eq!(summary["pairs_seen"], 2);
    assert_eq!(summary["pairs_emitted"], 1);
    assert_eq!(summary["invalid_pairs"], 1);
    assert_eq!(
        summary["invalid_fastq_samples"][0]["kind"],
        "paired_header_mismatch"
    );
    assert_eq!(summary["invalid_fastq_samples"][0]["left_header"], "bad1/1");
    assert_eq!(
        summary["invalid_fastq_samples"][0]["right_header"],
        "bad2/2"
    );
    assert_eq!(summary["invalid_fastq_samples"][0]["pairs_seen"], 1);
}
