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
fn warn_drop_invalid_fastq_drops_bad_record_and_keeps_streaming() {
    let temp = tempdir().expect("tempdir should be created");
    let input = temp.path().join("reads.fastq");
    let summary = temp.path().join("summary.json");
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
            "--summary",
            summary.to_str().expect("summary path should be UTF-8"),
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
        b"@good1\nACGT\n+\nIIII\n@good2\nTGCA\n+\nJJJJ\n"
    );

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("dropping invalid FASTQ record"),
        "stderr did not include invalid-record warning: {stderr}"
    );
    assert!(
        stderr.contains("bad"),
        "stderr did not include invalid record header: {stderr}"
    );

    let summary_json = fs::read_to_string(summary).expect("summary should be readable");
    let summary: Value = serde_json::from_str(&summary_json).expect("summary should be JSON");
    assert_eq!(summary["reads_seen"], 3);
    assert_eq!(summary["reads_emitted"], 2);
    assert_eq!(summary["invalid_reads"], 1);
    assert_eq!(
        summary["invalid_fastq_samples"][0]["kind"],
        "sequence_quality_length_mismatch"
    );
    assert_eq!(summary["invalid_fastq_samples"][0]["header"], "bad");
    assert_eq!(summary["invalid_fastq_samples"][0]["sequence_len"], 4);
    assert_eq!(summary["invalid_fastq_samples"][0]["quality_len"], 1);
    assert_eq!(summary["invalid_fastq_samples_truncated"], false);
}

#[test]
fn invalid_fastq_report_writes_every_invalid_event_as_jsonl() {
    let temp = tempdir().expect("tempdir should be created");
    let input = temp.path().join("reads.fastq");
    let report = temp.path().join("invalid-fastq.jsonl");
    fs::write(
        &input,
        b"@bad1\nAAAA\n+\nI\n@good\nTGCA\n+\nJJJJ\n@bad2\nCCCC\n+\nJJ\n",
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
            "silent-drop",
            "--invalid-fastq-report",
            report.to_str().expect("report path should be UTF-8"),
            "-qqq",
        ])
        .output()
        .expect("nuclease should run");

    assert!(
        output.status.success(),
        "nuclease failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(output.stdout, b"@good\nTGCA\n+\nJJJJ\n");

    let report = fs::read_to_string(report).expect("invalid FASTQ report should be readable");
    let events = report
        .lines()
        .map(|line| serde_json::from_str::<Value>(line).expect("event should be JSON"))
        .collect::<Vec<_>>();
    assert_eq!(events.len(), 2);
    assert_eq!(events[0]["header"], "bad1");
    assert_eq!(events[0]["sequence_len"], 4);
    assert_eq!(events[0]["quality_len"], 1);
    assert_eq!(events[1]["header"], "bad2");
    assert_eq!(events[1]["quality_len"], 2);
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
        output.status.success() || stderr.contains("FASTQ parser failed"),
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
