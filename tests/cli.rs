//! CLI contract tests for `nuclease`.

use std::{
    fs,
    process::{Command, Stdio},
};

use serde_json::Value;
use tempfile::tempdir;

fn nuclease() -> Command {
    Command::new(env!("CARGO_BIN_EXE_nuclease"))
}

#[cfg(unix)]
fn nuclease_with_closed_stderr() -> Command {
    let mut command = Command::new("/bin/sh");
    command
        .args(["-c", "exec 2>&-; exec \"$@\"", "nuclease-shell"])
        .arg(env!("CARGO_BIN_EXE_nuclease"));
    command
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
            "--in",
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
fn missing_local_input_exits_unavailable() {
    let temp = tempdir().expect("tempdir should be created");
    let input = temp.path().join("missing.fastq");

    let output = nuclease()
        .args([
            "--in",
            input.to_str().expect("fixture path should be UTF-8"),
            "--passthrough",
            "-qqq",
        ])
        .output()
        .expect("nuclease should run");

    assert_eq!(output.status.code(), Some(66));
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("failed to open local FASTQ input"));
    assert!(stderr.contains(input.to_str().expect("fixture path should be UTF-8")));
    assert!(stderr.contains("check that the path exists and is readable"));
}

#[test]
fn required_output_create_failure_exits_io() {
    let temp = tempdir().expect("tempdir should be created");
    let input = temp.path().join("reads.fastq");
    let output_path = temp.path().join("missing-parent").join("clean.fastq");
    fs::write(&input, b"@read1\nACGT\n+\nIIII\n").expect("fixture FASTQ should be writable");

    let output = nuclease()
        .args([
            "--in",
            input.to_str().expect("fixture path should be UTF-8"),
            "--passthrough",
            "--out",
            output_path.to_str().expect("output path should be UTF-8"),
            "-qqq",
        ])
        .output()
        .expect("nuclease should run");

    assert_eq!(output.status.code(), Some(74));
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("failed to create required output"));
    assert!(stderr.contains(output_path.to_str().expect("output path should be UTF-8")));
    assert!(stderr.contains("check that the parent directory exists and is writable"));
}

#[test]
fn required_invalid_input_report_create_failure_exits_io() {
    let temp = tempdir().expect("tempdir should be created");
    let input = temp.path().join("reads.fastq");
    let report = temp
        .path()
        .join("missing-parent")
        .join("invalid-input.jsonl");
    fs::write(&input, b"@read1\nACGT\n+\nIIII\n").expect("fixture FASTQ should be writable");

    let output = nuclease()
        .args([
            "--in",
            input.to_str().expect("fixture path should be UTF-8"),
            "--passthrough",
            "--on-invalid-input",
            "skip",
            "--invalid-input-report",
            report.to_str().expect("report path should be UTF-8"),
            "-qqq",
        ])
        .output()
        .expect("nuclease should run");

    assert_eq!(output.status.code(), Some(74));
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("failed to create requested invalid-input JSONL report"));
    assert!(stderr.contains(report.to_str().expect("report path should be UTF-8")));
}

#[cfg(unix)]
#[test]
fn closed_downstream_pipe_exits_io_with_destination_context() {
    let temp = tempdir().expect("tempdir should be created");
    let input = temp.path().join("reads.fastq");
    let reads = b"@read1\nACGT\n+\nIIII\n".repeat(100_000);
    fs::write(&input, reads).expect("fixture FASTQ should be writable");

    let mut child = nuclease()
        .args([
            "--in",
            input.to_str().expect("fixture path should be UTF-8"),
            "--passthrough",
            "-qqq",
        ])
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("nuclease should start");
    drop(
        child
            .stdout
            .take()
            .expect("nuclease stdout should be piped"),
    );

    let output = child
        .wait_with_output()
        .expect("nuclease should report its process status");

    assert_eq!(output.status.code(), Some(74));
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("required output stdout was closed before nuclease finished writing"));
    assert!(stderr.contains("ensure the downstream process consumes the complete output"));
}

#[cfg(unix)]
#[test]
fn closed_stderr_does_not_change_selected_error_status() {
    let temp = tempdir().expect("tempdir should be created");
    let input = temp.path().join("malformed.fastq");
    fs::write(&input, b"@bad\nACGT\n+\nI\n").expect("fixture FASTQ should be writable");

    let status = nuclease_with_closed_stderr()
        .args([
            "--in",
            input.to_str().expect("fixture path should be UTF-8"),
            "--passthrough",
            "-qqq",
        ])
        .stdout(Stdio::null())
        .status()
        .expect("nuclease should run with closed stderr");

    assert_eq!(status.code(), Some(65));
}

#[cfg(unix)]
#[test]
fn plain_progress_is_best_effort_when_stderr_is_closed() {
    let temp = tempdir().expect("tempdir should be created");
    let input = temp.path().join("reads.fastq");
    let output = temp.path().join("clean.fastq");
    fs::write(&input, b"@read1\nACGT\n+\nIIII\n").expect("fixture FASTQ should be writable");

    let status = nuclease_with_closed_stderr()
        .args([
            "--in",
            input.to_str().expect("fixture path should be UTF-8"),
            "--passthrough",
            "--out",
            output.to_str().expect("output path should be UTF-8"),
            "--progress-every",
            "1",
            "-qq",
        ])
        .stdout(Stdio::null())
        .status()
        .expect("nuclease should run with closed stderr");

    assert!(status.success());
    assert_eq!(
        fs::read(output).expect("output should be readable"),
        fs::read(input).expect("input should be readable")
    );
}

#[cfg(unix)]
#[test]
fn human_summary_is_best_effort_when_stderr_is_closed() {
    let temp = tempdir().expect("tempdir should be created");
    let input = temp.path().join("reads.fastq");
    let output = temp.path().join("clean.fastq");
    fs::write(&input, b"@read1\nACGT\n+\nIIII\n").expect("fixture FASTQ should be writable");

    let status = nuclease_with_closed_stderr()
        .args([
            "--in",
            input.to_str().expect("fixture path should be UTF-8"),
            "--passthrough",
            "--out",
            output.to_str().expect("output path should be UTF-8"),
            "--progress-every",
            "0",
        ])
        .stdout(Stdio::null())
        .status()
        .expect("nuclease should run with closed stderr");

    assert!(status.success());
    assert!(output.exists());
}

#[cfg(unix)]
#[test]
fn tracing_is_best_effort_when_stderr_is_closed() {
    let temp = tempdir().expect("tempdir should be created");
    let input = temp.path().join("reads.fastq");
    let output = temp.path().join("clean.fastq");
    fs::write(&input, b"@read1\nACGT\n+\nIIII\n").expect("fixture FASTQ should be writable");

    let status = nuclease_with_closed_stderr()
        .args([
            "--in",
            input.to_str().expect("fixture path should be UTF-8"),
            "--passthrough",
            "--out",
            output.to_str().expect("output path should be UTF-8"),
            "--progress-every",
            "0",
            "-v",
        ])
        .stdout(Stdio::null())
        .status()
        .expect("nuclease should run with closed stderr");

    assert!(status.success());
    assert!(output.exists());
}

#[test]
fn adapter_preset_defaults_to_no_adapter_trimming() {
    let temp = tempdir().expect("tempdir should be created");
    let input = temp.path().join("reads.fastq");
    let fixture = b"@read1\nACGTAGATCGGAAG\n+\nIIIIIIIIIIIIII\n";
    fs::write(&input, fixture).expect("fixture FASTQ should be writable");

    let output = nuclease()
        .args([
            "--in",
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
    assert_eq!(output.stdout, fixture);
}

#[test]
fn illumina_truseq_adapter_preset_trims_suffix_overlap() {
    let temp = tempdir().expect("tempdir should be created");
    let input = temp.path().join("reads.fastq");
    fs::write(&input, b"@read1\nACGTAGATCGGAAG\n+\nIIIIIIIIIIIIII\n")
        .expect("fixture FASTQ should be writable");

    let output = nuclease()
        .args([
            "--in",
            input.to_str().expect("fixture path should be UTF-8"),
            "--adapter-preset",
            "illumina-truseq",
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
            "--in",
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
fn mgi_dnbseq_adapter_preset_trims_suffix_overlap() {
    let temp = tempdir().expect("tempdir should be created");
    let input = temp.path().join("reads.fastq");
    fs::write(&input, b"@read1\nACGTAAGTCGGAGG\n+\nIIIIIIIIIIIIII\n")
        .expect("fixture FASTQ should be writable");

    let output = nuclease()
        .args([
            "--in",
            input.to_str().expect("fixture path should be UTF-8"),
            "--adapter-preset",
            "mgi-dnbseq",
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
fn bin_qualities_without_value_defaults_to_five_bins() {
    let temp = tempdir().expect("tempdir should be created");
    let input = temp.path().join("reads.fastq");
    fs::write(&input, b"@read1\nACG\n+\n!5I\n").expect("fixture FASTQ should be writable");

    let output = nuclease()
        .args([
            "--in",
            input.to_str().expect("fixture path should be UTF-8"),
            "--bin-qualities",
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
    assert_eq!(output.stdout, b"@read1\nACG\n+\n%5E\n");
}

#[test]
fn bin_qualities_accepts_explicit_count_values() {
    let temp = tempdir().expect("tempdir should be created");
    let input = temp.path().join("reads.fastq");
    fs::write(&input, b"@read1\nACG\n+\n!5I\n").expect("fixture FASTQ should be writable");

    let output = nuclease()
        .args([
            "--in",
            input.to_str().expect("fixture path should be UTF-8"),
            "--bin-qualities",
            "2",
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
    assert_eq!(output.stdout, b"@read1\nACG\n+\n+??\n");
}

#[test]
fn bin_qualities_accepts_equals_form() {
    let temp = tempdir().expect("tempdir should be created");
    let input = temp.path().join("reads.fastq");
    fs::write(&input, b"@read1\nACG\n+\n!5I\n").expect("fixture FASTQ should be writable");

    let output = nuclease()
        .args([
            "--in",
            input.to_str().expect("fixture path should be UTF-8"),
            "--bin-qualities=3",
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
    assert_eq!(output.stdout, b"@read1\nACG\n+\n'5B\n");
}

#[test]
fn bin_qualities_rejects_unsupported_count_at_parse_time() {
    let temp = tempdir().expect("tempdir should be created");
    let input = temp.path().join("reads.fastq");
    fs::write(&input, b"@read1\nACGT\n+\nIIII\n").expect("fixture FASTQ should be writable");

    let output = nuclease()
        .args([
            "--in",
            input.to_str().expect("fixture path should be UTF-8"),
            "--bin-qualities",
            "4",
            "-qqq",
        ])
        .output()
        .expect("nuclease should run");

    assert!(
        !output.status.success(),
        "unsupported bin count should fail"
    );
    assert!(
        String::from_utf8_lossy(&output.stderr)
            .contains("quality bin count must be one of 2, 3, or 5"),
        "stderr should explain supported bin counts: {}",
        String::from_utf8_lossy(&output.stderr)
    );
}

#[test]
fn bin_qualities_runs_after_mean_quality_filtering() {
    let temp = tempdir().expect("tempdir should be created");
    let input = temp.path().join("reads.fastq");
    fs::write(&input, b"@read1\nACGT\n+\n!!!!\n").expect("fixture FASTQ should be writable");

    let output = nuclease()
        .args([
            "--in",
            input.to_str().expect("fixture path should be UTF-8"),
            "--bin-qualities",
            "2",
            "--min-length",
            "1",
            "--trim-min-q",
            "0",
            "--min-mean-q",
            "10",
            "-qqq",
        ])
        .output()
        .expect("nuclease should run");

    assert!(
        output.status.success(),
        "nuclease failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(output.stdout, b"");
}

#[test]
fn inert_long_read_catalogs_are_not_cli_adapter_presets() {
    let temp = tempdir().expect("tempdir should be created");
    let input = temp.path().join("reads.fastq");
    let fixture = b"@read1\nACGTAAGAAAGTTGTC\n+\nIIIIIIIIIIIIIIII\n";
    fs::write(&input, fixture).expect("fixture FASTQ should be writable");

    let output = nuclease()
        .args([
            "--in",
            input.to_str().expect("fixture path should be UTF-8"),
            "--adapter-preset",
            "ont-native-barcoding-v14",
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
        !output.status.success(),
        "inert preset should not be accepted"
    );
}

#[test]
fn passthrough_emits_validated_reads_without_preprocessing() {
    let temp = tempdir().expect("tempdir should be created");
    let input = temp.path().join("reads.fastq");
    let summary = temp.path().join("summary.json");
    let fixture = b"@short-adapter\nACGTAGATCGGAAG\n+\nIIIIIIIIIIIIII\n";

    fs::write(&input, fixture).expect("fixture FASTQ should be writable");

    let output = nuclease()
        .args([
            "--in",
            input.to_str().expect("fixture path should be UTF-8"),
            "--passthrough",
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
    assert_eq!(output.stdout, fixture);

    let summary_json = fs::read_to_string(summary).expect("summary should be readable");
    let summary: Value = serde_json::from_str(&summary_json).expect("summary should be JSON");
    assert_eq!(summary["reads_seen"], 1);
    assert_eq!(summary["reads_emitted"], 1);
    assert_eq!(summary["reads_rejected"], 0);
    assert_eq!(
        summary["transform_breakdown"].as_array().map(Vec::len),
        Some(0)
    );
}

#[test]
fn passthrough_rejects_merge_pairs() {
    let temp = tempdir().expect("tempdir should be created");
    let input = temp.path().join("reads.fastq");
    fs::write(&input, b"@read1\nACGT\n+\nIIII\n").expect("fixture FASTQ should be writable");

    let output = nuclease()
        .args([
            "--in",
            input.to_str().expect("fixture path should be UTF-8"),
            "--passthrough",
            "--merge-pairs",
            "-qqq",
        ])
        .output()
        .expect("nuclease should run");

    assert!(
        !output.status.success(),
        "passthrough and merge-pairs should conflict"
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("cannot be used with"),
        "stderr did not explain argument conflict: {stderr}"
    );
}

#[test]
fn skip_invalid_input_requires_a_durable_report() {
    let temp = tempdir().expect("tempdir should be created");
    let input = temp.path().join("reads.fastq");
    fs::write(&input, b"@read1\nACGT\n+\nIIII\n").expect("fixture FASTQ should be writable");

    let output = nuclease()
        .args([
            "--in",
            input.to_str().expect("fixture path should be UTF-8"),
            "--on-invalid-input",
            "skip",
            "-qqq",
        ])
        .output()
        .expect("nuclease should run");

    assert!(
        !output.status.success(),
        "skip without a report should fail"
    );
    assert_eq!(
        output.status.code(),
        Some(2),
        "Clap should own the missing required report relationship"
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("--invalid-input-report"),
        "stderr did not require a durable invalid-input report: {stderr}"
    );
}

#[test]
fn skip_invalid_input_does_not_recover_terminal_parser_error() {
    let temp = tempdir().expect("tempdir should be created");
    let input = temp.path().join("reads.fastq");
    let report = temp.path().join("invalid-input.jsonl");
    fs::write(
        &input,
        b"@good1\nACGT\n+\nIIII\n@bad\nAAAA\n+\nI\n@good2\nTGCA\n+\nJJJJ\n",
    )
    .expect("fixture FASTQ should be writable");

    let output = nuclease()
        .args([
            "--in",
            input.to_str().expect("fixture path should be UTF-8"),
            "--min-length",
            "4",
            "--trim-min-q",
            "0",
            "--min-mean-q",
            "0",
            "--on-invalid-input",
            "skip",
            "--invalid-input-report",
            report.to_str().expect("report path should be UTF-8"),
            "-qqq",
        ])
        .output()
        .expect("nuclease should run");

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        !output.status.success(),
        "terminal parser errors should remain fatal under skip"
    );
    assert_eq!(
        output.status.code(),
        Some(65),
        "malformed local FASTQ should use the stable data-error status"
    );
    assert!(
        stderr.contains("record parser rejected malformed input"),
        "stderr did not include controlled parser diagnostic: {stderr}"
    );
    assert!(
        stderr.contains("admission_policy=skip"),
        "stderr did not include active admission policy: {stderr}"
    );
    assert!(
        stderr.contains("parser_error_kind=unequal_lengths"),
        "stderr did not include needletail error kind: {stderr}"
    );
    assert!(
        !stderr.contains("The application panicked"),
        "parser-level error should not surface as a panic: {stderr}"
    );
}

#[test]
fn invalid_input_report_writes_terminal_record_parse_failure() {
    let temp = tempdir().expect("tempdir should be created");
    let input = temp.path().join("reads.fastq");
    let report = temp.path().join("invalid-input.jsonl");
    fs::write(&input, b"@bad1\nAAAA\n+\nI\n@good\nTGCA\n+\nJJJJ\n")
        .expect("fixture FASTQ should be writable");

    let output = nuclease()
        .args([
            "--in",
            input.to_str().expect("fixture path should be UTF-8"),
            "--min-length",
            "4",
            "--trim-min-q",
            "0",
            "--min-mean-q",
            "0",
            "--on-invalid-input",
            "skip",
            "--invalid-input-report",
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
    assert_eq!(output.status.code(), Some(65));

    let report = fs::read_to_string(report).expect("invalid-input report should be readable");
    let events = report
        .lines()
        .map(|line| serde_json::from_str::<Value>(line).expect("event should be JSON"))
        .collect::<Vec<_>>();
    assert_eq!(events.len(), 1);
    assert_eq!(events[0]["kind"], "record_parse_failure");
    assert_eq!(events[0]["mate"], "single");
    assert_eq!(events[0]["continued"], false);
    assert_eq!(events[0]["parser_kind"], "unequal_lengths");
    assert_eq!(events[0]["line"], 1);
    assert!(
        events[0]["message"]
            .as_str()
            .is_some_and(|message| message.contains("quality length is 1")),
        "event did not preserve parser error message: {:?}",
        events[0]
    );
}

#[test]
fn malformed_fastq_returns_a_controlled_unequal_lengths_error() {
    let temp = tempdir().expect("tempdir should be created");
    let input = temp.path().join("malformed.fastq");
    let report = temp.path().join("invalid-input.jsonl");
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
            "--in",
            input.to_str().expect("fixture path should be UTF-8"),
            "--min-length",
            "1",
            "--trim-min-q",
            "0",
            "--min-mean-q",
            "0",
            "--on-invalid-input",
            "skip",
            "--invalid-input-report",
            report.to_str().expect("report path should be UTF-8"),
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
        !output.status.success(),
        "terminal parser error should fail"
    );
    assert_eq!(output.status.code(), Some(65));
    assert!(
        stderr.contains("record parser rejected malformed input"),
        "malformed input should fail with controlled parser diagnostics: {stderr}"
    );
    assert!(
        stderr.contains("parser_error_kind=unequal_lengths"),
        "malformed input should preserve the parser error kind: {stderr}"
    );
}

#[test]
fn unsupported_fasta_is_not_an_admission_event() {
    let temp = tempdir().expect("tempdir should be created");
    let input = temp.path().join("sequences.fasta");
    let report = temp.path().join("invalid-input.jsonl");
    fs::write(&input, b">sequence-1\nACGT\n").expect("FASTA fixture should be writable");

    let output = nuclease()
        .args([
            "--in",
            input.to_str().expect("fixture path should be UTF-8"),
            "--on-invalid-input",
            "skip",
            "--invalid-input-report",
            report.to_str().expect("report path should be UTF-8"),
            "-qqq",
        ])
        .output()
        .expect("nuclease should run");

    assert!(!output.status.success(), "unsupported FASTA should fail");
    assert_eq!(output.status.code(), Some(65));
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("unsupported input format"),
        "stderr did not classify the source format: {stderr}"
    );
    assert!(
        !stderr.contains("did not provide quality scores"),
        "FASTA should fail before quality extraction: {stderr}"
    );
    assert_eq!(
        fs::read_to_string(report).expect("invalid-input report should be readable"),
        "",
        "source-format failures must not become admission events"
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
fn interleaved_paired_input_streams_pairs_and_writes_summary() {
    let temp = tempdir().expect("tempdir should be created");
    let input = temp.path().join("reads_interleaved.fastq");
    let summary = temp.path().join("summary.json");
    fs::write(
        &input,
        b"@read1/1\nAAAA\n+\nIIII\n@read1/2\nTTTT\n+\nKKKK\n@read2/1\nCCCC\n+\nJJJJ\n@read2/2\nGGGG\n+\nLLLL\n",
    )
    .expect("interleaved fixture should be writable");

    let output = nuclease()
        .args([
            "--in",
            input.to_str().expect("input path should be UTF-8"),
            "--paired",
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
fn interleaved_paired_input_can_merge_pairs() {
    let temp = tempdir().expect("tempdir should be created");
    let input = temp.path().join("reads_interleaved.fastq");
    fs::write(
        &input,
        concat!(
            "@read-1/1\n",
            "ACGTTGCAGTACGATCGTACGGAATTCGCCGATGACTGACCTAGGTCAGTACGATC\n",
            "+\n",
            "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n",
            "@read-1/2\n",
            "GATCGTACTGACCTAGGTCAGTCATCGGCGAATTCCGTACGATCGTACTGCAACGT\n",
            "+\n",
            "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n",
        ),
    )
    .expect("interleaved fixture should be writable");

    let output = nuclease()
        .args([
            "--in",
            input.to_str().expect("input path should be UTF-8"),
            "--paired",
            "--merge-pairs",
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
    let stdout = String::from_utf8(output.stdout).expect("FASTQ output should be UTF-8");
    assert!(
        stdout.starts_with("@read-1\n"),
        "merged read should use normalized pair id as header: {stdout}"
    );
    assert_eq!(
        stdout.matches('\n').count(),
        4,
        "merged interleaved pair should emit one FASTQ record: {stdout}"
    );
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
            "--in",
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
    assert_eq!(
        output.status.code(),
        Some(2),
        "post-resolution semantic usage should match Clap usage status"
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
    assert_eq!(
        output.status.code(),
        Some(2),
        "Clap should own statically expressible output conflicts"
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("--merge-pairs")
            && stderr.contains("cannot be used with")
            && stderr.contains("--out1")
            && stderr.contains("--out2"),
        "stderr did not explain merge-pairs output constraint: {stderr}"
    );
}

#[test]
fn invalid_ena_accession_is_a_clap_value_error() {
    let output = nuclease()
        .args(["--ena", "PRJNA1247874"])
        .output()
        .expect("nuclease should run");

    assert_eq!(output.status.code(), Some(2));
    assert!(
        String::from_utf8_lossy(&output.stderr)
            .contains("use an SRR, ERR, or DRR run accession followed by digits")
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
    assert_eq!(output.status.code(), Some(65));
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
        stderr.contains("present_mate: left") && stderr.contains("header: read2/1"),
        "stderr did not identify the unmatched record: {stderr}"
    );
    assert!(
        stderr.contains("local-paired:"),
        "stderr did not include input source label: {stderr}"
    );
}

#[test]
fn skip_drains_and_reports_every_unmatched_left_record() {
    let temp = tempdir().expect("tempdir should be created");
    let input1 = temp.path().join("reads_1.fastq");
    let input2 = temp.path().join("reads_2.fastq");
    let report = temp.path().join("invalid-input.jsonl");
    let summary = temp.path().join("summary.json");
    fs::write(
        &input1,
        b"@read1/1\nAAAA\n+\nIIII\n@read2/1\nCCCC\n+\nJJJJ\n@read3/1\nGGGG\n+\nLLLL\n",
    )
    .expect("read 1 fixture should be writable");
    fs::write(&input2, b"@read1/2\nTTTT\n+\nKKKK\n").expect("read 2 fixture should be writable");

    let output = nuclease()
        .args([
            "--in1",
            input1.to_str().expect("read 1 path should be UTF-8"),
            "--in2",
            input2.to_str().expect("read 2 path should be UTF-8"),
            "--min-length",
            "4",
            "--trim-min-q",
            "0",
            "--min-mean-q",
            "0",
            "--on-invalid-input",
            "skip",
            "--invalid-input-report",
            report.to_str().expect("report path should be UTF-8"),
            "--summary",
            summary.to_str().expect("summary path should be UTF-8"),
            "-qqq",
        ])
        .output()
        .expect("nuclease should run");

    assert!(
        output.status.success(),
        "skip should drain the longer input: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(output.stderr, b"", "quiet mode should suppress warnings");
    assert_eq!(
        output.stdout,
        b"@read1/1\nAAAA\n+\nIIII\n@read1/2\nTTTT\n+\nKKKK\n"
    );

    let events = fs::read_to_string(report)
        .expect("invalid-input report should be readable")
        .lines()
        .map(|line| serde_json::from_str::<Value>(line).expect("event should be JSON"))
        .collect::<Vec<_>>();
    assert_eq!(events.len(), 2);
    assert_eq!(events[0]["kind"], "missing_mate");
    assert_eq!(events[0]["present_mate"], "left");
    assert_eq!(events[0]["header"], "read2/1");
    assert_eq!(events[1]["header"], "read3/1");
    assert!(events.iter().all(|event| event["continued"] == true));

    let summary: Value =
        serde_json::from_str(&fs::read_to_string(summary).expect("summary should be readable"))
            .expect("summary should be JSON");
    assert_eq!(summary["reads_seen"], 4);
    assert_eq!(summary["bases_seen"], 16);
    assert_eq!(summary["pairs_seen"], 1);
    assert_eq!(summary["invalid_pairs"], 0);
    assert_eq!(
        summary["admission_event_breakdown"][0]["code"],
        "missing_mate"
    );
    assert_eq!(summary["admission_event_breakdown"][0]["count"], 2);
}

#[test]
fn skip_drains_and_reports_every_unmatched_right_record() {
    let temp = tempdir().expect("tempdir should be created");
    let input1 = temp.path().join("reads_1.fastq");
    let input2 = temp.path().join("reads_2.fastq");
    let report = temp.path().join("invalid-input.jsonl");
    fs::write(&input1, b"@read1/1\nAAAA\n+\nIIII\n").expect("read 1 fixture should be writable");
    fs::write(
        &input2,
        b"@read1/2\nTTTT\n+\nKKKK\n@read2/2\nCCCC\n+\nJJJJ\n@read3/2\nGGGG\n+\nLLLL\n",
    )
    .expect("read 2 fixture should be writable");

    let output = nuclease()
        .args([
            "--in1",
            input1.to_str().expect("read 1 path should be UTF-8"),
            "--in2",
            input2.to_str().expect("read 2 path should be UTF-8"),
            "--passthrough",
            "--on-invalid-input",
            "skip",
            "--invalid-input-report",
            report.to_str().expect("report path should be UTF-8"),
            "-qqq",
        ])
        .output()
        .expect("nuclease should run");

    assert!(
        output.status.success(),
        "skip should drain the longer right input: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let events = fs::read_to_string(report)
        .expect("invalid-input report should be readable")
        .lines()
        .map(|line| serde_json::from_str::<Value>(line).expect("event should be JSON"))
        .collect::<Vec<_>>();
    assert_eq!(events.len(), 2);
    assert!(events.iter().all(|event| event["present_mate"] == "right"));
    assert_eq!(events[0]["header"], "read2/2");
    assert_eq!(events[1]["header"], "read3/2");
}

#[test]
fn skip_reports_a_trailing_interleaved_record() {
    let temp = tempdir().expect("tempdir should be created");
    let input = temp.path().join("reads.fastq");
    let report = temp.path().join("invalid-input.jsonl");
    let summary = temp.path().join("summary.json");
    fs::write(&input, b"@read1/1\nAAAA\n+\nIIII\n")
        .expect("interleaved fixture should be writable");

    let output = nuclease()
        .args([
            "--in",
            input.to_str().expect("input path should be UTF-8"),
            "--paired",
            "--passthrough",
            "--on-invalid-input",
            "skip",
            "--invalid-input-report",
            report.to_str().expect("report path should be UTF-8"),
            "--summary",
            summary.to_str().expect("summary path should be UTF-8"),
            "-qqq",
        ])
        .output()
        .expect("nuclease should run");

    assert!(
        output.status.success(),
        "skip should accept one trailing interleaved record: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(output.stdout, b"");

    let event: Value = serde_json::from_str(
        fs::read_to_string(report)
            .expect("invalid-input report should be readable")
            .trim(),
    )
    .expect("event should be JSON");
    assert_eq!(event["kind"], "missing_mate");
    assert_eq!(event["present_mate"], "left");
    assert_eq!(event["header"], "read1/1");
    assert_eq!(event["pairs_seen"], 0);
    assert_eq!(event["continued"], true);

    let summary: Value =
        serde_json::from_str(&fs::read_to_string(summary).expect("summary should be readable"))
            .expect("summary should be JSON");
    assert_eq!(summary["reads_seen"], 1);
    assert_eq!(summary["pairs_seen"], 0);
    assert_eq!(summary["pairs_emitted"], 0);
    assert_eq!(summary["pairs_rejected"], 0);
    assert_eq!(summary["invalid_pairs"], 0);
}

#[test]
fn paired_fastq_mate_id_mismatch_errors_by_default() {
    let temp = tempdir().expect("tempdir should be created");
    let input1 = temp.path().join("reads_1.fastq");
    let input2 = temp.path().join("reads_2.fastq");
    let report = temp.path().join("invalid-input.jsonl");
    fs::write(&input1, b"@read1/1\nAAAA\n+\nIIII\n").expect("read 1 fixture should be writable");
    fs::write(&input2, b"@other/2\nTTTT\n+\nKKKK\n").expect("read 2 fixture should be writable");

    let output = nuclease()
        .args([
            "--in1",
            input1.to_str().expect("read 1 path should be UTF-8"),
            "--in2",
            input2.to_str().expect("read 2 path should be UTF-8"),
            "--min-length",
            "4",
            "--trim-min-q",
            "0",
            "--min-mean-q",
            "0",
            "--invalid-input-report",
            report.to_str().expect("report path should be UTF-8"),
            "-qqq",
        ])
        .output()
        .expect("nuclease should run");

    assert!(!output.status.success(), "mismatched mate IDs should fail");
    assert_eq!(output.status.code(), Some(65));
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("paired input identifiers do not agree"),
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

    let event: Value = serde_json::from_str(
        fs::read_to_string(report)
            .expect("invalid-input report should be readable")
            .trim(),
    )
    .expect("event should be JSON");
    assert_eq!(event["kind"], "pair_identifier_mismatch");
    assert_eq!(event["continued"], false);
    assert_eq!(event["left_header"], "read1/1");
    assert_eq!(event["right_header"], "other/2");
}

#[test]
fn paired_identifier_mismatch_skip_continues_with_later_pairs() {
    let temp = tempdir().expect("tempdir should be created");
    let input1 = temp.path().join("reads_1.fastq");
    let input2 = temp.path().join("reads_2.fastq");
    let summary = temp.path().join("summary.json");
    let report = temp.path().join("invalid-input.jsonl");
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
            "--min-length",
            "4",
            "--trim-min-q",
            "0",
            "--min-mean-q",
            "0",
            "--on-invalid-input",
            "skip",
            "--invalid-input-report",
            report.to_str().expect("report path should be UTF-8"),
            "--summary",
            summary.to_str().expect("summary path should be UTF-8"),
        ])
        .output()
        .expect("nuclease should run");

    assert!(
        output.status.success(),
        "skip should continue after a pair identifier mismatch: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        output.stdout,
        b"@good/1\nCCCC\n+\nJJJJ\n@good/2\nGGGG\n+\nLLLL\n"
    );

    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("skipping paired input with mismatched identifiers"),
        "stderr did not warn about invalid pair: {stderr}"
    );

    let summary_json = fs::read_to_string(summary).expect("summary should be readable");
    let summary: Value = serde_json::from_str(&summary_json).expect("summary should be JSON");
    assert_eq!(summary["pairs_seen"], 2);
    assert_eq!(summary["pairs_emitted"], 1);
    assert_eq!(summary["invalid_pairs"], 1);
    assert_eq!(
        summary["admission_samples"][0]["kind"],
        "pair_identifier_mismatch"
    );
    assert_eq!(summary["admission_samples"][0]["left_header"], "bad1/1");
    assert_eq!(summary["admission_samples"][0]["right_header"], "bad2/2");
    assert_eq!(summary["admission_samples"][0]["pairs_seen"], 1);
    assert_eq!(
        summary["admission_event_breakdown"][0]["code"],
        "pair_identifier_mismatch"
    );

    let report = fs::read_to_string(report).expect("invalid-input report should be readable");
    let event: Value = serde_json::from_str(report.trim()).expect("event should be JSON");
    assert_eq!(event["kind"], "pair_identifier_mismatch");
    assert_eq!(event["continued"], true);
}
