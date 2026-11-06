//! Offline checks of the installer's release selection and installation paths.

#![cfg(unix)]

use std::{
    fs,
    os::unix::fs::{PermissionsExt, symlink},
    path::Path,
    process::{Command, Output},
};

use tempfile::{TempDir, tempdir};

fn executable(path: &Path, contents: &str) {
    fs::write(path, contents).expect("write executable");
    fs::set_permissions(path, fs::Permissions::from_mode(0o755)).expect("executable permissions");
}

fn installer(release_url: &str) -> (TempDir, Command) {
    let dir = tempdir().expect("installer test directory");
    let bin = dir.path().join("bin");
    let scratch = dir.path().join("scratch");
    fs::create_dir(&bin).expect("mock command directory");
    fs::create_dir(&scratch).expect("installer temporary directory");
    // Only these local utilities and the mocks are reachable through PATH.
    // GNU tar invokes gzip to extract compressed archives.
    for name in [
        "mktemp", "mkdir", "chmod", "mv", "rm", "tar", "gzip", "cp", "cat",
    ] {
        let program = ["/usr/bin", "/bin"]
            .into_iter()
            .map(|parent| Path::new(parent).join(name))
            .find(|path| path.is_file())
            .expect("installer utility");
        symlink(program, bin.join(name)).expect("utility link");
    }
    executable(
        &bin.join("uname"),
        r#"#!/bin/bash
set -euo pipefail
case "$1" in
  -s) printf '%s\n' "${TEST_OS:-Linux}" ;;
  -m) printf '%s\n' "${TEST_ARCH:-x86_64}" ;;
  *) exit 99 ;;
esac
"#,
    );
    executable(
        &bin.join("curl"),
        r#"#!/bin/bash
set -euo pipefail
case "${@: -1}" in
  https://github.com/nrminor/nuclease/releases/latest)
    if [[ "${0##*/}" == wget ]]; then
      printf '%s\n' "$RELEASE_URL" >&2
    else
      printf '%s\n' "$RELEASE_URL"
    fi
    exit "${LOOKUP_EXIT:-0}"
    ;;
  https://github.com/nrminor/nuclease/releases/download/*)
    printf '%s\n' "$@" > "$TEST_ROOT/download.args"
    if [[ -n "${BINARY_ARCHIVE:-}" ]]; then
      while [[ $# -gt 0 ]]; do
        if [[ "$1" == -o || "$1" == -O ]]; then cp "$BINARY_ARCHIVE" "$2"; exit 0; fi
        shift
      done
      exit 99
    fi
    exit 22
    ;;
  *) echo 'Unexpected network request' >&2; exit 99 ;;
esac
"#,
    );
    executable(
        &bin.join("cargo"),
        r#"#!/bin/bash
set -euo pipefail
printf '%s\n' "$@" > "$TEST_ROOT/cargo.args"
if [[ "${RUSTFLAGS+x}" == x ]]; then
  printf '%s' "$RUSTFLAGS" > "$TEST_ROOT/cargo.rustflags"
fi
if [[ "${CARGO_ENCODED_RUSTFLAGS+x}" == x ]]; then
  printf '%s' "$CARGO_ENCODED_RUSTFLAGS" > "$TEST_ROOT/cargo.encoded-rustflags"
fi
if [[ "${CARGO_EXIT:-0}" != 0 ]]; then exit "$CARGO_EXIT"; fi
while [[ $# -gt 0 ]]; do
  if [[ "$1" == --root ]]; then
    mkdir -p "$2/bin"
    printf '#!/bin/sh\nexit 0\n' > "$2/bin/nuclease"
    chmod +x "$2/bin/nuclease"
    exit 0
  fi
  shift
done
exit 99
"#,
    );
    let mut command = Command::new("/bin/bash");
    command
        .arg(Path::new(env!("CARGO_MANIFEST_DIR")).join("INSTALL.sh"))
        .current_dir(dir.path())
        .env_clear()
        .env("PATH", &bin)
        .env("HOME", dir.path().join("home"))
        .env("TMPDIR", scratch)
        .env("TEST_ROOT", dir.path())
        .env("RELEASE_URL", release_url);
    (dir, command)
}

fn succeeds(output: &Output) {
    assert!(
        output.status.success(),
        "status={}; stdout={}; stderr={}",
        output.status,
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr),
    );
}

#[test]
fn source_fallback_installs_the_selected_release_with_locked_dependencies() {
    let (dir, mut command) = installer("https://github.com/nrminor/nuclease/releases/tag/v0.8.0");
    succeeds(&command.output().expect("run installer"));
    let arguments = fs::read_to_string(dir.path().join("cargo.args")).expect("Cargo arguments");
    let arguments: Vec<_> = arguments.lines().collect();
    assert!(arguments.windows(2).any(|pair| pair == ["--tag", "v0.8.0"]));
    assert!(arguments.contains(&"--locked"));
    assert!(
        arguments
            .windows(2)
            .any(|pair| pair == ["--git", "https://github.com/nrminor/nuclease.git"])
    );
    assert!(dir.path().join("home/.local/bin/nuclease").is_file());
    let download =
        fs::read_to_string(dir.path().join("download.args")).expect("download arguments");
    assert!(download.contains("/releases/download/v0.8.0/nuclease-"));
}

#[test]
fn source_fallback_enables_avx2_on_x86_and_preserves_compiler_flags() {
    for (os, arch) in [
        ("Linux", "x86_64"),
        ("Darwin", "x86_64"),
        ("Linux", "amd64"),
    ] {
        for (plain, encoded, expected_plain, expected_encoded) in [
            (None, None, Some("-Ctarget-feature=+avx2"), None),
            (Some(""), None, Some("-Ctarget-feature=+avx2"), None),
            (
                Some("-C debuginfo=1"),
                None,
                Some("-C debuginfo=1 -Ctarget-feature=+avx2"),
                None,
            ),
            (
                Some("-C debuginfo=1"),
                Some("--cfg\u{1f}build_label=\"with spaces\""),
                Some("-C debuginfo=1"),
                Some("--cfg\u{1f}build_label=\"with spaces\"\u{1f}-Ctarget-feature=+avx2"),
            ),
            (
                Some("-C debuginfo=1"),
                Some(""),
                Some("-C debuginfo=1"),
                Some("-Ctarget-feature=+avx2"),
            ),
            (None, Some(""), None, Some("-Ctarget-feature=+avx2")),
        ] {
            let (dir, mut command) =
                installer("https://github.com/nrminor/nuclease/releases/tag/v0.8.0");
            command.env("TEST_OS", os).env("TEST_ARCH", arch);
            if let Some(flags) = plain {
                command.env("RUSTFLAGS", flags);
            }
            if let Some(flags) = encoded {
                command.env("CARGO_ENCODED_RUSTFLAGS", flags);
            }
            succeeds(&command.output().expect("installer"));
            for (file, expected) in [
                ("cargo.rustflags", expected_plain),
                ("cargo.encoded-rustflags", expected_encoded),
            ] {
                let path = dir.path().join(file);
                if let Some(expected) = expected {
                    assert_eq!(
                        fs::read_to_string(path).expect("compiler flags"),
                        expected,
                        "{os} {arch}: {file}",
                    );
                } else {
                    assert!(!path.exists(), "{os} {arch}: unexpected {file}");
                }
            }
        }
    }
}

#[test]
fn source_fallback_leaves_arm_compiler_flags_unchanged() {
    for (os, arch) in [("Linux", "aarch64"), ("Darwin", "arm64")] {
        for (plain, encoded) in [
            (None, None),
            (Some(""), Some("")),
            (Some("-C debuginfo=1"), None),
            (
                Some("-C debuginfo=1"),
                Some("--cfg\u{1f}build_label=\"with spaces\""),
            ),
        ] {
            let (dir, mut command) =
                installer("https://github.com/nrminor/nuclease/releases/tag/v0.8.0");
            command.env("TEST_OS", os).env("TEST_ARCH", arch);
            if let Some(flags) = plain {
                command.env("RUSTFLAGS", flags);
            }
            if let Some(flags) = encoded {
                command.env("CARGO_ENCODED_RUSTFLAGS", flags);
            }
            succeeds(&command.output().expect("installer"));
            for (file, expected) in [
                ("cargo.rustflags", plain),
                ("cargo.encoded-rustflags", encoded),
            ] {
                let path = dir.path().join(file);
                if let Some(expected) = expected {
                    assert_eq!(
                        fs::read_to_string(path).expect("compiler flags"),
                        expected,
                        "{os} {arch}: {file}",
                    );
                } else {
                    assert!(!path.exists(), "{os} {arch}: unexpected {file}");
                }
            }
        }
    }
}

#[test]
fn unavailable_or_unexpected_release_destination_stops_before_installation() {
    for (url, exit) in [
        (
            "https://github.com/nrminor/nuclease/releases/tag/v0.8.0",
            "22",
        ),
        ("", "0"),
        ("https://github.com/nrminor/nuclease/releases/latest", "0"),
        ("https://github.com/nrminor/nuclease/releases/tag/main", "0"),
        ("https://example.invalid/releases/tag/v0.8.0", "0"),
        ("https://github.com/other/project/releases/tag/v0.8.0", "0"),
        ("/nrminor/nuclease/releases/tag/v0.8.0", "0"),
        (
            "https://github.com/nrminor/nuclease/releases/tag/v0.8.0%2Bbuild",
            "0",
        ),
    ] {
        let (dir, mut command) = installer(url);
        let output = command
            .env("LOOKUP_EXIT", exit)
            .output()
            .expect("installer");
        assert!(!output.status.success(), "accepted {url}");
        assert!(
            String::from_utf8_lossy(&output.stderr).contains("Could not determine latest release"),
            "missing lookup diagnostic for {url}: {output:?}"
        );
        assert!(!dir.path().join("download.args").exists());
        assert!(!dir.path().join("cargo.args").exists());
        assert!(!dir.path().join("home/.local/bin/nuclease").exists());
    }
}

#[test]
fn wget_release_discovery_requires_a_successful_canonical_redirect() {
    for (headers, exit, accepted) in [
        (
            "  HTTP/1.1 302 Found\r\n  Location: https://github.com/nrminor/nuclease/releases/latest\r\n  HTTP/1.1 302 Found\r\n  lOcAtIoN: https://github.com/nrminor/nuclease/releases/tag/v0.8.0\r\n  HTTP/1.1 200 OK\r\n",
            "0",
            true,
        ),
        ("  HTTP/1.1 200 OK\n", "0", false),
        (
            "  Location: https://github.com/nrminor/nuclease/releases/tag/v0.8.0 unexpected\r\n",
            "0",
            false,
        ),
        (
            "  Location: /nrminor/nuclease/releases/tag/v0.8.0\n",
            "0",
            false,
        ),
        (
            "  Location: https://github.com/other/project/releases/tag/v0.8.0\n",
            "0",
            false,
        ),
        (
            "  Location: https://github.com/nrminor/nuclease/releases/tag/v0.8.0\n  Location: https://github.com/login\n",
            "0",
            false,
        ),
        (
            "  Location: https://github.com/nrminor/nuclease/releases/tag/v0.8.0\n  HTTP/1.1 404 Not Found\n",
            "8",
            false,
        ),
    ] {
        let (dir, mut command) = installer(headers);
        fs::rename(dir.path().join("bin/curl"), dir.path().join("bin/wget"))
            .expect("wget-only environment");
        let output = command
            .env("LOOKUP_EXIT", exit)
            .output()
            .expect("installer");
        if accepted {
            succeeds(&output);
            let arguments =
                fs::read_to_string(dir.path().join("cargo.args")).expect("Cargo arguments");
            let arguments: Vec<_> = arguments.lines().collect();
            assert!(arguments.windows(2).any(|pair| pair == ["--tag", "v0.8.0"]));
            assert!(arguments.contains(&"--locked"));
            assert!(dir.path().join("home/.local/bin/nuclease").is_file());
        } else {
            assert!(!output.status.success(), "accepted {headers}");
            assert!(
                String::from_utf8_lossy(&output.stderr)
                    .contains("Could not determine latest release"),
                "missing lookup diagnostic for {headers}: {output:?}"
            );
            assert!(!dir.path().join("download.args").exists());
            assert!(!dir.path().join("cargo.args").exists());
            assert!(!dir.path().join("home/.local/bin/nuclease").exists());
        }
    }
}

#[test]
fn source_fallback_preserves_conda_install_paths_with_spaces() {
    let (dir, mut command) =
        installer("https://github.com/nrminor/nuclease/releases/tag/v0.8.0-rc.1+build.2");
    let prefix = dir.path().join("conda environment");
    succeeds(
        &command
            .env("CONDA_PREFIX", &prefix)
            .output()
            .expect("installer"),
    );
    assert!(prefix.join("bin/nuclease").is_file());
    assert!(!dir.path().join("home/.local/bin/nuclease").exists());
    let arguments = fs::read_to_string(dir.path().join("cargo.args")).expect("Cargo arguments");
    let arguments: Vec<_> = arguments.lines().collect();
    assert!(
        arguments
            .windows(2)
            .any(|pair| pair == ["--tag", "v0.8.0-rc.1+build.2"])
    );
    assert!(arguments.contains(&"--locked"));
}

#[test]
fn source_build_failure_is_reported_without_installing() {
    let (dir, mut command) = installer("https://github.com/nrminor/nuclease/releases/tag/v0.8.0");
    let output = command.env("CARGO_EXIT", "1").output().expect("installer");
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("Failed to install from source"));
    assert!(!dir.path().join("home/.local/bin/nuclease").exists());
}

#[test]
fn available_release_binary_does_not_invoke_cargo() {
    let (dir, mut command) = installer("https://github.com/nrminor/nuclease/releases/tag/v0.8.0");
    let fixture = dir.path().join("fixture");
    fs::create_dir(&fixture).expect("fixture directory");
    let contents = "#!/bin/sh\nprintf 'nuclease 0.8.0\\n'\n";
    executable(&fixture.join("nuclease"), contents);
    let archive = dir.path().join("release.tar.gz");
    succeeds(
        &Command::new("tar")
            .arg("-czf")
            .arg(&archive)
            .arg("-C")
            .arg(&fixture)
            .arg("nuclease")
            .output()
            .expect("fixture archive"),
    );
    let output = command
        .env("BINARY_ARCHIVE", archive)
        .output()
        .expect("installer");
    succeeds(&output);
    assert!(
        !dir.path().join("cargo.args").exists(),
        "unexpected source fallback: {output:?}",
    );
    assert_eq!(
        fs::read_to_string(dir.path().join("home/.local/bin/nuclease")).expect("installed binary"),
        contents,
    );
    let download =
        fs::read_to_string(dir.path().join("download.args")).expect("download arguments");
    assert!(download.contains("/releases/download/v0.8.0/nuclease-"));
}
