# Nuclease project justfile.
# All repeating commands should be recipes here.

# Show the available recipes by default.
default:
    @just --list

# Interactively choose a recipe to run.
choose:
    @just --choose

# List all available recipes.
list:
    @just --list

alias l := list
alias ls := list
alias ch := choose

# Run all pre-commit checks.
check: fmt-check lint test doc-check
    @echo "✅ All checks passed"

# Run all checks on the full codebase.
check-all: fmt-check lint-all test-all doc-check
    @echo "✅ All checks passed on full codebase"

# Check formatting without modifying files.
fmt-check:
    cargo fmt --all -- --check

# Apply formatting fixes.
fmt:
    cargo fmt --all

alias f := fmt

# Run clippy with warnings treated as errors.
lint:
    cargo clippy --all-targets --all-features -- -D warnings

alias cl := lint

# Run clippy on all files.
lint-all:
    cargo clippy --all-targets --all-features -- -D warnings

# Lint shell scripts with shellcheck.
lint-shell:
    shellcheck INSTALL.sh

# Run tests with nextest.
test:
    cargo nextest run --all-features --no-tests=pass

alias tr := test

# Run all tests including ignored tests.
test-all:
    cargo nextest run --all-features --run-ignored all --no-tests=pass

# Run tests with verbose output.
test-verbose:
    cargo nextest run --all-features --no-capture --no-tests=pass

# Build the Rust crate in debug mode.
build:
    cargo build

alias b := build

# Install the nuclease binary from this checkout.
install:
    cargo install --path=.

alias i := install

# Build the Rust crate in release mode.
build-release:
    cargo build --release

alias r := build-release

# Typecheck the Rust crate.
check-compile:
    cargo check --all-targets --all-features

alias cr := check-compile

# Check that documentation builds without errors.
doc-check:
    cargo doc --no-deps --document-private-items

# Generate and open documentation.
doc:
    cargo doc --no-deps --open

# Compatibility alias for the nt-terroir-era Rust check recipe name.
rust-check-all: check

alias rca := rust-check-all

# Show current jj status.
status:
    jj status

# Show jj log.
log:
    jj log

# Clean build artifacts.
clean:
    cargo clean

# Update dependencies.
update:
    cargo update

# Run the application help as a small CLI smoke test.
smoke:
    cargo run -- --help

# Prepare a commit: run all checks, then show status.
prepare-commit: check
    @echo ""
    @echo "Ready to commit. Run: jj commit -m 'your message'"
    @jj status

# Prepare for push: run full checks.
prepare-push: check-all
    @echo ""
    @echo "Ready to push. Run: jj git push"
