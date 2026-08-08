//! Progress rendering for long-running streaming jobs.

use std::{
    io::{self, Write as _},
    time::Instant,
};

use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};

use crate::observer::RunObserver;

/// Progress rendering mode selected from CLI quietness and terminal capabilities.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ProgressMode {
    /// Continuously refresh a single stderr status line.
    Live,
    /// Emit plain periodic stderr lines without terminal cursor control.
    Plain,
    /// Disable progress reporting entirely.
    Off,
}

/// Emits progress updates after a configurable number of records have been seen.
pub struct ProgressReporter {
    mode: ProgressMode,
    every: u64,
    next_threshold: u64,
    started_at: Instant,
    bar: Option<ProgressBar>,
}

impl ProgressReporter {
    /// Construct a reporter that updates every `every` seen records in the selected mode.
    pub fn new(mode: ProgressMode, every: u64) -> Self {
        let bar = (mode == ProgressMode::Live).then(|| {
            let bar = ProgressBar::with_draw_target(None, ProgressDrawTarget::stderr());
            bar.set_style(
                ProgressStyle::default_spinner()
                    .tick_strings(&["⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏"])
                    .template("{spinner} {msg}")
                    .expect("progress template should be valid"),
            );
            bar
        });

        Self {
            mode,
            every,
            next_threshold: every,
            started_at: Instant::now(),
            bar,
        }
    }

    /// Emit a progress update if the next threshold has been crossed.
    pub fn maybe_report(&mut self, observer: &RunObserver) {
        if self.mode == ProgressMode::Off
            || self.every == 0
            || observer.reads_seen < self.next_threshold
        {
            return;
        }

        let message = render_message(observer, self.started_at.elapsed().as_secs_f64());

        match self.mode {
            ProgressMode::Live => {
                if let Some(bar) = &self.bar {
                    bar.set_message(message);
                    bar.tick();
                }
            }
            ProgressMode::Plain => {
                let _ = writeln!(io::stderr(), "{message}");
            }
            ProgressMode::Off => {}
        }

        self.next_threshold = self.next_threshold.saturating_add(self.every);
    }

    /// Clear any live terminal progress rendering before final summary output.
    pub fn finish(&self) {
        if let Some(bar) = &self.bar {
            bar.finish_and_clear();
        }
    }
}

fn render_message(observer: &RunObserver, elapsed_seconds: f64) -> String {
    let reads_per_second = rate(observer.reads_seen, elapsed_seconds);
    let bases_per_second = rate(observer.bases_seen, elapsed_seconds);
    let (retained_fraction, retained_unit) = if observer.pairs_seen > 0 {
        (
            fraction(observer.pairs_emitted, observer.pairs_seen),
            "pairs",
        )
    } else {
        (
            fraction(observer.reads_emitted, observer.reads_seen),
            "reads",
        )
    };
    let top_rejection = top_rejection_reason(observer).map_or_else(
        || "none".to_owned(),
        |(code, count)| format!("{code}:{count}"),
    );

    format!(
        "reads={} emitted={} rejected={} ({:.1}% {} retained) bases={} emitted_bases={} {:.0} reads/s {:.1} bases/s top_reject={}",
        observer.reads_seen,
        observer.reads_emitted,
        observer.reads_rejected,
        retained_fraction * 100.0,
        retained_unit,
        observer.bases_seen,
        observer.bases_emitted,
        reads_per_second,
        bases_per_second,
        top_rejection,
    )
}

fn top_rejection_reason(observer: &RunObserver) -> Option<(&'static str, u64)> {
    observer
        .rejection_counts
        .iter()
        .max_by(|(left_code, left_count), (right_code, right_count)| {
            left_count
                .cmp(right_count)
                .then_with(|| right_code.cmp(left_code))
        })
        .map(|(code, count)| (*code, *count))
}

fn fraction(numerator: u64, denominator: u64) -> f64 {
    if denominator == 0 {
        0.0
    } else {
        u64_to_f64(numerator) / u64_to_f64(denominator)
    }
}

fn rate(total: u64, elapsed_seconds: f64) -> f64 {
    if elapsed_seconds <= f64::EPSILON {
        0.0
    } else {
        u64_to_f64(total) / elapsed_seconds
    }
}

#[expect(
    clippy::cast_precision_loss,
    reason = "progress rates and fractions do not require integer-exact floating-point values"
)]
fn u64_to_f64(value: u64) -> f64 {
    value as f64
}

#[cfg(test)]
mod tests {
    use crate::observer::RunObserver;

    use super::{ProgressMode, ProgressReporter, rate, render_message};

    #[test]
    fn rate_is_zero_when_elapsed_is_zero() {
        assert!((rate(42, 0.0) - 0.0).abs() < f64::EPSILON);
    }

    #[test]
    fn rate_divides_total_by_elapsed() {
        assert!((rate(100, 4.0) - 25.0).abs() < f64::EPSILON);
    }

    #[test]
    fn merged_pair_progress_reports_pair_retention() {
        let mut observer = RunObserver::new("test".to_owned());
        observer.reads_seen = 2;
        observer.reads_emitted = 1;
        observer.pairs_seen = 1;
        observer.pairs_emitted = 1;

        let message = render_message(&observer, 1.0);

        assert!(
            message.contains("(100.0% pairs retained)"),
            "unexpected progress message: {message}"
        );
    }

    #[test]
    fn single_end_progress_reports_read_retention() {
        let mut observer = RunObserver::new("test".to_owned());
        observer.reads_seen = 4;
        observer.reads_emitted = 3;

        let message = render_message(&observer, 1.0);

        assert!(
            message.contains("(75.0% reads retained)"),
            "unexpected progress message: {message}"
        );
    }

    #[test]
    fn live_reporter_constructs_progress_bar() {
        let reporter = ProgressReporter::new(ProgressMode::Live, 100);
        assert_eq!(reporter.mode, ProgressMode::Live);
        assert!(reporter.bar.is_some());
    }
}
