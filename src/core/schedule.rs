// SPDX-License-Identifier: Apache-2.0
// Modified for this library; copyright and attribution notices are in NOTICE.
//! Phase scheduling and cleanup. Cheap cleanup surrounds
//! fast and medium exploration. A cycle ends after medium exploration.

use crate::{
    core::{execution::Executor, model::Model},
    postsolve::tape::Certificate,
    settings::{Progress, PropagationSettings, SparsificationSettings},
};
use std::time::{Duration, Instant};

#[derive(Clone, Copy, Debug)]
pub(crate) struct Limits {
    pub time: Duration,
    /// Maximum newly introduced coefficients in A and P per substitution.
    pub fill: usize,
    pub sparsify: bool,
    pub progress: Progress,
    pub propagation: PropagationSettings,
    pub sparsification: SparsificationSettings,
}

#[derive(Debug, Default)]
pub(crate) struct Stats {
    pub fast_phases: usize,
    pub medium_phases: usize,
    pub parallel_comparisons: usize,
    pub time_limit: bool,
}

impl Model {
    /// Repeat cheap cleanup until no model edits remain. Stable queue
    /// IDs allow later rules to create new singletons during this loop.
    fn cleanup(&mut self) -> Result<(), Certificate> {
        loop {
            let before = self.revision;
            if self.rules.fixed_variables {
                self.fixed_variables();
            }
            if self.rules.cones && !self.cones.is_empty() {
                self.simplify_cones()?;
            }
            if self.rules.empty_columns {
                self.empty_columns()?;
            }
            if self.rules.dual_fixing {
                self.simple_dual_fix()?;
            }
            if self.rules.singleton_rows {
                self.singleton_rows()?;
            }
            if self.rules.empty_rows {
                self.empty_rows()?;
            }
            if self.rules.cones && !self.cones.is_empty() {
                self.simplify_cones()?;
            }
            if self.revision == before {
                break;
            }
        }
        Ok(())
    }

    pub fn run(&mut self, limits: Limits, executor: &Executor) -> Result<Stats, Certificate> {
        let result = self.run_phases(limits, executor);
        if result
            .as_ref()
            .is_ok_and(|stats| !stats.time_limit && self.rules.redundant_bounds)
        {
            self.remove_redundant_bounds();
        }
        result.map_err(|mut certificate| {
            // The tape uses stable model indices, including when a
            // certificate ends exploration before final matrix packing.
            self.postsolve
                .recover(&mut certificate.point, certificate.mode);
            certificate
        })
    }

    fn run_phases(&mut self, limits: Limits, executor: &Executor) -> Result<Stats, Certificate> {
        let start = Instant::now();
        let mut stats = Stats::default();
        let mut fast = true;
        let mut cycle_size = self.work_size();
        let mut cycle_revision = self.revision;
        loop {
            if start.elapsed() >= limits.time {
                stats.time_limit = true;
                break;
            }
            self.cleanup()?;
            let before = self.work_size();
            let before_revision = self.revision;
            if fast {
                stats.fast_phases += 1;
                if self.rules.singleton_columns {
                    self.singleton_columns(limits.fill);
                }
                self.cleanup()?;
                if self.rules.doubleton_equalities {
                    self.doubleton_equalities(limits.fill);
                }
                if self.rules.short_equalities {
                    self.short_equalities(limits.fill, start + limits.time);
                }
                self.cleanup()?;
                // Repeat fast phases under the selected progress policy,
                // then try the more expensive medium rules.
                fast = significant_progress(
                    limits.progress,
                    before,
                    self.work_size(),
                    before_revision != self.revision,
                );
            } else {
                stats.medium_phases += 1;
                if self.rules.bound_propagation {
                    self.propagate_rounds(limits, start + limits.time)?;
                }
                if self.rules.dual_fixing {
                    self.coupled_dual_fix();
                }
                self.cleanup()?;
                if self.rules.short_equalities {
                    self.short_equalities(limits.fill, start + limits.time);
                }
                self.cleanup()?;
                if self.rules.parallel_rows {
                    stats.parallel_comparisons += self.parallel_rows(executor)?;
                }
                if self.rules.parallel_columns {
                    stats.parallel_comparisons += self.parallel_columns(executor)?;
                }
                self.cleanup()?;
                let after = self.work_size();
                if !significant_progress(
                    limits.progress,
                    cycle_size,
                    after,
                    cycle_revision != self.revision,
                ) {
                    break;
                }
                cycle_size = after;
                cycle_revision = self.revision;
                if matches!(limits.progress, Progress::AnyChange) {
                    // A pivot's degree or curvature can change without editing
                    // this equality itself. Revisit those candidates next cycle.
                    for (i, row) in self.rows.iter().enumerate() {
                        if matches!(row, crate::core::model::RowDomain::Linear(b) if b.equality()) {
                            match self.a.row(i).len() {
                                2 => self.queues.doubleton_rows.push(i),
                                3.. => self.queues.short_equalities.push(i),
                                _ => {}
                            }
                        }
                    }
                }
                fast = true;
            }
        }
        if !stats.time_limit && self.rules.equality_dependencies {
            let deadline = start + limits.time;
            if self.equality_dependencies(deadline)? > 0 {
                self.sparsify_cleanup(limits, deadline)?;
            }
            stats.time_limit = Instant::now() >= deadline;
        }
        if !stats.time_limit && limits.sparsify {
            let deadline = start + limits.time;
            if self.sparsify_rows(deadline, limits.sparsification) > 0 {
                self.sparsify_cleanup(limits, deadline)?;
            }
            stats.time_limit = Instant::now() >= deadline;
        }
        Ok(stats)
    }

    /// Bounds can expose rules without immediately shrinking the matrix.
    /// Extra rounds share a configurable work allowance and do not repeat
    /// global parallel scans.
    fn propagate_rounds(&mut self, limits: Limits, deadline: Instant) -> Result<(), Certificate> {
        let mut tightened = self.propagate_bounds()?;
        let mut work = limits
            .propagation
            .work_limit
            .resolve((self.a.nnz() / 4).max(256));
        for _ in 0..limits.propagation.additional_rounds {
            if tightened == 0 || Instant::now() >= deadline {
                break;
            }
            let cost = self
                .queues
                .changed_activities
                .iter()
                .try_fold(0usize, |cost, i| {
                    let cost = cost.saturating_add(self.a.row(i).len());
                    (cost <= work).then_some(cost)
                });
            let Some(cost) = cost else { break };
            if cost == 0 {
                break;
            }
            work -= cost;
            self.cleanup()?;
            if self.rules.singleton_columns {
                self.singleton_columns(limits.fill);
            }
            if self.rules.doubleton_equalities {
                self.doubleton_equalities(limits.fill);
            }
            if self.rules.short_equalities {
                self.short_equalities(limits.fill, deadline);
            }
            self.cleanup()?;
            tightened = self.propagate_bounds()?;
        }
        Ok(())
    }

    // Drain consequences even for size changes below the ordinary 5% cycle
    // threshold. Existing substitution limits still account for A and P fill.
    fn sparsify_cleanup(&mut self, limits: Limits, deadline: Instant) -> Result<(), Certificate> {
        while Instant::now() < deadline {
            let before = self.revision;
            self.cleanup()?;
            if self.rules.bound_propagation {
                self.propagate_bounds()?;
            }
            if self.rules.singleton_columns {
                self.singleton_columns(limits.fill);
            }
            if self.rules.doubleton_equalities {
                self.doubleton_equalities(limits.fill);
            }
            if self.rules.short_equalities {
                self.short_equalities(limits.fill, deadline);
            }
            self.cleanup()?;
            if self.revision == before {
                break;
            }
        }
        Ok(())
    }

    fn work_size(&self) -> usize {
        // In a QP, trading fewer A entries for a much denser Hessian need not
        // be progress. Include P when measuring progress between phases and cycles.
        self.a.nnz() + self.objective.p.nnz()
    }
}

fn significant_progress(policy: Progress, before: usize, after: usize, changed: bool) -> bool {
    match policy {
        Progress::AnyChange => changed,
        Progress::Nonzeros { minimum_reduction } => {
            let threshold = if minimum_reduction.is_nan() {
                0.05
            } else {
                minimum_reduction.clamp(0.0, 1.0)
            };
            before > 0 && (after as f64) < (1.0 - threshold) * (before as f64)
        }
    }
}
