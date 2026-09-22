// SPDX-License-Identifier: Apache-2.0
// Modified for this library; copyright and attribution notices are in NOTICE.
//! Rule families, organized by the model property they act on, and their
//! phase scheduling. Each family checks applicability and uses `Model`
//! mutations to update sparse storage, work queues, and recovery records
//! together. Cheap cleanup surrounds fast and medium exploration; a cycle
//! ends after medium exploration.

mod bound_shift;
mod bounds;
mod cones;
mod dependencies;
pub(crate) mod dominated_columns;
mod dual_fixing;
pub(crate) mod dual_propagation;
mod folding;
mod implied_free;
mod parallel;
mod rows;
mod sparsification;
mod substitution;
mod variables;

use crate::{
    executor::Executor,
    model::{Model, tape::Certificate},
    settings::Progress,
};
use std::time::{Duration, Instant};

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
            if self.settings.rules.fixed_variables {
                self.fixed_variables();
            }
            if self.settings.rules.cones && !self.cones.is_empty() {
                self.simplify_cones()?;
            }
            if self.settings.rules.empty_columns {
                self.empty_columns()?;
            }
            if self.settings.rules.dual_fixing {
                self.simple_dual_fix()?;
            }
            if self.settings.rules.singleton_rows {
                self.singleton_rows()?;
            }
            if self.settings.rules.empty_rows {
                self.empty_rows()?;
            }
            if self.settings.rules.cones && !self.cones.is_empty() {
                self.simplify_cones()?;
            }
            if self.revision == before {
                break;
            }
        }
        Ok(())
    }

    /// Run every phase within `time`, the remaining part of the budget.
    pub fn run(&mut self, time: Duration, executor: &Executor) -> Result<Stats, Certificate> {
        let result = self.run_phases(time, executor);
        result.map_err(|mut certificate| {
            // The tape uses stable model indices, including when a
            // certificate ends exploration before final matrix packing.
            self.postsolve
                .recover(&mut certificate.point, certificate.mode);
            certificate
        })
    }

    fn run_phases(&mut self, time: Duration, executor: &Executor) -> Result<Stats, Certificate> {
        let start = Instant::now();
        let mut stats = Stats::default();
        // Folding needs the original symmetry before asymmetric pivot choices.
        if self.settings.rules.lp_folding && start.elapsed() < time {
            // Conic callers encode variable bounds as singleton rows. Recover
            // those bounds before looking for symmetry, without choosing pivots.
            if self.settings.rules.singleton_rows
                && self.cones.is_empty()
                && self.objective.p.nnz() == 0
            {
                self.singleton_rows()?;
            }
            self.fold_lp(start + time);
        }
        let mut fast = true;
        let mut cycle_size = self.work_size();
        let mut cycle_revision = self.revision;
        loop {
            if start.elapsed() >= time {
                stats.time_limit = true;
                break;
            }
            self.cleanup()?;
            let before = self.work_size();
            let before_revision = self.revision;
            if fast {
                stats.fast_phases += 1;
                if self.settings.rules.singleton_columns {
                    self.singleton_columns();
                }
                self.cleanup()?;
                if self.settings.rules.implied_free_equalities {
                    self.implied_free_equalities(start + time);
                }
                if self.settings.rules.doubleton_equalities {
                    self.doubleton_equalities();
                }
                if self.settings.rules.short_equalities {
                    self.short_equalities(start + time);
                }
                self.cleanup()?;
                // Repeat fast phases under the selected progress policy,
                // then try the more expensive medium rules.
                fast = significant_progress(
                    self.settings.progress,
                    before,
                    self.work_size(),
                    before_revision != self.revision,
                );
            } else {
                stats.medium_phases += 1;
                if self.settings.rules.bound_propagation {
                    self.propagate_rounds(start + time)?;
                }
                if self.settings.rules.dual_fixing {
                    self.coupled_dual_fix();
                }
                self.cleanup()?;
                if self.settings.rules.implied_free_equalities {
                    self.implied_free_equalities(start + time);
                }
                if self.settings.rules.short_equalities {
                    self.short_equalities(start + time);
                }
                self.cleanup()?;
                if self.settings.rules.parallel_rows {
                    stats.parallel_comparisons += self.parallel_rows(executor)?;
                }
                if self.settings.rules.parallel_columns {
                    stats.parallel_comparisons += self.parallel_columns(executor)?;
                }
                if self.settings.rules.dominated_columns
                    && self.settings.dominated_columns.general_search
                {
                    self.dominated_columns()?;
                }
                self.cleanup()?;
                let after = self.work_size();
                if !significant_progress(
                    self.settings.progress,
                    cycle_size,
                    after,
                    cycle_revision != self.revision,
                ) {
                    break;
                }
                cycle_size = after;
                cycle_revision = self.revision;
                if matches!(self.settings.progress, Progress::AnyChange) {
                    // A pivot's degree or curvature can change without editing
                    // this equality itself. Revisit those candidates next cycle.
                    for (i, row) in self.rows.iter().enumerate() {
                        if matches!(row, crate::model::RowDomain::Linear(b) if b.equality()) {
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
        if !stats.time_limit && self.settings.rules.equality_dependencies {
            let deadline = start + time;
            if self.equality_dependencies(deadline)? > 0 {
                self.sparsify_cleanup(deadline)?;
            }
            stats.time_limit = Instant::now() >= deadline;
        }
        if !stats.time_limit && self.settings.rules.sparsification {
            let deadline = start + time;
            if self.sparsify_rows(deadline) > 0 {
                self.sparsify_cleanup(deadline)?;
            }
            stats.time_limit = Instant::now() >= deadline;
        }
        if !stats.time_limit && self.settings.rules.redundant_bounds {
            self.remove_redundant_bounds();
        }
        // Removing implied bounds exposes free pivots and one-sided columns.
        // Visit them before dual propagation, while the independent switches
        // continue to control every consequence rule.
        if !stats.time_limit && self.settings.rules.implied_free_equalities {
            let deadline = start + time;
            if self.implied_free_equalities(deadline) > 0 {
                self.substitution_cleanup(deadline)?;
            }
            stats.time_limit = Instant::now() >= deadline;
        }
        if !stats.time_limit && self.settings.rules.bound_shift {
            let deadline = start + time;
            if self.bound_shift(deadline) > 0 {
                self.substitution_cleanup(deadline)?;
            }
            stats.time_limit = Instant::now() >= deadline;
        }
        // Dual propagation runs once, on the final model: implied-free columns
        // are exposed only now, and one pass here finds what repeated passes
        // in the medium phases would, without their per-phase sweeps.
        if !stats.time_limit && self.settings.rules.dual_propagation {
            let deadline = start + time;
            if self.dual_propagation()? > 0 {
                // Drain the direct consequences only: new equalities feed the
                // substitution rules and fixed columns feed cleanup. A further
                // propagation round or bound sweep would cost more than the
                // few extra reductions it finds here.
                self.substitution_cleanup(deadline)?;
            }
            stats.time_limit = Instant::now() >= deadline;
        }
        Ok(stats)
    }

    /// Bounds can expose rules without immediately shrinking the matrix.
    /// Extra rounds share a configurable work allowance and do not repeat
    /// global parallel scans.
    fn propagate_rounds(&mut self, deadline: Instant) -> Result<(), Certificate> {
        let mut tightened = self.propagate_bounds()?;
        let mut work = self
            .settings
            .propagation
            .work_limit
            .resolve((self.a.nnz() / 4).max(256));
        for _ in 0..self.settings.propagation.additional_rounds {
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
            if self.settings.rules.singleton_columns {
                self.singleton_columns();
            }
            if self.settings.rules.implied_free_equalities {
                self.implied_free_equalities(deadline);
            }
            if self.settings.rules.doubleton_equalities {
                self.doubleton_equalities();
            }
            if self.settings.rules.short_equalities {
                self.short_equalities(deadline);
            }
            self.cleanup()?;
            tightened = self.propagate_bounds()?;
        }
        Ok(())
    }

    // Drain consequences even for size changes below the ordinary 5% cycle
    // threshold. Existing substitution limits still account for A and P fill.
    fn sparsify_cleanup(&mut self, deadline: Instant) -> Result<(), Certificate> {
        while Instant::now() < deadline {
            let before = self.revision;
            self.cleanup()?;
            if self.settings.rules.bound_propagation {
                self.propagate_bounds()?;
            }
            if self.settings.rules.singleton_columns {
                self.singleton_columns();
            }
            if self.settings.rules.implied_free_equalities {
                self.implied_free_equalities(deadline);
            }
            if self.settings.rules.doubleton_equalities {
                self.doubleton_equalities();
            }
            if self.settings.rules.short_equalities {
                self.short_equalities(deadline);
            }
            self.cleanup()?;
            if self.revision == before {
                break;
            }
        }
        Ok(())
    }

    /// Cleanup and substitution only, without propagation rounds.
    fn substitution_cleanup(&mut self, deadline: Instant) -> Result<(), Certificate> {
        while Instant::now() < deadline {
            let before = self.revision;
            self.cleanup()?;
            if self.settings.rules.singleton_columns {
                self.singleton_columns();
            }
            if self.settings.rules.implied_free_equalities {
                self.implied_free_equalities(deadline);
            }
            if self.settings.rules.doubleton_equalities {
                self.doubleton_equalities();
            }
            if self.settings.rules.short_equalities {
                self.short_equalities(deadline);
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
