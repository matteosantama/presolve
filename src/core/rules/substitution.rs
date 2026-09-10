// SPDX-License-Identifier: Apache-2.0
// Modified for this library; copyright and attribution notices are in NOTICE.
//! Substitute singleton columns and doubleton or short equalities.
//! Shared model mutations account for Hessian fill and gradient recovery.

use crate::{
    core::model::{Model, RowDomain},
    problem::Bounds,
};
use std::time::Instant;

impl Model {
    /// Equality substitution with bounded stable alternatives and sparse-work scoring.
    pub fn short_equalities(&mut self, max_fill: usize, deadline: Instant) {
        let options = self.equalities;
        let relative = if options.relative_pivot > 0.0 && options.relative_pivot <= 1.0 {
            options.relative_pivot
        } else {
            1.0
        };
        let mut work = options.work_limit.resolve(self.a.nnz().saturating_mul(2));
        let mut candidates = Vec::new();
        for i in self.queues.short_equalities.take_round() {
            if work == 0 || Instant::now() >= deadline {
                if matches!(options.work_limit, crate::settings::WorkLimit::Default) {
                    break;
                }
                self.queues.short_equalities.push(i);
                continue;
            }
            let RowDomain::Linear(bounds) = self.rows[i] else {
                continue;
            };
            if !bounds.equality() {
                continue;
            }
            self.equality_stats.rows_examined += 1;
            let row = self.a.row(i);
            let length = row.len();
            if length < 3 || length > options.max_row_length || options.max_pivot_attempts == 0 {
                self.equality_stats.structural_rejections += 1;
                continue;
            }
            let largest = row.iter().map(|(_, a)| a.abs()).fold(0.0, f64::max);
            candidates.clear();
            for (j, a) in row {
                let degree = self.a.column(j).len();
                if (options.require_free_variable && self.bounds[j] != Bounds::FREE)
                    || (options.require_linear_variable && !self.objective.p.column(j).is_empty())
                    || degree < options.min_column_length
                    || degree > options.max_column_length
                {
                    self.equality_stats.structural_rejections += 1;
                    continue;
                }
                if a.abs() < largest && (relative == 1.0 || a.abs() / largest < relative) {
                    self.equality_stats.pivot_rejections += 1;
                    continue;
                }
                let diagonal = usize::from(self.objective.p.get(j, j) != 0.0);
                let quadratic = (self.objective.p.column(j).len() - diagonal)
                    .saturating_mul(length - 1)
                    .saturating_add(
                        diagonal.saturating_mul((length - 1).saturating_mul(length - 1)),
                    )
                    .saturating_mul(2);
                let cost = length.saturating_mul(degree).saturating_add(quadratic);
                let candidate = (
                    self.bounds[j] != Bounds::FREE,
                    if options.cost_aware { cost } else { degree },
                    j,
                    degree,
                    cost,
                );
                if options.max_pivot_attempts == 1 && !candidates.is_empty() {
                    if candidate < candidates[0] {
                        candidates[0] = candidate;
                    }
                } else {
                    candidates.push(candidate);
                }
            }
            if candidates.len() > options.max_pivot_attempts {
                candidates.select_nth_unstable(options.max_pivot_attempts);
                candidates.truncate(options.max_pivot_attempts);
            }
            candidates.sort_unstable();
            let mut deferred = false;
            for &(_, _, j, degree, estimate) in candidates.iter().take(options.max_pivot_attempts) {
                if Instant::now() >= deadline {
                    deferred = true;
                    break;
                }
                let cost = self.a.column(j).iter().fold(estimate, |cost, (k, _)| {
                    cost.saturating_add(self.a.row(k).len())
                });
                if cost > work {
                    self.equality_stats.work_rejections += 1;
                    deferred = true;
                    continue;
                }
                work -= cost;
                self.equality_stats.estimated_work =
                    self.equality_stats.estimated_work.saturating_add(cost);
                let effective = if self.bounds[j] == Bounds::FREE {
                    Bounds::FREE
                } else {
                    let implied = self.singleton_range(i, j, self.a.get(i, j), bounds.lower);
                    self.non_implied_bounds(j, implied)
                };
                let fill = if options.preserve_nonzeros {
                    let removed = degree.saturating_add(if effective == Bounds::FREE {
                        length - 1
                    } else {
                        0
                    });
                    max_fill.min(removed)
                } else {
                    max_fill
                };
                self.equality_stats.attempts += 1;
                if self.substitute(i, j, effective, fill) {
                    self.equality_stats.accepted += 1;
                    deferred = false;
                    break;
                }
                self.equality_stats.rejected_updates += 1;
                use crate::core::objective::SubstitutionFailure;
                match self.substitution_failure {
                    SubstitutionFailure::Numerical => self.equality_stats.numerical_rejections += 1,
                    SubstitutionFailure::ConstraintFill => {
                        self.equality_stats.constraint_fill_rejections += 1
                    }
                    SubstitutionFailure::QuadraticFill => {
                        self.equality_stats.quadratic_fill_rejections += 1
                    }
                    SubstitutionFailure::HessianGrowth => {
                        self.equality_stats.hessian_growth_rejections += 1
                    }
                    SubstitutionFailure::Deadline => self.equality_stats.deadline_rejections += 1,
                }
            }
            if deferred && !matches!(options.work_limit, crate::settings::WorkLimit::Default) {
                self.queues.short_equalities.push(i);
            }
        }
    }

    /// Bounds of the singleton value when its row is made tight at `rhs`.
    /// These prove which original bounds can be removed with the column.
    fn singleton_range(&mut self, row: usize, column: usize, pivot: f64, rhs: f64) -> Bounds {
        let residual = self.residual_activity(row, column);
        let (lower, upper) = if pivot > 0.0 {
            (residual.max.value(), residual.min.value())
        } else {
            (residual.min.value(), residual.max.value())
        };
        Bounds {
            lower: lower.map_or(f64::NEG_INFINITY, |v| (rhs - v) / pivot),
            upper: upper.map_or(f64::INFINITY, |v| (rhs - v) / pivot),
        }
    }

    fn non_implied_bounds(&self, column: usize, implied: Bounds) -> Bounds {
        let b = self.bounds[column];
        Bounds {
            lower: if implied.lower.is_finite() && implied.lower >= b.lower {
                f64::NEG_INFINITY
            } else {
                b.lower
            },
            upper: if implied.upper.is_finite() && implied.upper <= b.upper {
                f64::INFINITY
            } else {
                b.upper
            },
        }
    }

    pub fn singleton_columns(&mut self, max_fill: usize) {
        // A neighbour's bound may make this column implied free without
        // changing its own degree. Inspect each affected row once per round.
        for i in self.queues.singleton_activity_rows.take_round() {
            if self.rows[i] == RowDomain::Deleted {
                continue;
            }
            for (j, _) in self.a.row(i) {
                if self.a.column(j).len() == 1 {
                    self.queues.singleton_columns.push(j);
                }
            }
        }
        for j in self.queues.singleton_columns.take_round() {
            if self.deadline.is_some_and(|d| Instant::now() >= d) {
                break;
            }
            if !self.alive[j] {
                continue;
            }
            let column = self.a.column(j);
            if column.len() != 1 {
                continue;
            };
            let (i, a) = column.iter().next().unwrap();
            let RowDomain::Linear(b) = self.rows[i] else {
                continue;
            };
            if self.a.row(i).len() <= 1 {
                continue;
            }
            if b.equality() {
                // A small pivot can repeatedly amplify curvature and create
                // a huge objective constant that later cancels at the solution.
                // Defer curved singletons to a better-scaled doubleton pivot.
                if !self.objective.p.column(j).is_empty()
                    && self.a.row(i).iter().any(|(_, v)| v.abs() > a.abs())
                {
                    continue;
                }
                let implied = self.singleton_range(i, j, a, b.lower);
                let effective = self.non_implied_bounds(j, implied);
                // Defer columns with two non-implied bounds to the doubleton rule.
                if effective.lower.is_finite() && effective.upper.is_finite() {
                    continue;
                }
                self.substitute(i, j, effective, max_fill);
                continue;
            }
            // A linear singleton prefers a particular side of its row. A
            // coupled or curved objective may prefer an interior point, so
            // the LP inequality-to-equality rule does not apply to it.
            if !self.objective.p.column(j).is_empty() {
                continue;
            }
            let c = self.objective.c[j];
            let rhs = if c / a > 0.0 {
                b.lower
            } else if c / a < 0.0 {
                b.upper
            } else if b.lower.is_finite() {
                b.lower
            } else {
                b.upper
            };
            // Infinite improving sides are handled by simple_dual_fix, which
            // also constructs the recession certificate.
            if !rhs.is_finite() {
                continue;
            }
            let implied = self.singleton_range(i, j, a, rhs);
            let effective = self.non_implied_bounds(j, implied);
            let reachable = if c > 0.0 {
                !effective.lower.is_finite()
            } else if c < 0.0 {
                !effective.upper.is_finite()
            } else {
                effective == Bounds::FREE
            };
            if reachable {
                self.substitute_on_side(i, j, rhs, effective, max_fill);
            }
        }
    }

    pub fn doubleton_equalities(&mut self, max_fill: usize) {
        for i in self.queues.doubleton_rows.take_round() {
            if self.deadline.is_some_and(|d| Instant::now() >= d) {
                break;
            }
            if !matches!(self.rows[i],RowDomain::Linear(b) if b.equality()) {
                continue;
            }
            let row = self.a.row(i);
            if row.len() != 2 {
                continue;
            };
            let mut entries = row.iter();
            let (j, a) = entries.next().unwrap();
            let (k, b) = entries.next().unwrap();
            let nj = self.a.column(j).len();
            let nk = self.a.column(k).len();
            // Prefer singleton pivots, then integral substitution
            // ratio, then the shorter column. QP fill is checked afterwards.
            let integer = |v: f64| v.is_finite() && (v - v.round()).abs() <= 1e-12 * v.abs();
            let quadratic =
                !self.objective.p.column(j).is_empty() || !self.objective.p.column(k).is_empty();
            let remove_j = if quadratic && a.abs() != b.abs() {
                // LP's integral-ratio preference may repeatedly square an
                // amplifying factor in P. Use the larger pivot for QPs.
                a.abs() > b.abs()
            } else if nj == 1 && nk != 1 {
                true
            } else if nk == 1 && nj != 1 {
                false
            } else if integer(b / a) != integer(a / b) {
                integer(b / a)
            } else {
                nj < nk
            };
            let (first, second, ratio) = if remove_j {
                (j, k, b / a)
            } else {
                (k, j, a / b)
            };
            // Rust has no CSR shift limit; max_fill bounds newly allocated
            // coefficients in A and P; Hessian growth has a separate policy. Try the
            // other pivot if fill or arithmetic makes the preferred direction unsuitable.
            if !ratio.is_finite() || !(1e-7..=1e7).contains(&ratio.abs()) {
                continue;
            }
            if !self.substitute(i, first, self.bounds[first], max_fill) && !quadratic {
                self.substitute(i, second, self.bounds[second], max_fill);
            }
        }
    }
}
