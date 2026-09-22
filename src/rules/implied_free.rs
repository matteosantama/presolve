//! Equality aggregation whose pivot bounds are implied by that same equation.
//!
//! Removing a bounded pivot normally leaves its bounds as a new row. Here the
//! activity proof makes that row redundant, so each accepted pivot removes one
//! row as well as one column. Degree-two columns remain attractive in long rows:
//! their worst-case net constraint fill is `(r - 2) * (c - 2) - 2 = -2`.
//! The shared equality settings govern pivot stability, attempts, work, minimum
//! column degree, curvature, cost ranking, and fill preservation. Degree-two
//! pivots bypass maximum row/column length; explicit freedom is replaced by the
//! same-equation implication proof, independent of `require_free_variable`.

use crate::{
    model::{Model, RowDomain, activity::Activity, queues::Worklist},
    problem::Bounds,
};
use std::{cell::Cell, time::Instant};

impl Model {
    pub(super) fn implied_free_equalities(&mut self, deadline: Instant) -> usize {
        let options = self.settings.equalities;
        let relative = if options.relative_pivot > 0.0 && options.relative_pivot <= 1.0 {
            options.relative_pivot
        } else {
            1.0
        };
        let mut work = options.work_limit.resolve(self.a.nnz().saturating_mul(16));
        if work == 0 || options.max_pivot_attempts == 0 || Instant::now() >= deadline {
            return 0;
        }
        let mut sparse = Worklist::new(self.rows.len());
        let mut general = Worklist::new(self.rows.len());
        for (i, domain) in self.rows.iter().enumerate() {
            if i % 256 == 0 && Instant::now() >= deadline {
                return 0;
            }
            if !matches!(domain, RowDomain::Linear(b) if b.equality()) {
                continue;
            }
            let row = self.a.row(i);
            if row.len() < 2 || row.len() > work {
                continue;
            }
            work -= row.len();
            let sparse_pivot = row.iter().any(|(j, _)| self.a.column(j).len() <= 2);
            if sparse_pivot {
                sparse.push(i);
            } else {
                general.push(i);
            }
        }
        // Two degree buckets avoid a global sort. Affected rows are retried
        // locally; the shared matrix mutation invalidates their activity cache.
        let mut accepted = 0;
        let mut candidates = Vec::new();
        let mut affected = Vec::new();
        while let Some(i) = sparse.pop().or_else(|| general.pop()) {
            if work == 0 || Instant::now() >= deadline {
                break;
            }
            let RowDomain::Linear(bounds) = self.rows[i] else {
                continue;
            };
            if !bounds.equality() {
                continue;
            }
            let length = self.a.row(i).len();
            let scan_cost = length.saturating_mul(3);
            if length < 2 || scan_cost > work {
                continue;
            }
            work -= scan_cost;
            let activity = self.activity(i);
            let row = self.a.row(i);
            let largest = row.iter().map(|(_, a)| a.abs()).fold(0.0, f64::max);
            candidates.clear();
            for (j, a) in row {
                let degree = self.a.column(j).len();
                if degree < options.min_column_length {
                    continue;
                }
                if degree > 2
                    && (row.len() > options.max_row_length || degree > options.max_column_length)
                {
                    continue;
                }
                if options.require_linear_variable && !self.objective.p.column(j).is_empty() {
                    continue;
                }
                if a.abs() < largest && (relative == 1.0 || a.abs() / largest < relative) {
                    continue;
                }
                // Cached finite sums and infinity counts give constant-time
                // residuals. Charge the rare cancellation fallback explicitly.
                let remaining = Cell::new(work);
                let (lower, upper) = activity.implied(a, self.bounds[j], bounds, || {
                    if remaining.get() < length || Instant::now() >= deadline {
                        Activity::UNKNOWN
                    } else {
                        remaining.set(remaining.get() - length);
                        self.residual_activity(i, j)
                    }
                });
                work = remaining.get();
                let (lower, upper) = if a > 0.0 {
                    (lower, upper)
                } else {
                    (upper, lower)
                };
                let implied = Bounds {
                    lower: lower.unwrap_or(f64::NEG_INFINITY),
                    upper: upper.unwrap_or(f64::INFINITY),
                };
                if self.non_implied_bounds(j, implied) != Bounds::FREE {
                    continue;
                }
                let diagonal = usize::from(self.objective.p.get(j, j) != 0.0);
                let quadratic = (self.objective.p.column(j).len() - diagonal)
                    .saturating_mul(length - 1)
                    .saturating_add(
                        diagonal.saturating_mul((length - 1).saturating_mul(length - 1)),
                    )
                    .saturating_mul(2);
                let estimate = degree.saturating_mul(length).saturating_add(quadratic);
                let rank = if options.cost_aware { estimate } else { degree };
                candidates.push((rank, degree, j, estimate));
            }
            if candidates.len() > options.max_pivot_attempts {
                candidates.select_nth_unstable(options.max_pivot_attempts);
                candidates.truncate(options.max_pivot_attempts);
            }
            candidates.sort_unstable();
            for &(_, degree, j, estimate) in &candidates {
                let cost = self.a.column(j).iter().fold(estimate, |cost, (k, _)| {
                    cost.saturating_add(self.a.row(k).len())
                });
                if cost > work || Instant::now() >= deadline {
                    continue;
                }
                work -= cost;
                affected.clear();
                affected.extend(self.a.column(j).iter().map(|(k, _)| k));
                let fill = if options.preserve_nonzeros {
                    self.settings
                        .substitution_fill
                        .min(degree.saturating_add(length - 1))
                } else {
                    self.settings.substitution_fill
                };
                if self.substitute(i, j, Bounds::FREE, fill) {
                    accepted += 1;
                    for &k in &affected {
                        if k != i {
                            sparse.push(k);
                        }
                    }
                    break;
                }
            }
        }
        accepted
    }
}
