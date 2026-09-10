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
    /// Equality substitution with configurable structural and work limits.
    /// Retain the largest-magnitude pivot requirement to avoid amplification.
    pub fn short_equalities(&mut self, max_fill: usize, deadline: Instant) {
        let options = self.equalities;
        let mut work = options.work_limit.resolve(self.a.nnz().saturating_mul(2));
        for i in self.queues.short_equalities.take_round() {
            if work == 0 || Instant::now() >= deadline {
                if matches!(options.work_limit, crate::settings::WorkLimit::Default) {
                    break;
                }
                // A later cycle may have a fresh work allowance.
                self.queues.short_equalities.push(i);
                continue;
            }
            let RowDomain::Linear(bounds) = self.rows[i] else {
                continue;
            };
            if !bounds.equality() {
                continue;
            }
            let row = self.a.row(i);
            let length = row.len();
            if length < 3 || length > options.max_row_length {
                continue;
            }
            let largest = row.iter().map(|(_, a)| a.abs()).fold(0.0, f64::max);
            let pivot = row
                .iter()
                .filter_map(|(j, a)| {
                    let degree = self.a.column(j).len();
                    ((!options.require_free_variable || self.bounds[j] == Bounds::FREE)
                        && (!options.require_linear_variable
                            || self.objective.p.column(j).is_empty())
                        && degree >= options.min_column_length
                        && degree <= options.max_column_length
                        && a.abs() == largest)
                        .then_some((self.bounds[j] != Bounds::FREE, degree, j))
                })
                .min();
            let Some((_, degree, j)) = pivot else {
                continue;
            };
            let cost = self
                .a
                .column(j)
                .iter()
                .fold(length.saturating_mul(degree), |cost, (k, _)| {
                    cost.saturating_add(self.a.row(k).len())
                });
            if cost > work {
                if !matches!(options.work_limit, crate::settings::WorkLimit::Default) {
                    self.queues.short_equalities.push(i);
                }
                continue;
            }
            work -= cost;
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
            self.substitute(i, j, effective, fill);
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
