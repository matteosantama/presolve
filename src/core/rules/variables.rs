// SPDX-License-Identifier: Apache-2.0
// Modified for this library; copyright and attribution notices are in NOTICE.
//! Fix known variables and eliminate independent empty columns.

use crate::{core::model::Model, postsolve::tape::Certificate, problem::Bounds};

/// Beyond this Hessian degree the Schur complement is dense enough that the
/// no-growth check rejects the elimination anyway.
const MAX_ELIMINATION_DEGREE: usize = 16;

impl Model {
    pub fn fixed_variables(&mut self) {
        while let Some(j) = self.queues.fixed_columns.pop() {
            if self.alive[j] && self.bounds[j].equality() {
                self.fix(j, self.bounds[j].lower);
            }
        }
    }

    pub fn empty_columns(&mut self, max_fill: usize) -> Result<(), Certificate> {
        while let Some(j) = self.queues.empty_columns.pop() {
            if !self.alive[j] || !self.a.column(j).is_empty() {
                continue;
            }
            let Some(p) = self.objective.diagonal(j) else {
                // A coupled free column has an interior minimizer in closed
                // form; a bounded one would need a clipped, nonaffine rule.
                // A rejected elimination is retried only after the Hessian
                // changed, since the column's fill is a function of `P` alone.
                if self.rules.quadratic_elimination
                    && self.bounds[j] == Bounds::FREE
                    && self.objective.p.row(j).len() <= MAX_ELIMINATION_DEGREE
                    && self.elimination_rejected[j] != self.objective.p.revision + 1
                    && !self.eliminate_coupled(j, max_fill)
                {
                    self.elimination_rejected[j] = self.objective.p.revision + 1;
                }
                continue;
            };
            if p < 0.0 {
                continue;
            }
            let c = self.objective.c[j];
            let b = self.bounds[j];
            let value = if p > 0.0 {
                (-c / p).max(b.lower).min(b.upper)
            } else if c > 0.0 {
                b.lower
            } else if c < 0.0 {
                b.upper
            } else {
                0.0_f64.max(b.lower).min(b.upper)
            };
            if value.is_finite() {
                self.fix(j, value);
            } else if p == 0.0 && c != 0.0 {
                return Err(self.recession_certificate(j, -c.signum()));
            }
        }
        Ok(())
    }

    pub(super) fn recession_certificate(&self, column: usize, direction: f64) -> Certificate {
        self.dual_certificate([(column, direction)])
    }
}
