// SPDX-License-Identifier: Apache-2.0
// Modified for this library; copyright and attribution notices are in NOTICE.
//! Use objective derivatives and constraint locks to fix variables.

use crate::{
    core::model::Model,
    postsolve::tape::{Certificate, Side},
};

impl Model {
    pub fn simple_dual_fix(&mut self) -> Result<(), Certificate> {
        // Newly unlocked neighbours can appear while processing this round.
        // They remain queued for the next cleanup cycle.
        for j in self.queues.unlocked_columns.take_round() {
            if !self.alive[j] {
                continue;
            }
            let lock = self.locks[j];
            if lock.up > 0 && lock.down > 0 {
                continue;
            }
            let Some(p) = self.objective.diagonal(j) else {
                continue;
            };
            if p < 0.0 {
                continue;
            }
            let c = self.objective.c[j];
            let b = self.bounds[j];
            if p == 0.0 && c == 0.0 {
                if lock.down == 0 {
                    if b.lower.is_finite() {
                        self.fix(j, b.lower);
                    } else {
                        self.remove_unlocked(j);
                    }
                } else if lock.up == 0 {
                    if b.upper.is_finite() {
                        self.fix(j, b.upper);
                    } else {
                        self.remove_unlocked(j);
                    }
                }
                continue;
            }
            if p == 0.0
                && ((c > 0.0 && lock.down == 0 && !b.lower.is_finite())
                    || (c < 0.0 && lock.up == 0 && !b.upper.is_finite()))
            {
                return Err(self.recession_certificate(j, -c.signum()));
            }
            let gradient = |x: f64| p * x + c;
            if lock.down == 0
                && b.lower.is_finite()
                && gradient(b.lower).is_finite()
                && gradient(b.lower) >= 0.0
            {
                self.fix(j, b.lower);
            } else if lock.up == 0
                && b.upper.is_finite()
                && gradient(b.upper).is_finite()
                && gradient(b.upper) <= 0.0
            {
                self.fix(j, b.upper);
            }
        }
        Ok(())
    }

    /// Bound the coupled derivative only for unlocked, sparse Hessian rows.
    /// One scan per medium phase avoids tracking every Hessian neighbour on
    /// every bound change, which would be expensive for dense objectives.
    pub fn coupled_dual_fix(&mut self) {
        if self.objective.p.nnz() == 0 {
            return;
        }
        for j in 0..self.alive.len() {
            if !self.alive[j] || !(2..=33).contains(&self.objective.p.row(j).len()) {
                continue;
            }
            let lock = self.locks[j];
            let diagonal = self.objective.p.get(j, j);
            if (lock.down > 0 && lock.up > 0) || diagonal <= 0.0 {
                continue;
            }
            for side in [Side::Lower, Side::Upper] {
                let value = side.value(self.bounds[j]);
                if !value.is_finite()
                    || (side == Side::Lower && lock.down > 0)
                    || (side == Side::Upper && lock.up > 0)
                {
                    continue;
                }
                let mut gradient = self.objective.c[j];
                let mut magnitude = gradient.abs();
                for &(k, a) in self.objective.p.row(j) {
                    let x = if k == j {
                        value
                    } else if (side == Side::Lower) == (a > 0.0) {
                        self.bounds[k].lower
                    } else {
                        self.bounds[k].upper
                    };
                    let term = a * x;
                    gradient += term;
                    magnitude += term.abs();
                }
                // Require a resolved sign rather than rounding cancellation
                // into an optimality proof. Unbounded terms fail this check.
                let margin = 64.0 * f64::EPSILON * magnitude;
                let improving = match side {
                    Side::Lower => gradient > margin,
                    Side::Upper => gradient < -margin,
                };
                if gradient.is_finite() && improving && self.fix(j, value) {
                    break;
                }
            }
        }
    }
}
