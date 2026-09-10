// SPDX-License-Identifier: Apache-2.0
// Modified for this library; copyright and attribution notices are in NOTICE.
//! Propagate implied bounds and remove redundant variable bounds.

use crate::{
    core::{
        activity::Activity,
        model::{Model, RowDomain},
    },
    postsolve::tape::{Certificate, Equation, Side},
    problem::Bounds,
};
use std::sync::Arc;

impl Model {
    /// One medium-phase propagation round. Each bound proof owns the equation
    /// used to derive it; follow-up structural changes cannot invalidate tape.
    pub fn propagate_bounds(&mut self) -> Result<usize, Certificate> {
        let mut tightened = 0;
        for i in self.queues.changed_activities.take_round() {
            let RowDomain::Linear(original_bounds) = self.rows[i] else {
                continue;
            };
            if self.a.row(i).len() <= 1 {
                continue;
            }
            let activity = self.activity(i);
            let mut bounds = original_bounds;
            if bounds.upper.is_finite()
                && activity
                    .min
                    .value()
                    .is_some_and(|v| self.separated(v, bounds.upper))
            {
                return Err(self.row_certificate(i, -1.0));
            }
            if bounds.lower.is_finite()
                && activity
                    .max
                    .value()
                    .is_some_and(|v| self.separated(bounds.lower, v))
            {
                return Err(self.row_certificate(i, 1.0));
            }
            // Redundancy requires actual containment; a small interval may
            // carry a large quadratic objective cost.
            if bounds.lower.is_finite() && activity.min.value().is_some_and(|v| v >= bounds.lower) {
                bounds.lower = f64::NEG_INFINITY;
            }
            if bounds.upper.is_finite() && activity.max.value().is_some_and(|v| v <= bounds.upper) {
                bounds.upper = f64::INFINITY;
            }
            if bounds == Bounds::FREE {
                self.delete_row(i);
                continue;
            }
            if bounds != original_bounds {
                self.relax_row(i, bounds);
            }
            // Removing one term cannot make an extreme with two unbounded
            // contributions finite. Keep the row checks above this shortcut.
            if (!bounds.lower.is_finite() || activity.max.infinite >= 2)
                && (!bounds.upper.is_finite() || activity.min.infinite >= 2)
            {
                continue;
            }
            // Bound changes leave coefficients intact. Snapshot them only if
            // a tightening needs a proof, and share it for the rest of this row.
            let mut equation = None;
            let mut cursor = self.a.row(i).cursor();
            while let Some((j, a)) = cursor.next(&self.a) {
                let act = self.activity(i);
                let (min_term, max_term) = Activity::terms(a, self.bounds[j]);
                let lower = if bounds.lower.is_finite() {
                    act.max
                        .excluding(max_term)
                        .or_else(|| {
                            // Any other infinite contribution prevents a finite bound.
                            (act.max.infinite == usize::from(!max_term.is_finite()))
                                .then(|| self.residual_activity(i, j).max.value())
                                .flatten()
                        })
                        .map(|v| (bounds.lower - v) / a)
                } else {
                    None
                };
                let upper = if bounds.upper.is_finite() {
                    act.min
                        .excluding(min_term)
                        .or_else(|| {
                            (act.min.infinite == usize::from(!min_term.is_finite()))
                                .then(|| self.residual_activity(i, j).min.value())
                                .flatten()
                        })
                        .map(|v| (bounds.upper - v) / a)
                } else {
                    None
                };
                if let Some(value) = lower {
                    tightened += usize::from(self.implied_bound(
                        j,
                        if a > 0.0 { Side::Lower } else { Side::Upper },
                        value,
                        |model| {
                            equation
                                .get_or_insert_with(|| model.equation(i).unwrap())
                                .clone()
                        },
                        true,
                    )?);
                }
                if let Some(value) = upper {
                    tightened += usize::from(self.implied_bound(
                        j,
                        if a > 0.0 { Side::Upper } else { Side::Lower },
                        value,
                        |model| {
                            equation
                                .get_or_insert_with(|| model.equation(i).unwrap())
                                .clone()
                        },
                        true,
                    )?);
                }
            }
        }
        Ok(tightened)
    }

    /// Remove redundant bounds after the main rule phases,
    /// so implied bounds help rules without becoming extra barrier terms.
    /// Recompute implications as bounds are removed to avoid circular proofs.
    pub fn remove_redundant_bounds(&mut self) {
        let mut columns: Vec<_> = (0..self.bounds.len())
            .filter(|&j| {
                self.alive[j]
                    && !self.a.column(j).is_empty()
                    && (self.bounds[j].lower.is_finite() || self.bounds[j].upper.is_finite())
            })
            .collect();
        columns.sort_unstable_by_key(|&j| (self.a.column(j).len(), j));
        let mut lower_candidates = Vec::new();
        for j in columns {
            let b = self.bounds[j];
            // Preserve the existing one-sided sweep before trying extra sides.
            // An early lower-bound removal can invalidate a later upper-bound
            // proof, leaving more bound rows than the original cleanup did.
            if b.upper.is_finite() {
                if b.lower.is_finite() {
                    lower_candidates.push(j);
                }
                self.remove_redundant_bound(j, Side::Upper);
            } else if b.lower.is_finite() {
                self.remove_redundant_bound(j, Side::Lower);
            }
        }
        for j in lower_candidates {
            self.remove_redundant_bound(j, Side::Lower);
        }
    }

    fn remove_redundant_bound(&mut self, j: usize, side: Side) {
        let b = self.bounds[j];
        let mut cursor = self.a.column(j).cursor();
        while let Some((i, a)) = cursor.next(&self.a) {
            let RowDomain::Linear(row) = self.rows[i] else {
                continue;
            };
            let act = self.activity(i);
            let (min, max) = Activity::terms(a, b);
            let from_lower = (side == Side::Lower) == (a > 0.0);
            let (rhs, extreme, term) = if from_lower {
                (row.lower, act.max, max)
            } else {
                (row.upper, act.min, min)
            };
            if !rhs.is_finite() || extreme.infinite != usize::from(!term.is_finite()) {
                continue;
            }
            let residual = extreme.excluding(term).or_else(|| {
                let residual = self.residual_activity(i, j);
                if from_lower {
                    residual.max.value()
                } else {
                    residual.min.value()
                }
            });
            let Some(bound) = residual.map(|v| (rhs - v) / a).filter(|v| v.is_finite()) else {
                continue;
            };
            let implied = match side {
                Side::Lower => bound >= b.lower,
                Side::Upper => bound <= b.upper,
            };
            if implied {
                self.relax_bound(j, side);
                break;
            }
        }
    }

    pub(super) fn implied_bound(
        &mut self,
        j: usize,
        side: Side,
        value: f64,
        proof: impl FnOnce(&Self) -> Arc<Equation>,
        propagation: bool,
    ) -> Result<bool, Certificate> {
        if !value.is_finite() || (propagation && value.abs() >= self.numerics.huge_bound) {
            return Ok(false);
        }
        let old = self.bounds[j];
        let opposite = match side {
            Side::Lower => old.upper,
            Side::Upper => old.lower,
        };
        let (lower, upper) = match side {
            Side::Lower => (value, opposite),
            Side::Upper => (opposite, value),
        };
        if lower > upper {
            if self.separated(lower, upper) {
                let equation = proof(self);
                let a = self.a.get(equation.row, j);
                let sign = match side {
                    Side::Lower => a.signum(),
                    Side::Upper => -a.signum(),
                };
                return Err(self.row_certificate(equation.row, sign));
            }
            return Ok(false);
        }
        if propagation && side.value(old).is_finite() && value != opposite {
            // Skip insignificant finite changes using a relative threshold
            // and a floor scaled by the configured feasibility tolerance.
            let gain = match side {
                Side::Lower => value - old.lower,
                Side::Upper => old.upper - value,
            };
            if gain
                <= (self.propagation.minimum_gain_factor * self.numerics.feasibility)
                    .max(self.propagation.minimum_relative_gain * side.value(old).abs())
            {
                return Ok(false);
            }
        }
        if (side == Side::Lower && value <= old.lower)
            || (side == Side::Upper && value >= old.upper)
        {
            return Ok(false);
        }
        let equation = proof(self);
        Ok(self.tighten_bound(j, side, value, equation))
    }

    pub(super) fn separated(&self, lower: f64, upper: f64) -> bool {
        lower > upper + self.numerics.feasibility * (1.0 + lower.abs().max(upper.abs()))
    }
}
