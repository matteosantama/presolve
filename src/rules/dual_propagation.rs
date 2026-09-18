// SPDX-License-Identifier: Apache-2.0
// Modified for this library; copyright and attribution notices are in NOTICE.
//! Propagate the dual constraints to find shift directions. A multiplier
//! bound proved by propagation is a conic combination of dual rows, and that
//! combination is a primal direction: moving along it keeps every row and
//! bound feasible and never increases the objective. When the direction also
//! drives a variable to its bound, or a row to one side, any feasible point
//! can slide there, so the reduction needs no optimal solution to exist. Every
//! direction is verified exactly on the current model before it is used.

use crate::{
    model::{
        Model, RowDomain,
        activity::Activity,
        queues::Worklist,
        tape::{Certificate, Side},
    },
    problem::Bounds,
};
use std::collections::BTreeMap;

/// One multiplier tightening: a bound of row `row` was derived from the dual
/// row of `column`, on the `constraint` side of its interval.
#[derive(Clone, Copy, Debug)]
struct Record {
    row: u32,
    column: u32,
    constraint: Side,
}

/// Scratch space reused across passes so a pass that proves nothing allocates nothing.
#[derive(Debug, Default)]
pub(crate) struct DualScratch {
    /// Multiplier interval per row, in native signs.
    y: Vec<Bounds>,
    queue: Worklist,
    /// Rows whose multiplier tightened in the current round. Their columns
    /// join the next round once per row rather than once per tightening.
    rows: Worklist,
    round: Vec<usize>,
    tightened: Vec<usize>,
    /// Column activity over `y` from the column's latest visit. Every
    /// tightening revisits the affected columns, so after a drained queue
    /// each stored activity reflects the final multipliers.
    column_activity: Vec<Option<Activity>>,
    /// Proof records in derivation order; sequence number is index plus one.
    records: Vec<Record>,
    /// Ascending sequence numbers of the records for each `(row, side)`.
    proofs: Vec<Vec<u32>>,
    /// Direction extraction: net dual-row multipliers and the accumulated
    /// row activity of the direction, both dense with touched-index lists.
    lambda: Vec<f64>,
    lambda_touched: Vec<usize>,
    ad: Vec<f64>,
    ad_touched: Vec<usize>,
    weights: BTreeMap<u32, f64>,
}

/// Allowed interval of a linear row's multiplier from its finite sides.
fn multiplier_range(bounds: Bounds) -> Bounds {
    Bounds {
        lower: if bounds.lower.is_finite() && !bounds.upper.is_finite() {
            0.0
        } else {
            f64::NEG_INFINITY
        },
        upper: if bounds.upper.is_finite() && !bounds.lower.is_finite() {
            0.0
        } else {
            f64::INFINITY
        },
    }
}

fn slot(row: usize, side: Side) -> usize {
    2 * row + usize::from(side == Side::Upper)
}

/// The bound of `y_k` that extremizes `a * y_k` on the given side of a dual row:
/// the `Upper` constraint bounds the sum from above, so each term takes its
/// minimum; the `Lower` constraint takes each term's maximum.
fn used_bound(constraint: Side, a: f64) -> Side {
    match (constraint, a > 0.0) {
        (Side::Upper, true) | (Side::Lower, false) => Side::Lower,
        (Side::Upper, false) | (Side::Lower, true) => Side::Upper,
    }
}

/// What a verified direction does.
#[derive(Clone, Copy)]
enum Shift {
    /// Slides `column` onto this bound.
    Fix(usize, Side),
    /// Slides the activity of `row` onto this side.
    Restrict(usize, Side),
}

impl Model {
    /// Interval for `Σ a_ij y_i` implied by the sign of the reduced cost.
    fn dual_row_range(&self, b: Bounds, j: usize) -> Option<Bounds> {
        let c = self.objective.c[j];
        match (b.lower.is_finite(), b.upper.is_finite()) {
            (true, false) => Some(Bounds {
                lower: f64::NEG_INFINITY,
                upper: c,
            }),
            (false, true) => Some(Bounds {
                lower: c,
                upper: f64::INFINITY,
            }),
            (false, false) => Some(Bounds::fixed(c)),
            (true, true) => None,
        }
    }

    /// One pass: derive multiplier intervals with proofs, then turn each
    /// conclusion into a direction and apply the ones that verify. Returns the
    /// number of rows restricted and columns fixed. The pass runs after
    /// redundant variable bounds are removed, so every remaining finite side
    /// genuinely constrains the dual; an implied side would only weaken a
    /// column's dual row.
    pub fn dual_propagation(&mut self) -> Result<usize, Certificate> {
        let m = self.rows.len();
        let n = self.bounds.len();
        // A problem whose live columns all carry curvature has no dual row.
        if !(0..n).any(|j| self.linear_column(j) && !self.a.column(j).is_empty()) {
            return Ok(0);
        }
        let mut scratch = std::mem::take(&mut self.dual_scratch);
        scratch.y.clear();
        scratch.y.extend(self.rows.iter().map(|row| match row {
            RowDomain::Linear(b) => multiplier_range(*b),
            _ => Bounds::FREE,
        }));
        scratch.queue.reset(n);
        scratch.rows.reset(m);
        scratch.column_activity.clear();
        scratch.column_activity.resize(n, None);
        scratch.records.clear();
        scratch.proofs.clear();
        scratch.proofs.resize(2 * m, Vec::new());
        let budget = self
            .settings
            .dual_propagation
            .work_limit
            .resolve(self.a.nnz().saturating_mul(4));
        let mut work = budget;
        for j in 0..n {
            if self.linear_column(j) && self.dual_row_range(self.bounds[j], j).is_some() {
                scratch.queue.push(j);
            }
        }
        let mut consistent = true;
        // A drained queue means every column was visited after its last
        // multiplier change; the work limit leaves stored activities stale.
        let mut complete = true;
        // Round-based processing: a column changed during a round is revisited
        // in the next one, so long implication chains cost one visit per round
        // instead of one per intermediate improvement. Columns of a tightened
        // row are queued when the round ends, in the order the rows first
        // tightened, which matches queuing them at each tightening.
        'propagate: loop {
            scratch.rows.swap_round(&mut scratch.tightened);
            for &i in &scratch.tightened {
                for (k, _) in self.a.row(i) {
                    scratch.queue.push(k);
                }
            }
            scratch.queue.swap_round(&mut scratch.round);
            if scratch.round.is_empty() {
                break;
            }
            for at in 0..scratch.round.len() {
                let j = scratch.round[at];
                let column = self.a.column(j);
                if column.len() > work {
                    complete = false;
                    break 'propagate;
                }
                work -= column.len();
                let Some(range) = self.dual_row_range(self.bounds[j], j) else {
                    continue;
                };
                if !self.linear_column(j) {
                    continue;
                }
                let activity = Activity::compute(column, &scratch.y, None);
                scratch.column_activity[j] = Some(activity);
                if (range.upper.is_finite()
                    && activity
                        .min
                        .value()
                        .is_some_and(|v| self.separated(v, range.upper)))
                    || (range.lower.is_finite()
                        && activity
                            .max
                            .value()
                            .is_some_and(|v| self.separated(range.lower, v)))
                {
                    consistent = false;
                    break 'propagate;
                }
                // Removing one term cannot make an extreme with two unbounded
                // contributions finite.
                if (!range.lower.is_finite() || activity.max.infinite >= 2)
                    && (!range.upper.is_finite() || activity.min.infinite >= 2)
                {
                    continue;
                }
                for (i, a) in column {
                    let (from_lower, from_upper) = activity.implied(a, scratch.y[i], range, || {
                        Activity::compute(column, &scratch.y, Some(i))
                    });
                    // A positive coefficient turns the lower constraint into a
                    // lower bound; a negative one flips the sides.
                    let (lower, upper) = if a > 0.0 {
                        ((from_lower, Side::Lower), (from_upper, Side::Upper))
                    } else {
                        ((from_upper, Side::Upper), (from_lower, Side::Lower))
                    };
                    for (side, (value, constraint)) in [(Side::Lower, lower), (Side::Upper, upper)]
                    {
                        let Some(value) = value else { continue };
                        match self.tighten_multiplier(&mut scratch.y[i], side, value) {
                            Some(true) => {
                                scratch.rows.push(i);
                                scratch.records.push(Record {
                                    row: i as u32,
                                    column: j as u32,
                                    constraint,
                                });
                                let sequence = scratch.records.len() as u32;
                                scratch.proofs[slot(i, side)].push(sequence);
                            }
                            Some(false) => {}
                            None => {
                                consistent = false;
                                break 'propagate;
                            }
                        }
                    }
                }
            }
        }
        scratch.queue.clear();
        scratch.rows.clear();
        let mut reductions = 0;
        let result = (|| {
            if !consistent {
                return Ok(());
            }
            // Extraction shares the pass allowance so long proofs cannot make
            // the conclusions cost more than the propagation did.
            let mut extraction = budget;
            for i in 0..m {
                let RowDomain::Linear(b) = self.rows[i] else {
                    continue;
                };
                if b.equality() || self.a.row(i).is_empty() {
                    continue;
                }
                let y = scratch.y[i];
                let side = if b.lower.is_finite() && self.separated(y.lower, 0.0) {
                    Side::Lower
                } else if b.upper.is_finite() && self.separated(0.0, y.upper) {
                    Side::Upper
                } else {
                    continue;
                };
                scratch.weights.clear();
                if let Some(&sequence) = scratch.proofs[slot(i, side)].last() {
                    scratch.weights.insert(sequence, 1.0);
                }
                if self.apply_direction(
                    &mut scratch,
                    None,
                    Shift::Restrict(i, side),
                    &mut extraction,
                )? {
                    reductions += 1;
                }
            }
            for j in 0..n {
                if !self.linear_column(j) || self.a.column(j).is_empty() {
                    continue;
                }
                let b = self.bounds[j];
                if b.equality() || !(b.lower.is_finite() || b.upper.is_finite()) {
                    continue;
                }
                let c = self.objective.c[j];
                // Fixing earlier columns and restricting rows leave this
                // column's entries and every multiplier unchanged.
                let activity = match scratch.column_activity[j] {
                    Some(activity) if complete => activity,
                    _ => Activity::compute(self.a.column(j), &scratch.y, None),
                };
                let side = if b.lower.is_finite()
                    && activity
                        .max
                        .value()
                        .is_some_and(|v| self.separated(c - v, 0.0))
                {
                    Side::Lower
                } else if b.upper.is_finite()
                    && activity
                        .min
                        .value()
                        .is_some_and(|v| self.separated(0.0, c - v))
                {
                    Side::Upper
                } else {
                    continue;
                };
                // `z_j > 0` bounded the dual row from above, using the bound
                // of each multiplier that maximizes its term; `z_j < 0` from
                // below. Each term's weight is its coefficient magnitude.
                scratch.weights.clear();
                for (i, a) in self.a.column(j) {
                    let used = used_bound(side, a);
                    if let Some(&sequence) = scratch.proofs[slot(i, used)].last() {
                        *scratch.weights.entry(sequence).or_insert(0.0) += a.abs();
                    }
                }
                if self.apply_direction(
                    &mut scratch,
                    Some((j, side)),
                    Shift::Fix(j, side),
                    &mut extraction,
                )? {
                    reductions += 1;
                }
            }
            Ok(())
        })();
        self.dual_scratch = scratch;
        result.map(|()| reductions)
    }

    /// Accumulate the direction of a proof, verify it exactly on the current
    /// model, and apply the shift. Returns whether a reduction was made.
    fn apply_direction(
        &mut self,
        scratch: &mut DualScratch,
        seed: Option<(usize, Side)>,
        shift: Shift,
        extraction: &mut usize,
    ) -> Result<bool, Certificate> {
        let n = self.bounds.len();
        let m = self.rows.len();
        scratch.lambda.resize(n, 0.0);
        scratch.ad.resize(m, 0.0);
        for &j in &scratch.lambda_touched {
            scratch.lambda[j] = 0.0;
        }
        for &i in &scratch.ad_touched {
            scratch.ad[i] = 0.0;
        }
        scratch.lambda_touched.clear();
        scratch.ad_touched.clear();
        // Walk the proof from the latest record backwards. Each record spends
        // its weight on its column's dual row and on the bounds that row used,
        // which are the latest records older than the record itself.
        while let Some((sequence, weight)) = scratch.weights.pop_last() {
            let record = scratch.records[sequence as usize - 1];
            let j = record.column as usize;
            let column = self.a.column(j);
            if column.len() > *extraction {
                return Ok(false);
            }
            *extraction -= column.len();
            let pivot = self.a.get(record.row as usize, j).abs();
            if pivot == 0.0 {
                return Ok(false);
            }
            let scale = weight / pivot;
            if scratch.lambda[j] == 0.0 {
                scratch.lambda_touched.push(j);
            }
            scratch.lambda[j] += if record.constraint == Side::Lower {
                scale
            } else {
                -scale
            };
            for (k, a) in column {
                if k == record.row as usize {
                    continue;
                }
                let used = used_bound(record.constraint, a);
                let proofs = &scratch.proofs[slot(k, used)];
                let at = proofs.partition_point(|&s| s < sequence);
                if at > 0 {
                    *scratch.weights.entry(proofs[at - 1]).or_insert(0.0) += scale * a.abs();
                }
            }
        }
        // d = -λ, plus the unit step onto the fixed column's bound.
        let mut target = None;
        if let Some((q, side)) = seed {
            if scratch.lambda[q] == 0.0 {
                scratch.lambda_touched.push(q);
            }
            scratch.lambda[q] += if side == Side::Lower { 1.0 } else { -1.0 };
            target = Some(q);
        }
        let mut cost = 0.0;
        for at in 0..scratch.lambda_touched.len() {
            let j = scratch.lambda_touched[at];
            let d = -scratch.lambda[j];
            if d == 0.0 {
                continue;
            }
            if !self.linear_column(j) {
                return Ok(false);
            }
            let b = self.bounds[j];
            let allowed = match (b.lower.is_finite(), b.upper.is_finite()) {
                (true, true) => target == Some(j),
                (true, false) => d > 0.0 || target == Some(j),
                (false, true) => d < 0.0 || target == Some(j),
                (false, false) => true,
            };
            if !allowed || !d.is_finite() {
                return Ok(false);
            }
            cost += self.objective.c[j] * d;
            for (i, a) in self.a.column(j) {
                if scratch.ad[i] == 0.0 {
                    scratch.ad_touched.push(i);
                }
                scratch.ad[i] += a * d;
            }
        }
        if !cost.is_finite() || cost > 0.0 {
            return Ok(false);
        }
        let mut drives = false;
        for &i in &scratch.ad_touched {
            let change = scratch.ad[i];
            if !change.is_finite() {
                return Ok(false);
            }
            let RowDomain::Linear(b) = self.rows[i] else {
                return Ok(false);
            };
            let toward = match shift {
                Shift::Restrict(row, Side::Lower) if row == i => {
                    drives = change < 0.0;
                    true
                }
                Shift::Restrict(row, Side::Upper) if row == i => {
                    drives = change > 0.0;
                    true
                }
                _ => false,
            };
            let allowed = match (b.lower.is_finite(), b.upper.is_finite()) {
                (true, true) => change == 0.0,
                (true, false) => change >= 0.0 || toward,
                (false, true) => change <= 0.0 || toward,
                (false, false) => true,
            };
            if !allowed {
                return Ok(false);
            }
        }
        if let Shift::Fix(q, side) = shift {
            let d = -scratch.lambda[q];
            drives = match side {
                Side::Lower => d < 0.0,
                Side::Upper => d > 0.0,
            };
        }
        if !drives {
            // The direction never reaches the target, so it is a recession
            // direction; with negative cost it proves unboundedness.
            if cost < 0.0 {
                let ray = scratch
                    .lambda_touched
                    .iter()
                    .map(|&j| (j, -scratch.lambda[j]));
                return Err(self.dual_certificate(ray));
            }
            return Ok(false);
        }
        match shift {
            Shift::Fix(q, side) => Ok(self.fix(q, side.value(self.bounds[q]))),
            Shift::Restrict(i, side) => {
                self.restrict_row(i, side);
                Ok(true)
            }
        }
    }

    /// Improve one multiplier bound. `None` reports a contradiction.
    fn tighten_multiplier(&self, y: &mut Bounds, side: Side, value: f64) -> Option<bool> {
        if !value.is_finite() {
            return Some(false);
        }
        let old = side.value(*y);
        let opposite = match side {
            Side::Lower => y.upper,
            Side::Upper => y.lower,
        };
        let (lower, upper) = match side {
            Side::Lower => (value, opposite),
            Side::Upper => (opposite, value),
        };
        if lower > upper {
            return if self.separated(lower, upper) {
                None
            } else {
                Some(false)
            };
        }
        let gain = match side {
            Side::Lower => value - old,
            Side::Upper => old - value,
        };
        if gain <= 0.0 {
            return Some(false);
        }
        if old.is_finite()
            && gain
                <= (self.settings.propagation.minimum_gain_factor
                    * self.settings.numerics.feasibility)
                    .max(self.settings.propagation.minimum_relative_gain * old.abs())
        {
            return Some(false);
        }
        match side {
            Side::Lower => y.lower = value,
            Side::Upper => y.upper = value,
        }
        Some(true)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn multiplier_ranges_follow_native_signs() {
        let lower_only = Bounds {
            lower: 1.0,
            upper: f64::INFINITY,
        };
        let upper_only = Bounds {
            lower: f64::NEG_INFINITY,
            upper: 1.0,
        };
        assert_eq!(
            multiplier_range(lower_only),
            Bounds {
                lower: 0.0,
                upper: f64::INFINITY
            }
        );
        assert_eq!(
            multiplier_range(upper_only),
            Bounds {
                lower: f64::NEG_INFINITY,
                upper: 0.0
            }
        );
        assert_eq!(multiplier_range(Bounds::fixed(1.0)), Bounds::FREE);
        assert_eq!(
            multiplier_range(Bounds {
                lower: 0.0,
                upper: 1.0
            }),
            Bounds::FREE
        );
    }

    #[test]
    fn used_bounds_extremize_each_term() {
        assert_eq!(used_bound(Side::Upper, 2.0), Side::Lower);
        assert_eq!(used_bound(Side::Upper, -2.0), Side::Upper);
        assert_eq!(used_bound(Side::Lower, 2.0), Side::Upper);
        assert_eq!(used_bound(Side::Lower, -2.0), Side::Lower);
    }
}
