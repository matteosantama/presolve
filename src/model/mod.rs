// SPDX-License-Identifier: Apache-2.0
// Modified for this library; copyright and attribution notices are in NOTICE.
//! Mutable working model and shared rule transformations. Every mutation
//! updates sparse views, locks, dirty activities, and rule queues together.
//! Rule families in `crate::rules` drive the mutations; the recovery tape
//! records them for `crate::postsolve` to reverse.
//! Native dual signs satisfy:
//! `P x + c = A^T y + z`, with positive multipliers on lower bounds.

pub(crate) mod activity;
pub(crate) mod objective;
pub(crate) mod queues;
pub(crate) mod tape;

use crate::{
    matrix::{linked::LinkedMatrix, sparse::Entries},
    model::{
        activity::{Activity, Extreme, Locks},
        objective::Objective,
        queues::Queues,
        tape::{Certificate, Equation, Point, Recovery, RecoveryTape, Rule, Side},
    },
    problem::Bounds,
};
use std::{sync::Arc, time::Instant};

/// Working row domains; conic coordinates remain distinct from ranged linear rows.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) enum RowDomain {
    Linear(Bounds),
    Cone { rhs: f64, block: usize },
    Deleted,
}

pub(crate) struct Model {
    pub a: LinkedMatrix,
    pub objective: Objective,
    pub rows: Vec<RowDomain>,
    pub bounds: Vec<Bounds>,
    pub alive: Vec<bool>,
    pub locks: Vec<Locks>,
    pub queues: Queues,
    pub postsolve: RecoveryTape,
    /// Cached row activities; a stale entry carries the `STALE` marker
    /// instead of an `Option` discriminant, keeping the array at 32 bytes
    /// per row for the scattered accesses made on every bound change.
    activities: Vec<Activity>,
    /// Previous coefficients of the row being replaced, kept between calls
    /// so row edits do not allocate a fresh copy each time.
    row_scratch: Entries,
    /// One `RowKind` byte per row, refreshed by `changed_row`. A bound change
    /// visits every incident row's queue eligibility; the byte replaces reads
    /// of the wider domain and list-header arrays for that decision.
    row_kinds: Vec<u8>,
    /// Snapshot of each linear row for the tape, shared until the row's
    /// coefficients or sides change. Repeated propagation rounds over an
    /// unchanged row then reuse one copy instead of taking one per round.
    equations: Vec<Option<Arc<Equation>>>,
    pub revision: usize,
    /// Run configuration, applied once by `configure`.
    pub settings: crate::settings::Settings,
    pub dual_scratch: crate::rules::dual_propagation::DualScratch,
    pub dominated_scratch: crate::rules::dominated_columns::DominatedScratch,
    pub equality_stats: crate::result::EqualityStats,
    pub substitution_failure: self::objective::SubstitutionFailure,
    /// Hessian revision plus one at which a coupled column's elimination was
    /// last rejected, so cleanup does not retry it on every visit.
    pub elimination_rejected: Vec<usize>,
    pub deadline: Option<Instant>,
    pub cones: Vec<crate::problem::Cone>,
    pub cone_rows: Vec<Vec<usize>>,
    pub changed_cones: crate::model::queues::Worklist,
}

/// A row never has this many infinite terms, so the count marks a cached
/// activity that must be recomputed before use.
const STALE: usize = usize::MAX;

/// Queue eligibility of a row on a bound change, derived from its domain
/// and length. Every domain or length change ends in `changed_row`, which
/// keeps the classification current.
mod row_kind {
    pub const OTHER: u8 = 0;
    pub const CONE: u8 = 1;
    pub const DOUBLETON_EQUALITY: u8 = 2;
    pub const LONGER_EQUALITY: u8 = 3;
}

fn stale() -> Activity {
    Activity {
        min: Extreme {
            sum: 0.0,
            infinite: STALE,
        },
        max: Extreme::default(),
    }
}

/// Shift a cached activity for a bound change of one of its terms, marking
/// it stale when the incremental update would lose precision.
fn shift_cached(activity: &mut Activity, a: f64, old: Bounds, new: Bounds) {
    if activity.min.infinite != STALE && !activity.replace_bound(a, old, new) {
        *activity = stale();
    }
}

/// Translate a right-hand side under x_k = offset + ... . Infinite sides
/// remain infinite; overflow of a formerly finite side rejects the rule.
pub(crate) fn shifted(domain: RowDomain, shift: f64) -> Option<RowDomain> {
    if !shift.is_finite() {
        return None;
    }
    let side = |value: f64| {
        if value.is_finite() {
            (value - shift).is_finite().then_some(value - shift)
        } else {
            Some(value)
        }
    };
    Some(match domain {
        RowDomain::Linear(b) => RowDomain::Linear(Bounds {
            lower: side(b.lower)?,
            upper: side(b.upper)?,
        }),
        RowDomain::Cone { rhs, block } => RowDomain::Cone {
            rhs: side(rhs)?,
            block,
        },
        RowDomain::Deleted => return None,
    })
}

impl Model {
    pub fn set_cones(&mut self, cones: Vec<crate::problem::Cone>) {
        if cones.is_empty() {
            return;
        }
        self.cone_rows = vec![Vec::new(); cones.len()];
        for (i, row) in self.rows.iter().enumerate() {
            if let RowDomain::Cone { block, .. } = *row {
                self.cone_rows[block].push(i);
                self.changed_cones.push(block);
            }
        }
        self.cones = cones;
    }

    pub fn from_parts(
        a: LinkedMatrix,
        objective: Objective,
        rows: Vec<RowDomain>,
        bounds: Vec<Bounds>,
    ) -> Self {
        let n = bounds.len();
        let m = rows.len();
        let mut model = Self {
            a,
            objective,
            rows,
            bounds,
            alive: vec![true; n],
            locks: vec![Locks::default(); n],
            queues: Queues::new(m, n),
            postsolve: RecoveryTape::default(),
            activities: vec![stale(); m],
            row_scratch: Vec::new(),
            row_kinds: vec![row_kind::OTHER; m],
            equations: vec![None; m],
            revision: 0,
            settings: crate::settings::Settings::default(),
            dual_scratch: Default::default(),
            dominated_scratch: Default::default(),
            equality_stats: crate::result::EqualityStats::default(),
            substitution_failure: self::objective::SubstitutionFailure::Numerical,
            elimination_rejected: vec![0; n],
            deadline: None,
            cones: vec![],
            cone_rows: vec![],
            changed_cones: crate::model::queues::Worklist::new(0),
        };
        for i in 0..m {
            for (j, a) in model.a.row(i) {
                model.locks[j].add(Locks::contribution(a, model.rows[i]));
            }
            model.changed_row(i);
        }
        for j in 0..n {
            model.column_changed(j);
            if model.bounds[j].equality() {
                model.queues.fixed_columns.push(j);
            }
        }
        model
    }

    /// A linear column: alive, absent from the Hessian, and absent from every
    /// conic row, so its reduced cost is `c_j - Σ a_ij y_i` over linear rows.
    pub fn linear_column(&self, j: usize) -> bool {
        self.alive[j]
            && self.objective.p.column(j).is_empty()
            && (self.cones.is_empty()
                || self
                    .a
                    .column(j)
                    .iter()
                    .all(|(i, _)| matches!(self.rows[i], RowDomain::Linear(_))))
    }

    /// Apply the caller's settings once; rules read them from `self.settings`.
    pub fn configure(&mut self, settings: &crate::settings::Settings, deadline: Option<Instant>) {
        self.settings = settings.clone();
        self.settings.propagation = settings.propagation.sanitized();
        self.deadline = deadline;
    }

    pub fn add_variable(&mut self, bounds: Bounds) -> usize {
        let j = self.a.add_column();
        self.elimination_rejected.push(0);
        self.objective.p.add_variable();
        self.objective.c.push(0.0);
        self.bounds.push(bounds);
        self.alive.push(true);
        self.locks.push(Locks::default());
        self.revision += 1;
        j
    }

    pub fn replace_rows_batch(&mut self, updates: Vec<(usize, Entries, RowDomain)>) {
        let mut changed = Vec::with_capacity(updates.len());
        let mut matrix_updates = Vec::with_capacity(updates.len());
        for (i, entries, domain) in updates {
            for (j, a) in self.a.row(i) {
                self.locks[j].remove(Locks::contribution(a, self.rows[i]));
            }
            for &(j, a) in &entries {
                self.locks[j].add(Locks::contribution(a, domain));
            }
            self.rows[i] = domain;
            changed.push(i);
            matrix_updates.push((i, entries));
        }
        for j in self.a.replace_rows(matrix_updates) {
            if self.alive[j] {
                self.column_changed(j);
            }
        }
        // Domain changes also affect locks of coefficients that stayed equal.
        for i in changed {
            for (j, _) in self.a.row(i) {
                self.queues.unlocked_columns.push(j);
            }
            self.changed_row(i);
            self.revision += 1;
        }
    }

    pub fn equation(&mut self, row: usize) -> Option<Arc<Equation>> {
        let RowDomain::Linear(bounds) = self.rows[row] else {
            return None;
        };
        let equation = self.equations[row].get_or_insert_with(|| {
            Arc::new(Equation {
                row,
                entries: self.a.row(row).to_vec(),
                bounds,
            })
        });
        Some(Arc::clone(equation))
    }

    #[inline]
    fn changed_row(&mut self, row: usize) {
        let domain = self.rows[row];
        if let RowDomain::Cone { block, .. } = domain {
            self.changed_cones.push(block);
        }
        self.activities[row] = stale();
        self.equations[row] = None;
        let length = self.a.row(row).len();
        let equality = matches!(domain, RowDomain::Linear(b) if b.equality());
        self.row_kinds[row] = match domain {
            RowDomain::Cone { .. } => row_kind::CONE,
            RowDomain::Linear(_) if equality && length == 2 => row_kind::DOUBLETON_EQUALITY,
            RowDomain::Linear(_) if equality && length >= 3 => row_kind::LONGER_EQUALITY,
            _ => row_kind::OTHER,
        };
        if domain != RowDomain::Deleted {
            self.queues.row_changed(row, length, equality);
        }
    }

    /// Structural edits are infrequent relative to activity queries. Cache
    /// activities until an incident bound or coefficient changes; recomputing
    /// dirty rows also avoids cumulative subtract/add cancellation error.
    #[inline]
    pub fn activity(&mut self, row: usize) -> Activity {
        if self.activities[row].min.infinite == STALE {
            self.activities[row] = Activity::compute(self.a.row(row), &self.bounds, None);
        }
        self.activities[row]
    }

    pub fn residual_activity(&self, row: usize, column: usize) -> Activity {
        // This direct path is the cancellation fallback for the propagation
        // rule, which first tries subtracting from cached extremes.
        Activity::compute(self.a.row(row), &self.bounds, Some(column))
    }

    pub(super) fn replace_row(&mut self, row: usize, entries: &[(usize, f64)], domain: RowDomain) {
        let mut old = std::mem::take(&mut self.row_scratch);
        self.a.replace_row_into(row, entries, &mut old);
        let previous = std::mem::replace(&mut self.rows[row], domain);
        // Both supports are sorted, so merge them. A column whose lock
        // contribution is unchanged would only be decremented and then
        // incremented by the same amount, so it is left alone. Removed and
        // retained columns are queued first, in the old support's order; a
        // retained column is already queued by then, so the second pass
        // only visits columns that are new to the row.
        let mut next = 0;
        for &(j, a) in &old {
            let removed = Locks::contribution(a, previous);
            while next < entries.len() && entries[next].0 < j {
                let (k, b) = entries[next];
                self.locks[k].add(Locks::contribution(b, domain));
                next += 1;
            }
            if next < entries.len() && entries[next].0 == j {
                let added = Locks::contribution(entries[next].1, domain);
                if added != removed {
                    self.locks[j].remove(removed);
                    self.locks[j].add(added);
                }
                next += 1;
            } else {
                self.locks[j].remove(removed);
            }
        }
        for &(k, b) in &entries[next..] {
            self.locks[k].add(Locks::contribution(b, domain));
        }
        for &(j, _) in &old {
            if self.alive[j] {
                self.column_changed(j);
            }
        }
        let mut next = 0;
        for &(j, _) in entries {
            while next < old.len() && old[next].0 < j {
                next += 1;
            }
            if next < old.len() && old[next].0 == j {
                next += 1;
                continue;
            }
            if self.alive[j] {
                self.column_changed(j);
            }
        }
        self.row_scratch = old;
        self.changed_row(row);
        self.revision += 1;
    }

    /// Remove a row's coefficients and retire its domain, without a record;
    /// the caller pushes the rule that explains the row's multiplier.
    pub(super) fn clear_row(&mut self, row: usize) {
        debug_assert!(self.rows[row] != RowDomain::Deleted);
        self.replace_row(row, &[], RowDomain::Deleted);
    }

    /// Delete a row whose multiplier is simply zeroed on recovery.
    pub fn delete_row(&mut self, row: usize) {
        self.clear_row(row);
        self.postsolve.rules.push(Rule::DeletedRow(row));
    }

    fn column_changed(&mut self, j: usize) {
        self.queues.column_changed(j, self.a.column(j).len());
    }

    /// A zero point in working coordinates, for certificates and probes.
    pub(crate) fn point(&self) -> Point {
        Point::zeros(self.bounds.len(), self.rows.len())
    }

    /// Farkas certificate with multipliers `y` on rows and `z` on columns.
    pub(super) fn primal_certificate(
        &self,
        y: impl IntoIterator<Item = (usize, f64)>,
        z: impl IntoIterator<Item = (usize, f64)>,
    ) -> Certificate {
        let mut point = self.point();
        for (i, v) in y {
            point.y[i] = v;
        }
        for (j, v) in z {
            point.z[j] = v;
        }
        Certificate {
            mode: Recovery::PrimalInfeasibility,
            point,
        }
    }

    /// Recession direction with the given nonzero components.
    pub(super) fn dual_certificate(
        &self,
        x: impl IntoIterator<Item = (usize, f64)>,
    ) -> Certificate {
        let mut point = self.point();
        for (j, v) in x {
            point.x[j] = v;
        }
        Certificate {
            mode: Recovery::DualInfeasibility,
            point,
        }
    }

    fn set_bounds(&mut self, column: usize, bounds: Bounds) {
        let old = std::mem::replace(&mut self.bounds[column], bounds);
        if bounds.equality() {
            self.queues.fixed_columns.push(column);
        }
        let require_free = self.settings.equalities.require_free_variable;
        for (i, a) in self.a.column(column) {
            shift_cached(&mut self.activities[i], a, old, bounds);
            self.queues.changed_activities.push(i);
            // The kind byte decides the remaining queues, so the wider domain
            // is only read for a cone row's block.
            match self.row_kinds[i] {
                row_kind::CONE => {
                    if let RowDomain::Cone { block, .. } = self.rows[i] {
                        self.changed_cones.push(block);
                    }
                }
                row_kind::LONGER_EQUALITY if !require_free => {
                    self.queues.short_equalities.push(i);
                }
                row_kind::DOUBLETON_EQUALITY => self.queues.doubleton_rows.push(i),
                _ => {}
            }
        }
        self.column_changed(column);
        self.revision += 1;
    }

    /// Only the final redundant-bound pass relaxes bounds, after every rule
    /// that drains a work queue has finished. Cached activities and the
    /// revision must stay consistent for the remaining columns of that pass;
    /// queue entries would never be consumed, so `set_bounds` is not used.
    pub fn relax_bound(&mut self, column: usize, side: Side) {
        let mut bounds = self.bounds[column];
        match side {
            Side::Lower => bounds.lower = f64::NEG_INFINITY,
            Side::Upper => bounds.upper = f64::INFINITY,
        }
        let old = std::mem::replace(&mut self.bounds[column], bounds);
        for (i, a) in self.a.column(column) {
            shift_cached(&mut self.activities[i], a, old, bounds);
        }
        self.revision += 1;
    }

    /// The rule checks feasibility and the implication proof before
    /// calling this primitive. Reversal attributes the new bound to its row.
    pub fn tighten_bound(
        &mut self,
        column: usize,
        side: Side,
        value: f64,
        equation: Arc<Equation>,
    ) -> bool {
        if !value.is_finite() {
            return false;
        }
        let mut bounds = self.bounds[column];
        let old = side.value(bounds);
        match side {
            Side::Lower if value > old => bounds.lower = value,
            Side::Upper if value < old => bounds.upper = value,
            _ => return false,
        }
        self.postsolve.rules.push(Rule::TightenedBound {
            column,
            equation,
            side,
            old,
        });
        self.set_bounds(column, bounds);
        true
    }

    /// Fixing changes coupled linear coefficients as well as row sides.
    /// No state changes if any transformed coefficient would be nonfinite.
    pub fn fix(&mut self, column: usize, value: f64) -> bool {
        assert!(self.alive[column]);
        let Some(updates) = self
            .a
            .column(column)
            .iter()
            .map(|(i, a)| shifted(self.rows[i], a * value).map(|domain| (i, domain)))
            .collect::<Option<Vec<_>>>()
        else {
            return false;
        };
        let Some(gradient) = self
            .objective
            .substitute(column, value, &[], 0, false, None)
        else {
            return false;
        };
        for &(j, _) in &gradient.terms {
            self.column_changed(j);
        }
        self.alive[column] = false;
        let bounds = self.bounds[column];
        let entries = self.a.remove_column(column);
        for ((i, domain), &(row, a)) in updates.into_iter().zip(&entries) {
            debug_assert_eq!(i, row);
            self.locks[column].remove(Locks::contribution(a, self.rows[i]));
            self.rows[i] = domain;
            // Keep the cached activity warm: the removed term is the column's
            // old contribution, so subtracting it is the same incremental
            // update a bound change makes. Rules that fix many columns of
            // shared long rows would otherwise recompute those rows each time.
            let mut kept = self.activities[i];
            shift_cached(&mut kept, a, bounds, Bounds::fixed(0.0));
            self.changed_row(i);
            self.activities[i] = kept;
        }
        self.postsolve.rules.push(Rule::Fixed {
            column,
            value,
            gradient,
            entries,
        });
        self.revision += 1;
        true
    }

    /// Equality substitution shared by singleton columns and doubleton rows.
    /// `effective_bounds` may omit sides proved implied by the rule. Any
    /// remaining sides become a ranged row; a doubleton's resulting singleton
    /// can then be converted to a variable bound.
    pub fn substitute(
        &mut self,
        row: usize,
        column: usize,
        effective_bounds: Bounds,
        max_fill: usize,
    ) -> bool {
        let Some(equation) = self.equation(row) else {
            return false;
        };
        if !equation.bounds.equality() || !self.alive[column] {
            return false;
        }
        self.substitute_equation(column, equation, effective_bounds, max_fill)
    }

    /// The singleton rule has proved this inequality side can be active
    /// at an optimum. Keep the change transactional with the substitution.
    pub fn substitute_on_side(
        &mut self,
        row: usize,
        column: usize,
        side: f64,
        effective_bounds: Bounds,
        max_fill: usize,
    ) -> bool {
        if !matches!(self.rows[row], RowDomain::Linear(_))
            || !side.is_finite()
            || !self.alive[column]
        {
            return false;
        }
        // The tape needs the row at the chosen side, not the shared snapshot.
        let equation = Arc::new(Equation {
            row,
            entries: self.a.row(row).to_vec(),
            bounds: Bounds::fixed(side),
        });
        self.substitute_equation(column, equation, effective_bounds, max_fill)
    }

    fn substitute_equation(
        &mut self,
        column: usize,
        equation: Arc<Equation>,
        effective_bounds: Bounds,
        max_fill: usize,
    ) -> bool {
        use self::objective::SubstitutionFailure;
        self.substitution_failure = SubstitutionFailure::Numerical;
        let row = equation.row;
        let pivot = self.a.get(row, column);
        if pivot == 0.0 {
            return false;
        }
        let offset = equation.bounds.lower / pivot;
        let slopes: Entries = equation
            .entries
            .iter()
            .filter(|&&(j, _)| j != column)
            .map(|&(j, a)| (j, -a / pivot))
            .collect();
        if !offset.is_finite() || slopes.iter().any(|&(_, v)| !v.is_finite()) {
            return false;
        }
        let remaining: Entries = equation
            .entries
            .iter()
            .copied()
            .filter(|&(j, _)| j != column)
            .collect();
        let retained = effective_bounds != Bounds::FREE;
        let domain = if retained {
            let (l, u) = if pivot > 0.0 {
                (effective_bounds.upper, effective_bounds.lower)
            } else {
                (effective_bounds.lower, effective_bounds.upper)
            };
            let l_new = equation.bounds.lower - pivot * l;
            let u_new = equation.bounds.upper - pivot * u;
            if (l.is_finite() && !l_new.is_finite()) || (u.is_finite() && !u_new.is_finite()) {
                return false;
            }
            RowDomain::Linear(Bounds {
                lower: l_new,
                upper: u_new,
            })
        } else {
            RowDomain::Deleted
        };
        let other_rows: Entries = self
            .a
            .column(column)
            .iter()
            .filter(|&(i, _)| i != row)
            .collect();
        let mut fill = 0;
        let mut updates = Vec::with_capacity(other_rows.len());
        for &(i, a) in &other_rows {
            if self
                .deadline
                .is_some_and(|deadline| Instant::now() >= deadline)
            {
                self.substitution_failure = SubstitutionFailure::Deadline;
                return false;
            }
            let Some(domain) = shifted(self.rows[i], a * offset) else {
                return false;
            };
            // Merge the existing row and sorted substitution slopes directly.
            // Each old coefficient is read once, with no matrix searches.
            let mut original = self
                .a
                .row(i)
                .iter()
                .filter(|&(j, _)| j != column)
                .peekable();
            let mut changes = slopes.iter().copied().peekable();
            let mut entries = Vec::with_capacity(self.a.row(i).len() + slopes.len());
            while original.peek().is_some() || changes.peek().is_some() {
                let old_column = original.peek().map_or(usize::MAX, |e| e.0);
                let changed_column = changes.peek().map_or(usize::MAX, |e| e.0);
                if old_column < changed_column {
                    entries.push(original.next().unwrap());
                    continue;
                }
                let (j, v) = changes.next().unwrap();
                let old = if old_column == j {
                    original.next().unwrap().1
                } else {
                    0.0
                };
                let change = a * v;
                let value = old + change;
                // Do not turn cancellation roundoff into a tiny pivot in a
                // later equality. Reject the whole substitution transaction;
                // exact cancellations remain valid structural zeros.
                if !value.is_finite()
                    || (value != 0.0 && value.abs() < 1e-10 * old.abs().max(change.abs()))
                {
                    return false;
                }
                if value != 0.0 && old == 0.0 {
                    fill += 1;
                    if fill > max_fill {
                        self.substitution_failure = SubstitutionFailure::ConstraintFill;
                        return false;
                    }
                }
                if value != 0.0 {
                    entries.push((j, value));
                }
            }
            updates.push((i, entries, domain));
        }
        let gradient = match self.objective.try_substitute(
            column,
            offset,
            &slopes,
            max_fill - fill,
            self.settings.allow_hessian_growth,
            self.deadline,
        ) {
            Ok(gradient) => gradient,
            Err(reason) => {
                self.substitution_failure = reason;
                return false;
            }
        };
        for &(j, _) in gradient.terms.iter().chain(&slopes) {
            self.column_changed(j);
        }
        self.alive[column] = false;
        for (i, entries, domain) in updates {
            self.replace_row(i, &entries, domain);
        }
        self.replace_row(row, if retained { &remaining } else { &[] }, domain);
        self.postsolve.rules.push(Rule::Substituted {
            column,
            equation,
            gradient,
            other_rows,
            retained,
        });
        self.revision += 1;
        true
    }

    /// Row sides change locks and eligibility, but not matrix coefficients.
    fn replace_row_bounds(&mut self, row: usize, bounds: Bounds) {
        let old = self.rows[row];
        let domain = RowDomain::Linear(bounds);
        for (j, a) in self.a.row(row) {
            self.locks[j].remove(Locks::contribution(a, old));
            self.locks[j].add(Locks::contribution(a, domain));
            if self.alive[j] {
                self.queues.column_changed(j, self.a.column(j).len());
            }
        }
        self.rows[row] = domain;
        // Preserve the original queue order and fresh-activity arithmetic.
        self.changed_row(row);
        self.revision += 1;
    }

    /// Change a parallel row's side after the rule verifies proportionality.
    /// Deletion of the source must follow this record, so reversal initializes
    /// its multiplier to zero before transferring the tightened side's dual.
    pub fn tighten_row(&mut self, row: usize, source: usize, ratio: f64, side: Side, value: f64) {
        let RowDomain::Linear(mut bounds) = self.rows[row] else {
            unreachable!()
        };
        match side {
            Side::Lower => bounds.lower = value,
            Side::Upper => bounds.upper = value,
        }
        self.replace_row_bounds(row, bounds);
        self.postsolve.rules.push(Rule::TightenedRow {
            row,
            source,
            ratio,
            side,
        });
    }

    /// Removing a redundant side needs no dual transformation: the reduced
    /// row's remaining multiplier already has a sign valid for the original.
    pub fn relax_row(&mut self, row: usize, bounds: Bounds) {
        self.replace_row_bounds(row, bounds);
    }

    /// Restrict a row to one of its finite sides. The caller has a direction
    /// that moves any feasible point onto that side without leaving the
    /// feasible set or increasing the objective, so the reduced equality's
    /// multiplier already has the side's sign and no record is needed.
    pub fn restrict_row(&mut self, row: usize, side: Side) {
        let RowDomain::Linear(bounds) = self.rows[row] else {
            unreachable!()
        };
        self.replace_row_bounds(row, Bounds::fixed(side.value(bounds)));
    }

    pub fn aggregate(&mut self, keep: usize, removed: usize, ratio: f64, tolerance: f64) -> bool {
        if !ratio.is_finite()
            || ratio == 0.0
            || !self.objective.parallel(keep, removed, ratio, tolerance)
        {
            return false;
        }
        let b = self.bounds[keep];
        let c = self.bounds[removed];
        let (l, u) = if ratio > 0.0 {
            (c.lower, c.upper)
        } else {
            (c.upper, c.lower)
        };
        let lower = b.lower + ratio * l;
        let upper = b.upper + ratio * u;
        if lower.is_nan()
            || upper.is_nan()
            || (b.lower.is_finite() && l.is_finite() && !lower.is_finite())
            || (b.upper.is_finite() && u.is_finite() && !upper.is_finite())
        {
            return false;
        }
        self.postsolve.rules.push(Rule::ParallelColumns {
            keep,
            removed,
            ratio,
            keep_bounds: b,
            removed_bounds: c,
        });
        self.set_bounds(keep, Bounds { lower, upper });
        self.alive[removed] = false;
        for (i, a) in self.a.remove_column(removed) {
            self.locks[removed].remove(Locks::contribution(a, self.rows[i]));
            self.changed_row(i);
        }
        for &(j, _) in self.objective.p.column(removed) {
            self.queues.column_changed(j, self.a.column(j).len());
        }
        self.objective.aggregate(removed);
        self.revision += 1;
        true
    }

    /// Minimize a free column that no row contains out of the objective. The
    /// stationarity condition `p_jj x_j + Σ p_jk x_k + c_j = 0` is an affine
    /// substitution, so the Hessian update is the same Schur complement the
    /// equality substitutions apply, under the same fill limit. Net Hessian
    /// growth is never allowed here: a chain of such eliminations is a dense
    /// factorization, which is the solver's job.
    pub fn eliminate_coupled(&mut self, column: usize, max_fill: usize) -> bool {
        debug_assert!(self.a.column(column).is_empty() && self.bounds[column] == Bounds::FREE);
        let diagonal = self.objective.p.get(column, column);
        if !diagonal.is_finite() || diagonal <= 0.0 {
            return false;
        }
        let offset = -self.objective.c[column] / diagonal;
        let slopes: Entries = self
            .objective
            .p
            .row(column)
            .iter()
            .filter(|&&(k, _)| k != column)
            .map(|&(k, p)| (k, -p / diagonal))
            .collect();
        if !offset.is_finite() || slopes.iter().any(|&(_, s)| !s.is_finite()) {
            return false;
        }
        if self
            .objective
            .try_substitute(column, offset, &slopes, max_fill, false, self.deadline)
            .is_err()
        {
            return false;
        }
        for &(k, _) in &slopes {
            self.column_changed(k);
        }
        self.alive[column] = false;
        self.postsolve.rules.push(Rule::Eliminated {
            column,
            offset,
            slopes,
        });
        self.revision += 1;
        true
    }

    /// Symbolic-infinity elimination, restricted to a flat quadratic
    /// direction. The rule supplies the unlocked direction proof.
    pub fn remove_unlocked(&mut self, column: usize) {
        assert!(self.objective.p.column(column).is_empty() && self.objective.c[column] == 0.0);
        let rows: Vec<_> = self
            .a
            .column(column)
            .iter()
            .map(|(i, _)| i)
            .collect::<Vec<_>>()
            .into_iter()
            .map(|i| {
                self.equation(i)
                    .expect("a cone row locks both directions")
                    .as_ref()
                    .clone()
            })
            .collect();
        self.alive[column] = false;
        for equation in &rows {
            self.clear_row(equation.row);
        }
        self.postsolve.rules.push(Rule::Unlocked {
            column,
            bounds: self.bounds[column],
            rows,
        });
        self.revision += 1;
    }
}
