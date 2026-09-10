// SPDX-License-Identifier: Apache-2.0
// Modified for this library; copyright and attribution notices are in NOTICE.
//! Mutable working model and shared rule transformations. Every mutation
//! updates sparse views, locks, dirty activities, and rule queues together.

use crate::{
    core::{
        activity::{Activity, Locks},
        objective::Objective,
        queues::Queues,
    },
    matrix::{linked::LinkedMatrix, sparse::Entries},
    postsolve::tape::{Equation, RecoveryTape, Rule, Side},
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
    activities: Vec<Option<Activity>>,
    pub revision: usize,
    pub rules: crate::settings::Rules,
    pub numerics: crate::settings::Numerics,
    pub propagation: crate::settings::PropagationSettings,
    pub equalities: crate::settings::EqualitySettings,
    pub dependencies: crate::settings::DependencySettings,
    pub equality_stats: crate::result::EqualityStats,
    pub substitution_failure: super::objective::SubstitutionFailure,
    pub allow_hessian_growth: bool,
    pub deadline: Option<Instant>,
    pub cones: Vec<crate::problem::Cone>,
    pub cone_rows: Vec<Vec<usize>>,
    pub changed_cones: crate::core::queues::Worklist,
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
            activities: vec![None; m],
            revision: 0,
            rules: crate::settings::Rules::default(),
            numerics: crate::settings::Numerics::default(),
            propagation: crate::settings::PropagationSettings::default(),
            equalities: crate::settings::EqualitySettings::default(),
            dependencies: crate::settings::DependencySettings::default(),
            equality_stats: crate::result::EqualityStats::default(),
            substitution_failure: super::objective::SubstitutionFailure::Numerical,
            allow_hessian_growth: false,
            deadline: None,
            cones: vec![],
            cone_rows: vec![],
            changed_cones: crate::core::queues::Worklist::new(0),
        };
        for i in 0..m {
            for (j, a) in model.a.row(i) {
                model.locks[j].add(Locks::contribution(a, model.rows[i]));
            }
            model.changed_row(i);
        }
        for j in 0..n {
            model.queues.column_changed(j, model.a.column(j).len());
            if model.bounds[j].equality() {
                model.queues.fixed_columns.push(j);
            }
        }
        model
    }

    pub fn add_variable(&mut self, bounds: Bounds) -> usize {
        let j = self.a.add_column();
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
                self.queues.column_changed(j, self.a.column(j).len());
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

    pub fn equation(&self, row: usize) -> Option<Arc<Equation>> {
        let RowDomain::Linear(bounds) = self.rows[row] else {
            return None;
        };
        Some(Arc::new(Equation {
            row,
            entries: self.a.row(row).to_vec(),
            bounds,
        }))
    }

    fn changed_row(&mut self, row: usize) {
        if let RowDomain::Cone { block, .. } = self.rows[row] {
            self.changed_cones.push(block);
        }
        self.activities[row] = None;
        if self.rows[row] != RowDomain::Deleted {
            self.queues.row_changed(
                row,
                self.a.row(row).len(),
                matches!(self.rows[row], RowDomain::Linear(b) if b.equality()),
            );
        }
    }

    /// Structural edits are infrequent relative to activity queries. Cache
    /// activities until an incident bound or coefficient changes; recomputing
    /// dirty rows also avoids cumulative subtract/add cancellation error.
    pub fn activity(&mut self, row: usize) -> Activity {
        *self.activities[row]
            .get_or_insert_with(|| Activity::compute(self.a.row(row), &self.bounds, None))
    }

    pub fn residual_activity(&mut self, row: usize, column: usize) -> Activity {
        // This direct path is the cancellation fallback for the propagation
        // rule, which first tries subtracting from cached extremes.
        Activity::compute(self.a.row(row), &self.bounds, Some(column))
    }

    pub(super) fn replace_row(&mut self, row: usize, entries: &[(usize, f64)], domain: RowDomain) {
        let old = self.a.replace_row(row, entries);
        for &(j, a) in &old {
            self.locks[j].remove(Locks::contribution(a, self.rows[row]));
        }
        self.rows[row] = domain;
        for &(j, a) in entries {
            self.locks[j].add(Locks::contribution(a, domain));
        }
        for &(j, _) in old.iter().chain(entries) {
            if self.alive[j] {
                self.queues.column_changed(j, self.a.column(j).len());
            }
        }
        self.changed_row(row);
        self.revision += 1;
    }

    pub fn delete_row(&mut self, row: usize) {
        assert!(matches!(self.rows[row], RowDomain::Linear(_)));
        self.replace_row(row, &[], RowDomain::Deleted);
        self.postsolve.rules.push(Rule::DeletedRow(row));
    }

    fn set_bounds(&mut self, column: usize, bounds: Bounds) {
        let old = std::mem::replace(&mut self.bounds[column], bounds);
        if bounds.equality() {
            self.queues.fixed_columns.push(column);
        }
        for (i, a) in self.a.column(column) {
            if self.activities[i]
                .as_mut()
                .is_some_and(|activity| !activity.replace_bound(a, old, bounds))
            {
                self.activities[i] = None;
            }
            self.queues.changed_activities.push(i);
            if let RowDomain::Cone { block, .. } = self.rows[i] {
                self.changed_cones.push(block);
            }
            self.queues.singleton_activity_rows.push(i);
            if !self.equalities.require_free_variable
                && matches!(self.rows[i], RowDomain::Linear(b) if b.equality())
                && self.a.row(i).len() >= 3
            {
                self.queues.short_equalities.push(i);
            }
            if matches!(self.rows[i],RowDomain::Linear(b) if b.equality())
                && self.a.row(i).len() == 2
            {
                self.queues.doubleton_rows.push(i);
            }
        }
        self.queues
            .column_changed(column, self.a.column(column).len());
        self.revision += 1;
    }

    pub fn relax_bound(&mut self, column: usize, side: Side) {
        let mut bounds = self.bounds[column];
        match side {
            Side::Lower => bounds.lower = f64::NEG_INFINITY,
            Side::Upper => bounds.upper = f64::INFINITY,
        }
        self.set_bounds(column, bounds);
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
            self.queues.column_changed(j, self.a.column(j).len());
        }
        self.alive[column] = false;
        let entries = self.a.remove_column(column);
        for ((i, domain), &(row, a)) in updates.into_iter().zip(&entries) {
            debug_assert_eq!(i, row);
            self.locks[column].remove(Locks::contribution(a, self.rows[i]));
            self.rows[i] = domain;
            self.changed_row(i);
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
        let Some(mut equation) = self.equation(row) else {
            return false;
        };
        if !side.is_finite() || !self.alive[column] {
            return false;
        }
        Arc::make_mut(&mut equation).bounds = Bounds::fixed(side);
        self.substitute_equation(column, equation, effective_bounds, max_fill)
    }

    fn substitute_equation(
        &mut self,
        column: usize,
        equation: Arc<Equation>,
        effective_bounds: Bounds,
        max_fill: usize,
    ) -> bool {
        use super::objective::SubstitutionFailure;
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
            self.allow_hessian_growth,
            self.deadline,
        ) {
            Ok(gradient) => gradient,
            Err(reason) => {
                self.substitution_failure = reason;
                return false;
            }
        };
        for &(j, _) in gradient.terms.iter().chain(&slopes) {
            self.queues.column_changed(j, self.a.column(j).len());
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

    /// Symbolic-infinity elimination, restricted to a flat quadratic
    /// direction. The rule supplies the unlocked direction proof.
    pub fn remove_unlocked(&mut self, column: usize) {
        assert!(self.objective.p.column(column).is_empty() && self.objective.c[column] == 0.0);
        let rows: Vec<_> = self
            .a
            .column(column)
            .iter()
            .map(|(i, _)| {
                self.equation(i)
                    .expect("a cone row locks both directions")
                    .as_ref()
                    .clone()
            })
            .collect();
        self.alive[column] = false;
        for equation in &rows {
            self.replace_row(equation.row, &[], RowDomain::Deleted);
        }
        self.postsolve.rules.push(Rule::Unlocked {
            column,
            bounds: self.bounds[column],
            rows,
        });
        self.revision += 1;
    }
}
