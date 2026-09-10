// SPDX-License-Identifier: Apache-2.0
// Modified for this library; copyright and attribution notices are in NOTICE.
//! Reverse applied rules for solutions and certificates. A removed column's
//! objective derivative is evaluated after recovering its primal value,
//! supporting coupled quadratic objectives and sequences of substitutions.

use crate::{core::objective::Gradient, matrix::sparse::Entries, problem::Bounds};
use std::sync::Arc;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Recovery {
    Solution,
    /// Native-sign Farkas multipliers: A^T y + z = 0.
    PrimalInfeasibility,
    /// A feasible recession direction with negative objective slope.
    DualInfeasibility,
}

impl Recovery {
    fn primal(self) -> bool {
        self != Self::PrimalInfeasibility
    }
    fn dual(self) -> bool {
        self != Self::DualInfeasibility
    }
    fn offset(self, value: f64) -> f64 {
        if self == Self::Solution { value } else { 0.0 }
    }
    fn bounds(self, bounds: Bounds) -> Bounds {
        if self == Self::Solution {
            bounds
        } else {
            bounds.recession()
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Side {
    Lower,
    Upper,
}

impl Side {
    pub fn used(self, multiplier: f64) -> bool {
        match self {
            Self::Lower => multiplier > 0.0,
            Self::Upper => multiplier < 0.0,
        }
    }
    pub fn value(self, bounds: Bounds) -> f64 {
        match self {
            Self::Lower => bounds.lower,
            Self::Upper => bounds.upper,
        }
    }
}

#[derive(Clone, Debug)]
pub(crate) struct Equation {
    pub row: usize,
    pub entries: Entries,
    pub bounds: Bounds,
}

impl Equation {
    fn coefficient(&self, column: usize) -> f64 {
        let index = self
            .entries
            .binary_search_by_key(&column, |&(j, _)| j)
            .expect("saved row contains the eliminated or tightened column");
        self.entries[index].1
    }
    fn activity(&self, x: &[f64], exclude: Option<usize>) -> f64 {
        self.entries
            .iter()
            .filter(|&&(j, _)| Some(j) != exclude)
            .map(|&(j, a)| a * x[j])
            .sum()
    }
}

/// An applied rule and the data needed to reverse its changes.
#[derive(Clone, Debug)]
pub(crate) enum Rule {
    /// Evaluate at this point in reverse traversal, before earlier variable
    /// aggregations/substitutions change the meaning of the stored coefficients.
    ConeSlack {
        row: usize,
        rhs: f64,
        entries: Entries,
    },
    SocToLinear {
        head: usize,
        tail: usize,
    },
    SocFace {
        head: usize,
        tail: Vec<usize>,
    },
    PsdZeroFace {
        rows: Vec<usize>,
        order: usize,
    },
    /// Replace target rows by a_i - alpha*a_ref, using a bounded activity
    /// variable t = a_ref*x when the reference is not an equality.
    RowCombination {
        reference: usize,
        targets: Entries,
        activity: Option<(usize, Entries)>,
    },
    Fixed {
        column: usize,
        value: f64,
        gradient: Gradient,
        entries: Entries,
    },
    /// Substitution in the objective and every other row is one reversible
    /// transaction. `retained` means the equality became the column's bounds.
    Substituted {
        column: usize,
        equation: Arc<Equation>,
        gradient: Gradient,
        other_rows: Entries,
        retained: bool,
    },
    DependentRow {
        row: usize,
        coefficients: Entries,
    },
    /// A_keep = ratio * A_removed. Transfer any multiplier still on a
    /// redundant parallel row before discarding its warm-start coordinate.
    MergedRow {
        keep: usize,
        removed: usize,
        ratio: f64,
    },
    DeletedRow(usize),
    TightenedBound {
        column: usize,
        equation: Arc<Equation>,
        side: Side,
        old: f64,
    },
    TightenedRow {
        row: usize,
        source: usize,
        /// A_row = ratio * A_source.
        ratio: f64,
        side: Side,
    },
    ParallelColumns {
        keep: usize,
        removed: usize,
        ratio: f64,
        keep_bounds: Bounds,
        removed_bounds: Bounds,
    },
    /// A flat objective permits movement indefinitely in one direction.
    /// Save the removed rows to recover a finite feasible value afterwards.
    Unlocked {
        column: usize,
        bounds: Bounds,
        rows: Vec<Equation>,
    },
}

/// Infeasibility proof or recession ray in stable working coordinates.
#[derive(Debug)]
pub(crate) struct Certificate {
    pub mode: Recovery,
    pub point: Point,
}

#[derive(Clone, Debug)]
pub(crate) struct Point {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub z: Vec<f64>,
}

impl Point {
    pub fn zeros(columns: usize, rows: usize) -> Self {
        Self {
            x: vec![0.0; columns],
            y: vec![0.0; rows],
            z: vec![0.0; columns],
        }
    }
}

#[derive(Clone, Debug, Default)]
pub(crate) struct RecoveryTape {
    pub rules: Vec<Rule>,
}

fn active(x: f64, bound: f64) -> bool {
    bound.is_finite() && (x - bound).abs() <= 1e-8 * (1.0 + x.abs().max(bound.abs()))
}

impl RecoveryTape {
    /// Indices are stable model IDs; auxiliary columns are appended. The
    /// adapter packs coordinates, so this tape also recovers early presolve rays.
    pub fn recover(&self, point: &mut Point, mode: Recovery) {
        self.recover_with_slacks(point, mode, &mut []);
    }
    pub fn recover_with_slacks(&self, point: &mut Point, mode: Recovery, slacks: &mut [f64]) {
        for rule in self.rules.iter().rev() {
            match rule {
                Rule::ConeSlack { row, rhs, entries } => {
                    if !mode.primal() || slacks.is_empty() {
                        continue;
                    }
                    slacks[*row] = mode.offset(*rhs)
                        - entries.iter().map(|&(j, a)| a * point.x[j]).sum::<f64>();
                }
                Rule::PsdZeroFace { rows, order } => {
                    if !mode.dual() {
                        continue;
                    }
                    for j in 0..*order {
                        point.y[rows[j * (j + 1) / 2 + j]] = 0.;
                    }
                    let mut k = 0;
                    for j in 0..*order {
                        for i in 0..=j {
                            if i != j {
                                let w = point.y[rows[k]].abs() * std::f64::consts::FRAC_1_SQRT_2;
                                point.y[rows[j * (j + 1) / 2 + j]] -= w;
                                point.y[rows[i * (i + 1) / 2 + i]] -= w;
                            }
                            k += 1;
                        }
                    }
                }
                Rule::SocToLinear { head, tail } => {
                    if !mode.dual() {
                        continue;
                    }
                    let a = point.y[*head];
                    let b = point.y[*tail];
                    point.y[*head] = a + b;
                    point.y[*tail] = a - b;
                }
                Rule::SocFace { head, tail } => {
                    if !mode.dual() {
                        continue;
                    }
                    point.y[*head] = -tail.iter().fold(0.0_f64, |norm, &i| norm.hypot(point.y[i]));
                }
                Rule::RowCombination {
                    reference, targets, ..
                } => {
                    if !mode.dual() {
                        continue;
                    }
                    for &(i, alpha) in targets {
                        point.y[*reference] -= alpha * point.y[i];
                    }
                }
                Rule::Fixed {
                    column,
                    value,
                    gradient,
                    entries,
                } => {
                    if mode.primal() {
                        point.x[*column] = mode.offset(*value);
                    }
                    if mode.dual() {
                        let g = if mode == Recovery::Solution {
                            gradient.evaluate(&point.x)
                        } else {
                            0.0
                        };
                        point.z[*column] =
                            g - entries.iter().map(|&(i, a)| a * point.y[i]).sum::<f64>();
                    }
                }
                Rule::Substituted {
                    column,
                    equation,
                    gradient,
                    other_rows,
                    retained,
                } => {
                    let a = equation.coefficient(*column);
                    if mode.primal() {
                        point.x[*column] = (mode.offset(equation.bounds.lower)
                            - equation.activity(&point.x, Some(*column)))
                            / a;
                    }
                    if mode.dual() {
                        let g = if mode == Recovery::Solution {
                            gradient.evaluate(&point.x)
                        } else {
                            0.0
                        };
                        let y = if *retained {
                            point.y[equation.row]
                        } else {
                            0.0
                        };
                        point.z[*column] = -a * y;
                        point.y[equation.row] = y
                            + (g - other_rows.iter().map(|&(i, v)| v * point.y[i]).sum::<f64>())
                                / a;
                    }
                }
                Rule::DependentRow { row, .. } => {
                    if mode.dual() {
                        point.y[*row] = 0.0;
                    }
                }
                Rule::DeletedRow(row) | Rule::MergedRow { removed: row, .. } => {
                    if mode.dual() {
                        point.y[*row] = 0.0;
                    }
                }
                Rule::TightenedBound {
                    column,
                    equation,
                    side,
                    old,
                } => {
                    if !mode.dual() {
                        continue;
                    }
                    let z = point.z[*column];
                    if !side.used(z) {
                        continue;
                    }
                    if mode == Recovery::Solution && active(point.x[*column], *old) {
                        continue;
                    }
                    // Interior-point iterates have nonzero multipliers even
                    // away from an active bound. Transfer those too: dropping
                    // them would lose stationarity in the original problem.
                    let theta = z / equation.coefficient(*column);
                    point.y[equation.row] += theta;
                    // At this point every column of the saved row has been
                    // restored. No sentinel-based omission of other multipliers.
                    for &(j, a) in &equation.entries {
                        point.z[j] -= a * theta;
                    }
                    point.z[*column] = 0.0;
                }
                Rule::TightenedRow {
                    row,
                    source,
                    ratio,
                    side,
                } => {
                    if !mode.dual() || !side.used(point.y[*row]) {
                        continue;
                    }
                    point.y[*source] += ratio * point.y[*row];
                    point.y[*row] = 0.0;
                }
                Rule::ParallelColumns {
                    keep,
                    removed,
                    ratio,
                    keep_bounds,
                    removed_bounds,
                } => {
                    if mode.primal() {
                        let b = mode.bounds(*keep_bounds);
                        let c = mode.bounds(*removed_bounds);
                        let sum = point.x[*keep];
                        let anchor = if b.lower.is_finite() {
                            b.lower
                        } else if b.upper.is_finite() {
                            b.upper
                        } else {
                            0.0
                        };
                        let other = ((sum - anchor) / ratio).max(c.lower).min(c.upper);
                        point.x[*removed] = other;
                        point.x[*keep] = sum - ratio * other;
                    }
                    if mode.dual() {
                        point.z[*removed] = ratio * point.z[*keep];
                    }
                }
                Rule::Unlocked {
                    column,
                    bounds,
                    rows,
                } => {
                    if mode.primal() {
                        let mut b = mode.bounds(*bounds);
                        for equation in rows {
                            let a = equation.coefficient(*column);
                            let residual = equation.activity(&point.x, Some(*column));
                            let row_bounds = mode.bounds(equation.bounds);
                            let (l, u) = if a > 0.0 {
                                (row_bounds.lower, row_bounds.upper)
                            } else {
                                (row_bounds.upper, row_bounds.lower)
                            };
                            b.lower = b.lower.max((l - residual) / a);
                            b.upper = b.upper.min((u - residual) / a);
                        }
                        point.x[*column] = 0.0_f64.max(b.lower).min(b.upper);
                    }
                    if mode.dual() {
                        point.z[*column] = 0.0;
                        for equation in rows {
                            point.y[equation.row] = 0.0;
                        }
                    }
                }
            }
        }
    }

    /// Forward mapping of caller-supplied native primal and dual coordinates.
    /// These transformations preserve exact stationary starts where possible;
    /// the caller decides whether to use the resulting warm start.
    pub fn reduce_point(&self, point: &mut Point) {
        for rule in &self.rules {
            match rule {
                Rule::ConeSlack { .. } => (),
                Rule::PsdZeroFace { rows, order } => {
                    for j in 0..*order {
                        point.y[rows[j * (j + 1) / 2 + j]] = 0.;
                    }
                }
                Rule::SocToLinear { head, tail } => {
                    let a = point.y[*head];
                    let b = point.y[*tail];
                    point.y[*head] = 0.5 * (a + b);
                    point.y[*tail] = 0.5 * (a - b);
                }
                Rule::SocFace { head, .. } => {
                    point.y[*head] = 0.;
                }
                Rule::RowCombination {
                    reference,
                    targets,
                    activity,
                } => {
                    if let Some((j, terms)) = activity {
                        point.x[*j] = terms.iter().map(|&(k, a)| a * point.x[k]).sum();
                        point.z[*j] = point.y[*reference];
                    }
                    for &(i, alpha) in targets {
                        point.y[*reference] += alpha * point.y[i];
                    }
                }
                Rule::Fixed { column, value, .. } => {
                    point.x[*column] = *value;
                    point.z[*column] = 0.0;
                }
                Rule::Substituted {
                    column,
                    equation,
                    retained,
                    ..
                } => {
                    let a = equation.coefficient(*column);
                    let z = point.z[*column];
                    if *retained {
                        point.y[equation.row] = -z / a;
                    } else {
                        point.y[equation.row] = 0.0;
                        for &(j, b) in &equation.entries {
                            if j != *column {
                                point.z[j] -= b * z / a;
                            }
                        }
                    }
                    point.z[*column] = 0.0;
                }
                Rule::TightenedBound {
                    column, equation, ..
                } => {
                    if equation.entries.len() == 1 {
                        point.z[*column] += equation.coefficient(*column) * point.y[equation.row];
                        point.y[equation.row] = 0.0;
                    }
                }
                Rule::DependentRow { row, coefficients } => {
                    for &(j, a) in coefficients {
                        point.y[j] += a * point.y[*row];
                    }
                    point.y[*row] = 0.0;
                }
                Rule::DeletedRow(row) => point.y[*row] = 0.0,
                Rule::MergedRow {
                    keep,
                    removed,
                    ratio,
                } => {
                    point.y[*keep] += point.y[*removed] / ratio;
                    point.y[*removed] = 0.;
                }
                Rule::TightenedRow {
                    row,
                    source,
                    ratio,
                    side,
                } => {
                    let y = point.y[*source] / ratio;
                    if side.used(y) {
                        point.y[*row] += y;
                        point.y[*source] = 0.0;
                    }
                }
                Rule::ParallelColumns {
                    keep,
                    removed,
                    ratio,
                    ..
                } => {
                    point.x[*keep] += ratio * point.x[*removed];
                    point.z[*removed] = 0.0;
                }
                Rule::Unlocked { column, rows, .. } => {
                    point.z[*column] = 0.0;
                    for equation in rows {
                        point.y[equation.row] = 0.0;
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn inactive_bound_multipliers_are_transferred_for_interior_point_iterates() {
        // x0 + 2*x1 >= 5 and x1 <= 2 imply x0 >= 1. Neither
        // bound need be active at an approximate interior-point solution.
        let tape = RecoveryTape {
            rules: vec![
                Rule::TightenedBound {
                    column: 1,
                    side: Side::Upper,
                    old: f64::INFINITY,
                    equation: Arc::new(Equation {
                        row: 1,
                        entries: vec![(1, -2.0)],
                        bounds: Bounds {
                            lower: -4.0,
                            upper: f64::INFINITY,
                        },
                    }),
                },
                Rule::DeletedRow(1),
                Rule::TightenedBound {
                    column: 0,
                    side: Side::Lower,
                    old: f64::NEG_INFINITY,
                    equation: Arc::new(Equation {
                        row: 0,
                        entries: vec![(0, 1.0), (1, 2.0)],
                        bounds: Bounds {
                            lower: 5.0,
                            upper: f64::INFINITY,
                        },
                    }),
                },
            ],
        };
        let mut p = Point {
            x: vec![1.0 + 3e-4, 2.0 - 1e-4],
            y: vec![0.0; 2],
            z: vec![1e-4, 0.0],
        };
        let complementarity = (p.x[0] - 1.0) * p.z[0];
        tape.recover(&mut p, Recovery::Solution);
        assert_eq!(p.z, [0.0, 0.0]);
        assert_eq!(p.y, [1e-4, 1e-4]);
        assert_eq!(p.y[0], 1e-4);
        assert_eq!(2.0 * p.y[0] - 2.0 * p.y[1], 0.0);
        let recovered = (p.x[0] + 2.0 * p.x[1] - 5.0) * p.y[0] + (4.0 - 2.0 * p.x[1]) * p.y[1];
        assert!((recovered - complementarity).abs() < 1e-18);
    }

    #[test]
    fn coupled_substitution_recovers_stationarity_with_a_retained_bound() {
        // x0 + 2*x1 = 5; x0 >= 1. The objective is
        // .5*x^T [[2,1],[1,4]]*x + [0,-10]^T*x.
        // x=(1,2), native y=-.5, z=(4.5,0) satisfies the original KKT system.
        let tape = RecoveryTape {
            rules: vec![Rule::Substituted {
                column: 0,
                equation: Arc::new(Equation {
                    row: 0,
                    entries: vec![(0, 1.0), (1, 2.0)],
                    bounds: Bounds::fixed(5.0),
                }),
                gradient: Gradient {
                    constant: 0.0,
                    terms: vec![(0, 2.0), (1, 1.0)],
                },
                other_rows: vec![],
                retained: true,
            }],
        };
        let mut p = Point {
            x: vec![0.0, 2.0],
            y: vec![-4.5],
            z: vec![0.0, 0.0],
        };
        tape.recover(&mut p, Recovery::Solution);
        assert_eq!(p.x, [1.0, 2.0]);
        assert_eq!(p.y, [-0.5]);
        assert_eq!(p.z, [4.5, 0.0]);
        assert_eq!(2.0 * p.x[0] + p.x[1], p.y[0] + p.z[0]);
        assert_eq!(p.x[0] + 4.0 * p.x[1] - 10.0, 2.0 * p.y[0] + p.z[1]);
    }

    #[test]
    fn propagated_bounds_restore_farkas_multipliers_without_a_primal_point() {
        // x0 + 2*x1 >= 5, x1 <= 2 implies x0 >= 1.
        let tape = RecoveryTape {
            rules: vec![Rule::TightenedBound {
                column: 0,
                side: Side::Lower,
                old: f64::NEG_INFINITY,
                equation: Arc::new(Equation {
                    row: 0,
                    entries: vec![(0, 1.0), (1, 2.0)],
                    bounds: Bounds {
                        lower: 5.0,
                        upper: f64::INFINITY,
                    },
                }),
            }],
        };
        let mut p = Point {
            x: vec![f64::NAN; 2],
            y: vec![0.0],
            z: vec![3.0, 0.0],
        };
        tape.recover(&mut p, Recovery::PrimalInfeasibility);
        assert_eq!(p.y, [3.0]);
        assert_eq!(p.z, [0.0, -6.0]);
    }

    #[test]
    fn aggregation_recovers_solution_and_recession_bounds_for_both_signs() {
        for ratio in [-2.0, 2.0] {
            let bounds = Bounds {
                lower: 1.0,
                upper: f64::INFINITY,
            };
            let tape = RecoveryTape {
                rules: vec![Rule::ParallelColumns {
                    keep: 0,
                    removed: 1,
                    ratio,
                    keep_bounds: bounds,
                    removed_bounds: bounds,
                }],
            };
            let mut p = Point {
                x: vec![3.0 + ratio * 2.0, 0.0],
                y: vec![],
                z: vec![0.0; 2],
            };
            tape.recover(&mut p, Recovery::Solution);
            assert!(p.x.iter().all(|&x| bounds.contains(x)));
            assert_eq!(p.x[0] + ratio * p.x[1], 3.0 + ratio * 2.0);
            let mut p = Point {
                x: vec![ratio, 0.0],
                y: vec![],
                z: vec![0.0; 2],
            };
            tape.recover(&mut p, Recovery::DualInfeasibility);
            assert!(p.x.iter().all(|&x| x >= 0.0));
            assert_eq!(p.x[0] + ratio * p.x[1], ratio);
        }
    }
    #[test]
    fn cone_slacks_are_recovered_before_earlier_variable_aggregation() {
        let tape = RecoveryTape {
            rules: vec![
                Rule::ParallelColumns {
                    keep: 0,
                    removed: 1,
                    ratio: 2.,
                    keep_bounds: Bounds::FREE,
                    removed_bounds: Bounds::FREE,
                },
                Rule::ConeSlack {
                    row: 0,
                    rhs: 5.,
                    entries: vec![(0, 3.)],
                },
            ],
        };
        let mut point = Point {
            x: vec![7., 0.],
            y: vec![0.],
            z: vec![0.; 2],
        };
        let mut slacks = [0.];
        tape.recover_with_slacks(&mut point, Recovery::Solution, &mut slacks);
        assert_eq!(point.x[0] + 2. * point.x[1], 7.);
        assert_eq!(slacks[0], 5. - 3. * point.x[0] - 6. * point.x[1]);
        assert_eq!(slacks, [-16.]);
    }
}
