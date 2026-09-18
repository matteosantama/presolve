// SPDX-License-Identifier: Apache-2.0
// Modified for this library; copyright and attribution notices are in NOTICE.
//! Sparse quadratic objective updates under affine transformations.
//! Every affine substitution applies P' = T^T P T, c' = T^T(c+P d).

use crate::matrix::sparse::{Entries, SymmetricMatrix};
use std::time::Instant;

#[derive(Clone, Copy, Debug)]
pub(crate) enum SubstitutionFailure {
    Numerical,
    ConstraintFill,
    QuadraticFill,
    HessianGrowth,
    Deadline,
}

#[derive(Clone, Debug)]
pub(crate) struct Gradient {
    pub constant: f64,
    pub terms: Entries,
}

impl Gradient {
    pub fn evaluate(&self, x: &[f64]) -> f64 {
        self.terms
            .iter()
            .fold(self.constant, |sum, &(j, a)| sum + a * x[j])
    }
}

#[derive(Clone, Debug)]
struct Affected {
    column: usize,
    curvature: f64,
    slope: f64,
    has_slope: bool,
}

#[derive(Clone, Debug, Default)]
pub(crate) struct Scratch {
    affected: Vec<Affected>,
    curvature_indices: Vec<usize>,
    slope_indices: Vec<usize>,
    linear: Entries,
    quadratic: Vec<(usize, usize, f64)>,
}
impl Scratch {
    fn clear(&mut self) {
        self.affected.clear();
        self.curvature_indices.clear();
        self.slope_indices.clear();
        self.linear.clear();
        self.quadratic.clear();
    }
}

#[derive(Clone, Debug)]
pub(crate) struct Objective {
    pub p: SymmetricMatrix,
    pub c: Vec<f64>,
    pub constant: f64,
    pub scratch: Scratch,
}

impl Objective {
    pub fn gradient(&self, column: usize) -> Gradient {
        Gradient {
            constant: self.c[column],
            terms: self.p.row(column).to_vec(),
        }
    }

    pub fn diagonal(&self, column: usize) -> Option<f64> {
        self.p
            .row(column)
            .iter()
            .all(|&(j, _)| j == column)
            .then(|| self.p.get(column, column))
    }

    /// Eliminate x_k = offset + sum slope_j*x_j. All arithmetic and fill checks
    /// happen before mutation. Hessian growth is controlled independently of fill.
    pub fn substitute(
        &mut self,
        k: usize,
        offset: f64,
        slopes: &[(usize, f64)],
        max_fill: usize,
        allow_growth: bool,
        deadline: Option<Instant>,
    ) -> Option<Gradient> {
        self.try_substitute(k, offset, slopes, max_fill, allow_growth, deadline)
            .ok()
    }

    pub fn try_substitute(
        &mut self,
        k: usize,
        offset: f64,
        slopes: &[(usize, f64)],
        max_fill: usize,
        allow_growth: bool,
        deadline: Option<Instant>,
    ) -> Result<Gradient, SubstitutionFailure> {
        assert!(slopes.iter().all(|&(j, _)| j != k));
        debug_assert!(slopes.windows(2).all(|pair| pair[0].0 < pair[1].0));
        let mut scratch = std::mem::take(&mut self.scratch);
        scratch.clear();
        let result = (|| {
            let gradient = self.gradient(k);
            let diagonal = self.p.get(k, k);
            let constant = self.constant + self.c[k] * offset + 0.5 * diagonal * offset * offset;
            if !constant.is_finite() || !offset.is_finite() {
                return Err(SubstitutionFailure::Numerical);
            }
            let curved = !self.p.column(k).is_empty();
            // Merge the two sorted supports once. The LP path stages only its
            // linear coefficients, without constructing quadratic scratch data.
            let mut p = self
                .p
                .column(k)
                .iter()
                .copied()
                .filter(|&(j, _)| j != k)
                .peekable();
            let mut slopes = slopes.iter().copied().peekable();
            while p.peek().is_some() || slopes.peek().is_some() {
                let pj = p.peek().map_or(usize::MAX, |e| e.0);
                let sj = slopes.peek().map_or(usize::MAX, |e| e.0);
                let j = pj.min(sj);
                let curvature = if pj == j { p.next().unwrap().1 } else { 0.0 };
                let slope = if sj == j {
                    slopes.next().unwrap().1
                } else {
                    0.0
                };
                let mut linear = self.c[j];
                if pj == j {
                    linear += offset * curvature;
                }
                if sj == j {
                    linear += (self.c[k] + diagonal * offset) * slope;
                }
                scratch.linear.push((j, linear));
                if curved {
                    let at = scratch.affected.len();
                    scratch.affected.push(Affected {
                        column: j,
                        curvature,
                        slope,
                        has_slope: sj == j,
                    });
                    if pj == j {
                        scratch.curvature_indices.push(at);
                    }
                    if sj == j {
                        scratch.slope_indices.push(at);
                    }
                }
            }
            let mut added = 0;
            let mut removed = 2 * self.p.column(k).len() - usize::from(diagonal != 0.0);
            let mut visits = 0usize;
            let mut update =
                |j: usize, l: usize, missing: bool| -> Result<(), SubstitutionFailure> {
                    visits += 1;
                    if visits.is_multiple_of(1024) && deadline.is_some_and(|d| Instant::now() >= d)
                    {
                        return Err(SubstitutionFailure::Deadline);
                    }
                    let a = &scratch.affected[j];
                    let b = &scratch.affected[l];
                    let old = self.p.get(a.column, b.column);
                    // Each unordered candidate occurs once in each pass. Reject
                    // excessive fill before staging updates to existing positions.
                    if (old == 0.0) != missing {
                        return Ok(());
                    }
                    let value = old
                        + a.curvature * b.slope
                        + a.slope * b.curvature
                        + diagonal * a.slope * b.slope;
                    if !value.is_finite() {
                        return Err(SubstitutionFailure::Numerical);
                    }
                    let count = if j == l { 1 } else { 2 };
                    if old == 0.0 && value != 0.0 {
                        added += count;
                        if added > max_fill {
                            return Err(SubstitutionFailure::QuadraticFill);
                        }
                    } else if old != 0.0 && value == 0.0 {
                        removed += count;
                    }
                    if value != old {
                        scratch.quadratic.push((a.column, b.column, value));
                    }
                    Ok(())
                };
            for missing in [true, false] {
                for (j, a) in scratch.affected.iter().enumerate() {
                    let use_slopes = a.curvature != 0.0 || (diagonal != 0.0 && a.has_slope);
                    let use_curvature = a.has_slope;
                    if use_slopes && use_curvature {
                        for l in j..scratch.affected.len() {
                            update(j, l, missing)?;
                        }
                    } else {
                        let indices = if use_slopes {
                            &scratch.slope_indices
                        } else if use_curvature {
                            &scratch.curvature_indices
                        } else {
                            continue;
                        };
                        let start = indices.partition_point(|&l| l < j);
                        for &l in &indices[start..] {
                            update(j, l, missing)?;
                        }
                    }
                }
            }
            if !allow_growth && added > removed {
                return Err(SubstitutionFailure::HessianGrowth);
            }
            if scratch.linear.iter().any(|e| !e.1.is_finite()) {
                return Err(SubstitutionFailure::Numerical);
            }
            self.p.remove_variable(k);
            self.c[k] = 0.0;
            for &(j, c) in &scratch.linear {
                self.c[j] = c;
            }
            for &(j, l, p) in &scratch.quadratic {
                self.p.set(j, l, p);
            }
            self.constant = constant;
            Ok(gradient)
        })();
        self.scratch = scratch;
        result
    }

    /// Parallel-column aggregation is exact for a QP only when its objective
    /// also depends on x_j and x_k through x_j + ratio*x_k.
    pub fn parallel(&self, j: usize, k: usize, ratio: f64, tolerance: f64) -> bool {
        let a = self.c[k];
        let b = ratio * self.c[j];
        a.is_finite()
            && b.is_finite()
            && (a - b).abs() <= tolerance * a.abs().max(b.abs())
            && self.parallel_curvature(j, k, ratio, tolerance)
    }

    /// P*(-ratio*e_j + e_k) = 0 makes transfers along parallel A columns
    /// linear in the objective, even if the individual variables have curvature.
    pub fn parallel_curvature(&self, j: usize, k: usize, ratio: f64, tolerance: f64) -> bool {
        let near = |a: f64, b: f64| {
            a.is_finite() && b.is_finite() && (a - b).abs() <= tolerance * a.abs().max(b.abs())
        };
        self.p
            .column(j)
            .iter()
            .all(|&(r, a)| near(self.p.get(r, k), ratio * a))
            && self
                .p
                .column(k)
                .iter()
                .all(|&(r, a)| near(a, ratio * self.p.get(r, j)))
    }

    pub fn aggregate(&mut self, column: usize) {
        self.p.remove_variable(column);
        self.c[column] = 0.0;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[allow(clippy::needless_range_loop)] // Dense reference algebra uses explicit matrix indices.
    fn packed_updates_match_dense_affine_transformations_and_clear_failed_scratch() {
        let n = 6;
        for seed in 0..24 {
            let k = seed % n;
            let dense: Vec<f64> = (0..n)
                .flat_map(|i| {
                    (0..n).map(move |j| {
                        (0..3)
                            .map(|r| {
                                (((i + r + seed) % 5) as f64 - 2.)
                                    * (((j + r + seed) % 5) as f64 - 2.)
                            })
                            .sum()
                    })
                })
                .collect();
            let mut objective = Objective {
                p: SymmetricMatrix::from_matrix(
                    &crate::matrix::test_matrix(n, n, dense.clone()).unwrap(),
                ),
                c: (0..n).map(|j| j as f64 - 2.).collect(),
                constant: 3.,
                scratch: Default::default(),
            };
            let original = objective.clone();
            let slopes: Vec<_> = (0..n)
                .filter(|&j| j != k && (j + seed) % 3 != 0)
                .map(|j| (j, ((j + seed) % 5) as f64 / 2. - 1.))
                .collect();
            let mut transform = vec![vec![0.; n]; n];
            for j in 0..n {
                if j != k {
                    transform[j][j] = 1.;
                }
            }
            for &(j, a) in &slopes {
                transform[k][j] = a;
            }
            objective
                .try_substitute(k, 2., &slopes, usize::MAX, true, None)
                .unwrap();
            for j in 0..n {
                let c: f64 = (0..n)
                    .map(|a| transform[a][j] * (original.c[a] + 2. * dense[a * n + k]))
                    .sum();
                assert_eq!(objective.c[j], c);
                for l in 0..n {
                    let mut expected = 0.;
                    for a in 0..n {
                        for b in 0..n {
                            expected += transform[a][j] * dense[a * n + b] * transform[b][l];
                        }
                    }
                    assert_eq!(objective.p.get(j, l), expected);
                }
            }
        }
        let mut objective = Objective {
            p: SymmetricMatrix::zeros(n),
            c: vec![1.; n],
            constant: 0.,
            scratch: Default::default(),
        };
        objective.p.set(0, 0, 1.);
        let slopes: Vec<_> = (1..n).map(|j| (j, 1.)).collect();
        assert!(matches!(
            objective.try_substitute(0, 2., &slopes, 0, true, None),
            Err(SubstitutionFailure::QuadraticFill)
        ));
        assert_eq!(objective.p.nnz(), 1);
        objective
            .try_substitute(0, 2., &slopes, usize::MAX, true, None)
            .unwrap();
        for j in 1..n {
            for l in 1..n {
                assert_eq!(objective.p.get(j, l), 1.);
            }
        }
        assert_eq!(objective.c[1..], [4.; 5]);
        assert_eq!(objective.constant, 4.);
    }

    #[test]
    fn wide_substitution_skips_zero_quadratic_work_and_rejects_fill_before_constructing_it() {
        let n = 4097;
        let slopes: Vec<_> = (1..n).map(|j| (j, 1.0)).collect();
        let mut objective = Objective {
            p: SymmetricMatrix::zeros(n),
            c: vec![0.0; n],
            constant: 0.0,
            scratch: Default::default(),
        };
        objective.c[0] = 2.0;
        assert!(
            objective
                .substitute(0, 3.0, &slopes, 0, false, None)
                .is_some()
        );
        assert_eq!(objective.constant, 6.0);
        assert_eq!(objective.p.nnz(), 0);
        assert!(objective.c[1..].iter().all(|&v| v == 2.0));
        objective.p.set(0, 0, 1.0);
        assert!(
            objective
                .substitute(0, 0.0, &slopes, 64, false, None)
                .is_none()
        );
        assert_eq!(objective.p.nnz(), 1);
    }

    fn value(o: &Objective, x: &[f64]) -> f64 {
        o.constant
            + o.c.iter().zip(x).map(|(c, x)| c * x).sum::<f64>()
            + 0.5
                * (0..x.len())
                    .flat_map(|i| o.p.row(i).iter().map(move |&(j, a)| a * x[i] * x[j]))
                    .sum::<f64>()
    }

    #[test]
    fn substitution_rejects_net_p_growth_transactionally() {
        for entries in [
            vec![1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
            // Existing diagonals do not pay for the new off-diagonal pair.
            vec![1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
        ] {
            let mut objective = Objective {
                p: SymmetricMatrix::from_matrix(
                    &crate::matrix::test_matrix(3, 3, entries).unwrap(),
                ),
                c: vec![2.0, 3.0, 4.0],
                constant: 5.0,
                scratch: Default::default(),
            };
            let original = objective.clone();
            assert!(
                objective
                    .substitute(0, 2.0, &[(1, -1.0), (2, -1.0)], 64, false, None)
                    .is_none()
            );
            assert_eq!(objective.c, original.c);
            assert_eq!(objective.constant, original.constant);
            for i in 0..3 {
                assert_eq!(objective.p.row(i), original.p.row(i));
                assert_eq!(objective.p.column(i), original.p.column(i));
            }
        }
    }

    #[test]
    fn permitted_growth_preserves_the_objective_and_deadlines_are_transactional() {
        let mut objective = Objective {
            p: SymmetricMatrix::from_matrix(
                &crate::matrix::test_matrix(3, 3, vec![2., 1., 0., 1., 2., 0., 0., 0., 1.])
                    .unwrap(),
            ),
            c: vec![2., 3., 4.],
            constant: 5.,
            scratch: Default::default(),
        };
        let original = objective.clone();
        assert!(
            objective
                .substitute(0, 2., &[(1, -1.), (2, -1.)], usize::MAX, true, None)
                .is_some()
        );
        for y in [-2., 0., 3.] {
            for z in [-3., 0., 2.] {
                assert!(
                    (value(&original, &[2. - y - z, y, z]) - value(&objective, &[0., y, z])).abs()
                        < 1e-10
                );
            }
        }
        let mut wide = Objective {
            p: SymmetricMatrix::zeros(65),
            c: vec![1.; 65],
            constant: 7.,
            scratch: Default::default(),
        };
        wide.p.set(0, 0, 1.);
        let slopes: Vec<_> = (1..65).map(|j| (j, 1.)).collect();
        assert!(
            wide.substitute(0, 2., &slopes, usize::MAX, true, Some(Instant::now()))
                .is_none()
        );
        assert_eq!(wide.p.nnz(), 1);
        assert_eq!(wide.p.get(0, 0), 1.);
        assert_eq!(wide.c, vec![1.; 65]);
        assert_eq!(wide.constant, 7.);
    }

    #[test]
    fn substitution_allows_new_positions_without_net_p_growth() {
        for (entries, slopes, expected_nnz) in [
            // Move a diagonal entry without changing nnz.
            (
                vec![1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                vec![(1, 1.0)],
                1,
            ),
            // A new off-diagonal pair is paid for by removing the pivot row/column.
            (
                vec![4.0, 1.0, 1.0, 1.0, 2.0, 0.0, 1.0, 0.0, 2.0],
                vec![(1, -1.0), (2, -1.0)],
                4,
            ),
        ] {
            let mut objective = Objective {
                p: SymmetricMatrix::from_matrix(
                    &crate::matrix::test_matrix(3, 3, entries).unwrap(),
                ),
                c: vec![2.0, 3.0, 4.0],
                constant: 5.0,
                scratch: Default::default(),
            };
            let original = objective.clone();
            assert!(
                objective
                    .substitute(0, 2.0, &slopes, 64, false, None)
                    .is_some()
            );
            assert_eq!(objective.p.nnz(), expected_nnz);
            assert!(objective.p.nnz() <= original.p.nnz());
            for y in [-2.0, 0.0, 3.0] {
                for z in [-1.0, 0.0, 4.0] {
                    let mut x = [0.0, y, z];
                    x[0] = 2.0 + slopes.iter().map(|&(j, v)| v * x[j]).sum::<f64>();
                    assert_eq!(value(&original, &x), value(&objective, &[0.0, y, z]));
                }
            }
        }
    }

    #[test]
    fn cancellation_of_existing_entries_can_pay_for_new_entries() {
        // (x0 + x1)^2 becomes (x2 + x3)^2. Removing the pivot accounts for
        // three entries; cancelling the remaining x1 diagonal pays for the fourth.
        let mut objective = Objective {
            p: SymmetricMatrix::zeros(4),
            c: vec![0.0; 4],
            constant: 0.0,
            scratch: Default::default(),
        };
        for i in 0..2 {
            for j in 0..2 {
                objective.p.set(i, j, 1.0);
            }
        }
        assert!(
            objective
                .substitute(0, 0.0, &[(1, -1.0), (2, 1.0), (3, 1.0)], 64, false, None)
                .is_some()
        );
        assert_eq!(objective.p.nnz(), 4);
        assert!(objective.p.row(1).is_empty());
        for i in 2..4 {
            for j in 2..4 {
                assert_eq!(objective.p.get(i, j), 1.0);
            }
        }
    }

    #[test]
    fn substitution_allows_exact_cancellation_at_a_missing_entry() {
        let mut objective = Objective {
            p: SymmetricMatrix::from_matrix(
                &crate::matrix::test_matrix(
                    3,
                    3,
                    vec![4.0, -1.0, -1.0, -1.0, 2.0, 0.0, -1.0, 0.0, 2.0],
                )
                .unwrap(),
            ),
            c: vec![2.0, 3.0, 4.0],
            constant: 5.0,
            scratch: Default::default(),
        };
        let original = objective.clone();
        assert!(
            objective
                .substitute(0, 2.0, &[(1, 0.5), (2, 0.5)], 64, false, None)
                .is_some()
        );
        assert_eq!(objective.p.nnz(), 2);
        assert_eq!(objective.p.get(1, 2), 0.0);
        for y in [-2.0, 0.0, 3.0] {
            for z in [-1.0, 0.0, 4.0] {
                assert_eq!(
                    value(&original, &[2.0 + 0.5 * y + 0.5 * z, y, z]),
                    value(&objective, &[0.0, y, z])
                );
            }
        }
    }
}
