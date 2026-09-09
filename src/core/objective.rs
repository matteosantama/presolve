// SPDX-License-Identifier: Apache-2.0
// Modified for this library; copyright and attribution notices are in NOTICE.
//! Sparse quadratic objective updates under affine transformations.
//! Every affine substitution applies P' = T^T P T, c' = T^T(c+P d).

use crate::matrix::sparse::{Entries, SparseMatrix};
use std::collections::BTreeMap;

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
pub(crate) struct Objective {
    pub p: SparseMatrix,
    pub c: Vec<f64>,
    pub constant: f64,
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
    /// happen before mutation. The total number of P nonzeros must not increase.
    pub fn substitute(
        &mut self,
        k: usize,
        offset: f64,
        slopes: &[(usize, f64)],
        max_fill: usize,
    ) -> Option<Gradient> {
        assert!(slopes.iter().all(|&(j, _)| j != k));
        let gradient = self.gradient(k);
        let diagonal = self.p.get(k, k);
        let constant = self.constant + self.c[k] * offset + 0.5 * diagonal * offset * offset;
        if !constant.is_finite() || !offset.is_finite() {
            return None;
        }
        let mut linear = BTreeMap::new();
        let mut quadratic = BTreeMap::new();
        let slopes_by_column: BTreeMap<_, _> = slopes.iter().copied().collect();
        let mut added = 0;
        let mut removed = 2 * self.p.column(k).len() - usize::from(diagonal != 0.0);
        let mut update_quadratic = |j: usize, l: usize, check_fill_only: bool| -> Option<()> {
            let (j, l) = (j.min(l), j.max(l));
            let old = self.p.get(j, l);
            // Check missing coefficients before allocating updates to existing
            // entries, so excessive fill can be rejected before building a dense update.
            if check_fill_only && old != 0.0 {
                return Some(());
            }
            if quadratic.contains_key(&(j, l)) {
                return Some(());
            }
            let vj = slopes_by_column.get(&j).copied().unwrap_or(0.0);
            let vl = slopes_by_column.get(&l).copied().unwrap_or(0.0);
            let value = old + self.p.get(j, k) * vl + vj * self.p.get(k, l) + diagonal * vj * vl;
            if !value.is_finite() {
                return None;
            }
            let count = if j == l { 1 } else { 2 };
            if old == 0.0 && value != 0.0 {
                added += count;
                if added > max_fill {
                    return None;
                }
            } else if old != 0.0 && value == 0.0 {
                removed += count;
            }
            quadratic.insert((j, l), value);
            if j != l {
                quadratic.insert((l, j), value);
            }
            Some(())
        };
        for &(j, pjk) in self.p.column(k) {
            if j == k {
                continue;
            }
            *linear.entry(j).or_insert(self.c[j]) += offset * pjk;
        }
        for &(j, vj) in slopes {
            *linear.entry(j).or_insert(self.c[j]) += (self.c[k] + diagonal * offset) * vj;
        }
        for check_fill_only in [true, false] {
            for &(j, _) in self.p.column(k) {
                if j != k {
                    for &(l, _) in slopes {
                        update_quadratic(j, l, check_fill_only)?;
                    }
                }
            }
            if diagonal != 0.0 {
                for &(j, _) in slopes {
                    for &(l, _) in slopes {
                        update_quadratic(j, l, check_fill_only)?;
                    }
                }
            }
        }
        if added > removed {
            return None;
        }
        if linear
            .values()
            .chain(quadratic.values())
            .any(|v| !v.is_finite())
        {
            return None;
        }
        self.p.remove_row(k);
        self.p.remove_column(k);
        self.c[k] = 0.0;
        for (j, c) in linear {
            self.c[j] = c;
        }
        for ((j, l), p) in quadratic {
            self.p.set(j, l, p);
        }
        self.constant = constant;
        Some(gradient)
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
        self.p.remove_row(column);
        self.p.remove_column(column);
        self.c[column] = 0.0;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wide_substitution_skips_zero_quadratic_work_and_rejects_fill_before_constructing_it() {
        let n = 4097;
        let slopes: Vec<_> = (1..n).map(|j| (j, 1.0)).collect();
        let mut objective = Objective {
            p: SparseMatrix::zeros(n, n),
            c: vec![0.0; n],
            constant: 0.0,
        };
        objective.c[0] = 2.0;
        assert!(objective.substitute(0, 3.0, &slopes, 0).is_some());
        assert_eq!(objective.constant, 6.0);
        assert_eq!(objective.p.nnz(), 0);
        assert!(objective.c[1..].iter().all(|&v| v == 2.0));
        objective.p.set(0, 0, 1.0);
        assert!(objective.substitute(0, 0.0, &slopes, 64).is_none());
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
                p: SparseMatrix::from_matrix(&crate::matrix::test_matrix(3, 3, entries).unwrap()),
                c: vec![2.0, 3.0, 4.0],
                constant: 5.0,
            };
            let original = objective.clone();
            assert!(
                objective
                    .substitute(0, 2.0, &[(1, -1.0), (2, -1.0)], 64)
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
                p: SparseMatrix::from_matrix(&crate::matrix::test_matrix(3, 3, entries).unwrap()),
                c: vec![2.0, 3.0, 4.0],
                constant: 5.0,
            };
            let original = objective.clone();
            assert!(objective.substitute(0, 2.0, &slopes, 64).is_some());
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
            p: SparseMatrix::zeros(4, 4),
            c: vec![0.0; 4],
            constant: 0.0,
        };
        for i in 0..2 {
            for j in 0..2 {
                objective.p.set(i, j, 1.0);
            }
        }
        assert!(
            objective
                .substitute(0, 0.0, &[(1, -1.0), (2, 1.0), (3, 1.0)], 64)
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
            p: SparseMatrix::from_matrix(
                &crate::matrix::test_matrix(
                    3,
                    3,
                    vec![4.0, -1.0, -1.0, -1.0, 2.0, 0.0, -1.0, 0.0, 2.0],
                )
                .unwrap(),
            ),
            c: vec![2.0, 3.0, 4.0],
            constant: 5.0,
        };
        let original = objective.clone();
        assert!(
            objective
                .substitute(0, 2.0, &[(1, 0.5), (2, 0.5)], 64)
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
