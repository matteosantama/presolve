//! Exact dependencies among balanced incidence equalities.
//!
//! In a qualifying component each column occurs in exactly two eligible rows,
//! with opposite coefficients. A row is redundant when the right-hand sides
//! also sum exactly to zero. The homogeneous search runs first, preserving
//! dependencies that overlapping nonzero equalities could obscure.
use crate::{
    model::{Model, RowDomain, tape::Record},
    result::RuleId,
};
use std::time::Instant;
fn find(parent: &mut [usize], mut i: usize) -> usize {
    while parent[i] != i {
        parent[i] = parent[parent[i]];
        i = parent[i];
    }
    i
}
impl Model {
    pub(super) fn network_equalities(&mut self, deadline: Instant) -> usize {
        self.enter(RuleId::NetworkEqualities);
        let mut eligible = Vec::new();
        let (mut positive, mut negative) = (false, false);
        for (i, r) in self.rows.iter().enumerate() {
            if i % 1024 == 0 && Instant::now() >= deadline {
                return 0;
            }
            if matches!(r,RowDomain::Linear(b) if b.equality()) && !self.a.row(i).is_empty() {
                eligible.push(i);
                let RowDomain::Linear(bounds) = r else {
                    unreachable!()
                };
                positive |= bounds.lower > 0.0;
                negative |= bounds.lower < 0.0;
            }
        }
        // A zero sum containing a nonzero RHS requires both signs. If either
        // is absent, only the existing homogeneous search can remove a row.
        if !positive || !negative {
            if positive || negative {
                eligible
                    .retain(|&i| matches!(self.rows[i], RowDomain::Linear(b) if b.lower == 0.0));
            }
            return self.incidence_dependencies(eligible, deadline, false);
        }
        let homogeneous = eligible
            .iter()
            .copied()
            .filter(|&i| matches!(self.rows[i], RowDomain::Linear(b) if b.lower == 0.0))
            .collect();
        let mut removed = self.incidence_dependencies(homogeneous, deadline, false);
        if Instant::now() < deadline {
            eligible.retain(|&i| !self.a.row(i).is_empty());
            removed += self.incidence_dependencies(eligible, deadline, true);
        }
        removed
    }

    fn incidence_dependencies(
        &mut self,
        eligible: Vec<usize>,
        deadline: Instant,
        check_rhs: bool,
    ) -> usize {
        if eligible.len() < 2 {
            return 0;
        }
        // Most negative cases need only the domain scan. Allocate graph scratch
        // after finding enough rows for a possible dependence.
        let mut parent = vec![usize::MAX; self.rows.len()];
        let mut sizes = vec![1usize; self.rows.len()];
        for &i in &eligible {
            parent[i] = i;
        }
        let mut bad = vec![false; self.rows.len()];
        for j in 0..self.bounds.len() {
            if j % 256 == 0 && Instant::now() >= deadline {
                return 0;
            }
            let mut incidence = self
                .a
                .column(j)
                .iter()
                .filter(|&(i, _)| parent[i] != usize::MAX);
            let Some((i, a)) = incidence.next() else {
                continue;
            };
            let Some((k, b)) = incidence.next() else {
                bad[i] = true;
                continue;
            };
            if let Some((l, _)) = incidence.next() {
                bad[i] = true;
                bad[k] = true;
                bad[l] = true;
                for (r, _) in incidence {
                    bad[r] = true;
                }
                continue;
            }
            if a != -b {
                bad[i] = true;
                bad[k] = true;
                continue;
            }
            let mut x = find(&mut parent, i);
            let mut y = find(&mut parent, k);
            if x != y {
                if sizes[x] < sizes[y] {
                    std::mem::swap(&mut x, &mut y);
                }
                parent[y] = x;
                sizes[x] += sizes[y];
            }
        }
        let mut invalid = vec![false; self.rows.len()];
        for &i in &eligible {
            let r = find(&mut parent, i);
            invalid[r] |= bad[i];
        }
        let mut groups = std::collections::BTreeMap::<usize, Vec<usize>>::new();
        for i in eligible {
            let r = find(&mut parent, i);
            if !invalid[r] {
                groups.entry(r).or_default().push(i);
            }
        }
        let mut removed = 0;
        for group in groups.into_values() {
            if Instant::now() >= deadline {
                break;
            }
            if group.len() < 2 {
                continue;
            }
            // Certify every addition, not just a rounded final sum. This rejects
            // rounded-to-zero imbalance and conservatively skips uncertain sums.
            if check_rhs {
                let rhs = group.iter().try_fold(0.0, |sum, &i| {
                    let RowDomain::Linear(bounds) = self.rows[i] else {
                        unreachable!()
                    };
                    super::dependencies::exact_difference(sum, -bounds.lower)
                });
                if rhs != Some(0.0) {
                    continue;
                }
            }
            // Every column and the right-hand sides cancel exactly.
            // Drop the longest row to remove the most constraint entries.
            let row = *group
                .iter()
                .max_by_key(|&&i| (self.a.row(i).len(), std::cmp::Reverse(i)))
                .unwrap();
            let coefficients = group
                .into_iter()
                .filter(|&i| i != row)
                .map(|i| (i, -1.))
                .collect();
            self.clear_row(row);
            self.record(Record::DependentRow { row, coefficients });
            removed += 1;
        }
        removed
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        matrix::{linked::LinkedMatrix, sparse::SymmetricMatrix},
        model::{
            objective::Objective,
            tape::{Point, Recovery},
        },
        problem::Bounds,
    };
    use std::time::Duration;
    fn model(a: &[Vec<f64>], rhs: &[f64]) -> Model {
        let n = a[0].len();
        Model::from_parts(
            LinkedMatrix::from_columns(a.len(), n, |j| {
                a.iter()
                    .enumerate()
                    .filter_map(move |(i, r)| (r[j] != 0.).then_some((i, r[j])))
            }),
            Objective {
                p: SymmetricMatrix::zeros(n),
                c: vec![0.; n],
                constant: 0.,
                scratch: Default::default(),
            },
            rhs.iter()
                .map(|&b| RowDomain::Linear(Bounds::fixed(b)))
                .collect(),
            vec![Bounds::FREE; n],
        )
    }
    fn cycle() -> Vec<Vec<f64>> {
        vec![vec![1., 0., -1.], vec![-1., 1., 0.], vec![0., -1., 1.]]
    }
    #[test]
    fn incidence_removal_preserves_primal_and_warm_dual_gradient() {
        let a = cycle();
        let mut m = model(&a, &[0.; 3]);
        let deadline = Instant::now() + Duration::from_secs(1);
        assert_eq!(m.network_equalities(deadline), 1);
        assert_eq!(m.a.nnz(), 4);
        let original = Point {
            x: vec![2.; 3],
            y: vec![2., 3., 5.],
            z: vec![0.; 3],
        };
        let mut p = original.clone();
        m.postsolve.reduce_point(&mut p);
        m.postsolve.recover(&mut p, Recovery::Solution);
        assert_eq!(p.x, original.x);
        for (j, _) in a[0].iter().enumerate() {
            let gradient = |y: &[f64]| {
                a.iter()
                    .zip(y)
                    .map(|(row, multiplier)| row[j] * multiplier)
                    .sum::<f64>()
            };
            assert_eq!(gradient(&p.y), gradient(&original.y));
        }
    }
    #[test]
    fn balanced_nonzero_rhs_is_removed_but_rounded_balance_is_not() {
        let mut m = model(&cycle(), &[3., -1., -2.]);
        assert_eq!(
            m.network_equalities(Instant::now() + Duration::from_secs(1)),
            1
        );
        let mut m = model(&cycle(), &[1e16, 1., -1e16]);
        assert_eq!(
            m.network_equalities(Instant::now() + Duration::from_secs(1)),
            0
        );
    }
    #[test]
    fn homogeneous_dependencies_survive_overlapping_nonzero_equations() {
        let mut a = cycle();
        a.extend([vec![1., 0., 0.], vec![0., -1., 0.]]);
        let mut m = model(&a, &[0., 0., 0., 1., -1.]);
        assert_eq!(
            m.network_equalities(Instant::now() + Duration::from_secs(1)),
            1
        );
        assert!(m.a.row(0).is_empty());
        assert!(!m.a.row(3).is_empty());
        assert!(!m.a.row(4).is_empty());
    }
    #[test]
    fn mixed_rhs_balance_and_nonzero_rhs_dual_recovery() {
        let a = cycle();
        let mut m = model(&a, &[3., -3., 0.]);
        assert_eq!(
            m.network_equalities(Instant::now() + Duration::from_secs(1)),
            1
        );
        let original = Point {
            x: vec![3., 0., 0.],
            y: vec![2., 3., 5.],
            z: vec![0.; 3],
        };
        let mut restored = original.clone();
        m.postsolve.reduce_point(&mut restored);
        m.postsolve.recover(&mut restored, Recovery::Solution);
        assert_eq!(restored.x, original.x);
        for (j, _) in a[0].iter().enumerate() {
            let gradient = |y: &[f64]| {
                a.iter()
                    .zip(y)
                    .map(|(row, multiplier)| row[j] * multiplier)
                    .sum::<f64>()
            };
            assert_eq!(gradient(&restored.y), gradient(&original.y));
        }
        let dual_offset = |y: &[f64]| 3. * y[0] - 3. * y[1];
        assert_eq!(dual_offset(&restored.y), dual_offset(&original.y));
    }
    #[test]
    fn private_columns_wrong_signs_nonzero_rhs_and_expired_budget_are_rejected() {
        let mut cases = vec![(cycle(), vec![0., 0., 1.])];
        let mut a = cycle();
        a[2][1] = 1.;
        cases.push((a, vec![0.; 3]));
        let mut a = cycle();
        a[2][2] = 0.;
        cases.push((a, vec![0.; 3]));
        cases.push((vec![vec![1.], vec![1.], vec![-2.]], vec![0.; 3]));
        for (a, b) in cases {
            let mut m = model(&a, &b);
            let before = m.a.nnz();
            assert_eq!(
                m.network_equalities(Instant::now() + Duration::from_secs(1)),
                0
            );
            assert_eq!(m.a.nnz(), before);
        }
        let mut m = model(&cycle(), &[0.; 3]);
        assert_eq!(m.network_equalities(Instant::now()), 0);
    }
    #[test]
    fn components_are_independent_and_inequalities_do_not_enter_the_relation() {
        let c = cycle();
        let mut a = vec![vec![0.; 6]; 7];
        for block in 0..2 {
            for i in 0..3 {
                for j in 0..3 {
                    a[block * 3 + i][block * 3 + j] = c[i][j];
                }
            }
        }
        a[6] = vec![1.; 6];
        let mut m = model(&a, &[0.; 7]);
        m.rows[6] = RowDomain::Linear(Bounds {
            lower: f64::NEG_INFINITY,
            upper: 10.,
        });
        assert_eq!(
            m.network_equalities(Instant::now() + Duration::from_secs(1)),
            2
        );
        assert_eq!(m.a.row(6).len(), 6);
    }
}
