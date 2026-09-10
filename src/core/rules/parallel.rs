// SPDX-License-Identifier: Apache-2.0
// Modified for this library; copyright and attribution notices are in NOTICE.
//! Detect parallel rows and columns. Support and coefficient hashes
//! propose groups; sparse comparisons establish proportionality before edits.

use crate::{
    core::{
        execution::Executor,
        model::{Model, RowDomain},
    },
    postsolve::tape::{Certificate, Point, Recovery, Rule, Side},
    problem::Bounds,
};

// Avoid scheduling small scans. Work estimates count constraint entries and,
// for columns, full symmetric Hessian entries. These are initial cutoffs, not
// guarantees of a speedup on every matrix or machine.
const MIN_PARALLEL_ITEMS: usize = 1024;
const MIN_PARALLEL_NONZEROS: usize = 32 * 1024;

/// djb2 support hash and rounded, sign-normalized coefficient hash.
/// Division by the maximum avoids overflow from forming its reciprocal first.
fn fingerprint<I: Iterator<Item = (usize, f64)> + Clone>(entries: I) -> (u32, u32) {
    let max = entries.clone().map(|(_, v)| v.abs()).fold(0.0, f64::max);
    let sign = entries.clone().next().unwrap().1.signum();
    entries.fold((5381u32, 5381u32), |(support, coefficients), (j, a)| {
        let quantized = (1e6 * (a / max) * sign).round() as i32 as u32;
        (
            support.wrapping_mul(33).wrapping_add(j as u32),
            coefficients.wrapping_mul(33).wrapping_add(quantized),
        )
    })
}

/// Returns A_base / A_other, plus whether equality is exact in floating-point
/// arithmetic. The latter is required before two inequalities become equality.
fn proportional(
    base: crate::matrix::linked::View<'_>,
    other: crate::matrix::linked::View<'_>,
    tolerance: f64,
) -> Option<(f64, bool)> {
    if base.is_empty() || base.len() != other.len() {
        return None;
    }
    let ratio = base.iter().next().unwrap().1 / other.iter().next().unwrap().1;
    if !ratio.is_finite() || ratio == 0.0 {
        return None;
    }
    let mut exact = true;
    for ((j, a), (k, b)) in base.iter().zip(other) {
        let scaled = ratio * b;
        if j != k
            || !scaled.is_finite()
            || (a - scaled).abs() > tolerance * a.abs().max(scaled.abs())
        {
            return None;
        }
        exact &= a == scaled;
    }
    Some((ratio, exact))
}

fn candidate_groups<K: Ord + Send>(
    count: usize,
    nonzeros: usize,
    executor: &Executor,
    key: impl Fn(usize) -> Option<K> + Sync,
) -> Vec<Vec<usize>> {
    let enough_work = count >= MIN_PARALLEL_ITEMS && nonzeros >= MIN_PARALLEL_NONZEROS;
    let entry = |i| key(i).map(|key| (key, i));
    // The index is a unique tiebreaker, making group and member order identical
    // regardless of worker scheduling or sorting algorithm.
    let entries = executor.filter_map_sorted(count, enough_work, entry);
    let mut out = Vec::new();
    let mut start = 0;
    while start < entries.len() {
        let mut end = start + 1;
        while end < entries.len() && entries[end].0 == entries[start].0 {
            end += 1;
        }
        if end - start > 1 {
            out.push(entries[start..end].iter().map(|(_, i)| *i).collect());
        }
        start = end;
    }
    out
}

impl Model {
    pub fn parallel_rows(&mut self, executor: &Executor) -> Result<usize, Certificate> {
        let groups = candidate_groups(self.rows.len(), self.a.nnz(), executor, |i| {
            (matches!(self.rows[i], RowDomain::Linear(_)) && self.a.row(i).len() > 1)
                .then(|| fingerprint(self.a.row(i).iter()))
        });
        let mut comparisons = 0;
        let mut unmatched = Vec::new();
        for group in groups {
            let base = group[0];
            unmatched.clear();
            // Preserve the usual linear-time path and its representative order.
            for &other in &group[1..] {
                comparisons += 1;
                if !self.merge_parallel_rows(base, other)? {
                    unmatched.push((other, 0.0));
                }
            }
            if unmatched.len() < 2 {
                continue;
            }
            // A coarse fingerprint can contain several proportional classes.
            // Canonical ordering brings their members together without an
            // all-pairs search. Normalize once per row; avoid reciprocals that
            // overflow for small coefficients. Only collision bins pay this cost.
            for (i, scale) in &mut unmatched {
                let row = self.a.row(*i);
                *scale = row
                    .iter()
                    .map(|(_, a)| a.abs())
                    .fold(0.0, f64::max)
                    .copysign(row.iter().next().unwrap().1);
            }
            unmatched.sort_unstable_by(|&(i, a), &(j, b)| {
                let left = self.a.row(i);
                let right = self.a.row(j);
                left.len()
                    .cmp(&right.len())
                    .then_with(|| {
                        left.iter()
                            .zip(right)
                            .find_map(|((k, x), (l, y))| {
                                let order = k.cmp(&l).then_with(|| (x / a).total_cmp(&(y / b)));
                                (!order.is_eq()).then_some(order)
                            })
                            .unwrap_or(std::cmp::Ordering::Equal)
                    })
                    .then_with(|| i.cmp(&j))
            });
            let mut base = unmatched[0].0;
            for &(other, _) in &unmatched[1..] {
                comparisons += 1;
                if !self.merge_parallel_rows(base, other)? {
                    base = other;
                }
            }
        }
        Ok(comparisons)
    }

    /// Authoritative comparison and mutation, shared by both discovery paths.
    fn merge_parallel_rows(&mut self, base: usize, other: usize) -> Result<bool, Certificate> {
        let Some((ratio, exact)) =
            proportional(self.a.row(base), self.a.row(other), self.numerics.parallel)
        else {
            return Ok(false);
        };
        let RowDomain::Linear(b) = self.rows[base] else {
            return Ok(false);
        };
        let RowDomain::Linear(c) = self.rows[other] else {
            return Ok(false);
        };
        let (cl, cu) = if ratio > 0.0 {
            (c.lower, c.upper)
        } else {
            (c.upper, c.lower)
        };
        let lower = ratio * cl;
        let upper = ratio * cu;
        if (cl.is_finite() && !lower.is_finite()) || (cu.is_finite() && !upper.is_finite()) {
            return Ok(false);
        }
        let intersection = Bounds {
            lower: b.lower.max(lower),
            upper: b.upper.min(upper),
        };
        if intersection.lower > intersection.upper {
            let gap = intersection.lower - intersection.upper;
            if !exact
                || gap
                    <= self.numerics.feasibility
                        * (1.0 + intersection.lower.abs().max(intersection.upper.abs()))
            {
                return Ok(false);
            }
            let mut point = Point::zeros(self.bounds.len(), self.rows.len());
            point.y[base] = if b.lower > upper { 1.0 } else { -1.0 };
            point.y[other] = -ratio * point.y[base];
            return Err(Certificate {
                mode: Recovery::PrimalInfeasibility,
                point,
            });
        }
        if intersection.equality() && !b.equality() && !exact {
            return Ok(false);
        }
        if intersection.lower > b.lower {
            self.tighten_row(base, other, ratio, Side::Lower, intersection.lower);
        }
        if intersection.upper < b.upper {
            self.tighten_row(base, other, ratio, Side::Upper, intersection.upper);
        }
        self.replace_row(other, &[], RowDomain::Deleted);
        self.postsolve.rules.push(Rule::MergedRow {
            keep: base,
            removed: other,
            ratio,
        });
        Ok(true)
    }

    pub fn parallel_columns(&mut self, executor: &Executor) -> Result<usize, Certificate> {
        let groups = candidate_groups(
            self.bounds.len(),
            self.a.nnz().saturating_add(self.objective.p.nnz()),
            executor,
            |j| {
                (self.alive[j] && !self.a.column(j).is_empty()).then(|| {
                    let p = self.objective.p.column(j);
                    (
                        fingerprint(self.a.column(j).iter()),
                        (!p.is_empty()).then(|| fingerprint(p.iter().copied())),
                    )
                })
            },
        );
        let mut comparisons = 0;
        for group in groups {
            for (at, &j) in group.iter().enumerate() {
                if !self.alive[j] {
                    continue;
                }
                for &k in &group[at + 1..] {
                    if !self.alive[k] {
                        continue;
                    }
                    comparisons += 1;
                    let Some((ratio, exact)) =
                        proportional(self.a.column(k), self.a.column(j), self.numerics.parallel)
                    else {
                        continue;
                    };
                    // An approximate null direction can have quadratic cost
                    // at large x. Require exact A and P relations for aggregation
                    // and objective-dominance certificates in the QP port.
                    if !exact || !self.objective.parallel_curvature(j, k, ratio, 0.0) {
                        continue;
                    }
                    if self.aggregate(j, k, ratio, 0.0) {
                        continue;
                    }
                    let delta = self.objective.c[k] - ratio * self.objective.c[j];
                    if !delta.is_finite() || delta == 0.0 {
                        continue;
                    }
                    let dk = -delta.signum();
                    let dj = -ratio * dk;
                    let bj = if dj > 0.0 {
                        self.bounds[j].upper
                    } else {
                        self.bounds[j].lower
                    };
                    let bk = if dk > 0.0 {
                        self.bounds[k].upper
                    } else {
                        self.bounds[k].lower
                    };
                    if !bj.is_finite() && !bk.is_finite() {
                        let mut point = Point::zeros(self.bounds.len(), self.rows.len());
                        point.x[j] = dj;
                        point.x[k] = dk;
                        return Err(Certificate {
                            mode: Recovery::DualInfeasibility,
                            point,
                        });
                    }
                    if !bj.is_finite() {
                        self.fix(k, bk);
                    } else if !bk.is_finite() && self.fix(j, bj) {
                        break;
                    }
                }
            }
        }
        Ok(comparisons)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn discovery_uses_requested_threads_and_preserves_filtered_group_order() {
        let key = |i: usize| (i % 5 != 0).then_some(i % 17);
        let executor = Executor::new(3).unwrap();
        let serial = candidate_groups(2048, MIN_PARALLEL_NONZEROS, &Executor::Serial, key);
        let parallel = candidate_groups(2048, MIN_PARALLEL_NONZEROS, &executor, |i| {
            assert_eq!(rayon::current_num_threads(), 3);
            key(i)
        });
        assert_eq!(serial, parallel);

        // Both the slot and work thresholds must be met before using a pool.
        for (count, work) in [(1023, MIN_PARALLEL_NONZEROS), (2048, 32767)] {
            candidate_groups(count, work, &executor, |i| {
                assert!(rayon::current_thread_index().is_none());
                key(i)
            });
        }
    }
}
