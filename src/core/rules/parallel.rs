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
use wide::{f64x4, u32x4};

// Avoid scheduling small scans. Work estimates count constraint entries and,
// for columns, full symmetric Hessian entries. These are initial cutoffs, not
// guarantees of a speedup on every matrix or machine.
const MIN_PARALLEL_ITEMS: usize = 1024;
const MIN_PARALLEL_NONZEROS: usize = 32 * 1024;

/// djb2 support hash and rounded, sign-normalized coefficient hash.
/// Division by the maximum avoids overflow from forming its reciprocal first.
#[inline]
fn fingerprint<I: ExactSizeIterator<Item = (usize, f64)> + Clone>(mut entries: I) -> (u32, u32) {
    let max = entries.clone().map(|(_, v)| v.abs()).fold(0.0, f64::max);
    let sign = entries.clone().next().unwrap().1.signum();
    let scalar = |(support, coefficients): (u32, u32), (j, a): (usize, f64)| {
        let quantized = (1e6 * (a / max) * sign).round() as i32 as u32;
        (
            support.wrapping_mul(33).wrapping_add(j as u32),
            coefficients.wrapping_mul(33).wrapping_add(quantized),
        )
    };
    // Packing costs more than it saves on very short rows/columns.
    if entries.len() < 8 {
        return entries.fold((5381, 5381), scalar);
    }

    // Four interleaved djb2 streams advance by 33^4 per block. Weighting
    // them by [33^3, 33^2, 33, 1] at the end reconstructs the original hash
    // modulo 2^32, including the seed. This avoids a reduction in every block.
    const SEED: u32x4 = u32x4::new([0, 0, 0, 5381]);
    const STRIDE: u32x4 = u32x4::new([1_185_921; 4]);
    const WEIGHTS: u32x4 = u32x4::new([35_937, 1089, 33, 1]);
    const MILLION: f64x4 = f64x4::new([1e6; 4]);
    let mut support = SEED;
    let mut coefficients = support;
    let scale = f64x4::splat(max);
    let sign = f64x4::splat(sign);
    while entries.len() >= 4 {
        // Stack packing works for both contiguous Hessian adjacency and
        // linked matrix iterators without allocating or changing traversal.
        let block = std::array::from_fn::<_, 4, _>(|_| entries.next().unwrap());
        let values = f64x4::new(block.map(|(_, a)| a));
        let normalized = MILLION * (values / scale) * sign;
        // Keep Rust's ties-away-from-zero rounding in each lane.
        // Normalization bounds finite values by 1e6, so converting
        // via i64 preserves the original i32 result (and NaN still maps to 0)
        // while allowing SIMD conversion on targets with 64-bit float lanes.
        let quantized = normalized.to_array().map(|a| a.round() as i64 as u32);
        support = support * STRIDE + u32x4::new(block.map(|(j, _)| j as u32));
        coefficients = coefficients * STRIDE + u32x4::new(quantized);
    }
    let collapse = |lanes: u32x4| {
        (lanes * WEIGHTS)
            .to_array()
            .into_iter()
            .fold(0u32, u32::wrapping_add)
    };
    entries.fold((collapse(support), collapse(coefficients)), scalar)
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

    fn scalar_fingerprint<I: Iterator<Item = (usize, f64)> + Clone>(entries: I) -> (u32, u32) {
        let max = entries.clone().map(|(_, a)| a.abs()).fold(0.0, f64::max);
        let sign = entries.clone().next().unwrap().1.signum();
        entries.fold((5381u32, 5381u32), |(h, g), (j, a)| {
            let q = (1e6 * (a / max) * sign).round() as i32 as u32;
            (
                h.wrapping_mul(33).wrapping_add(j as u32),
                g.wrapping_mul(33).wrapping_add(q),
            )
        })
    }

    #[test]
    fn fingerprints_preserve_scalar_rounding_tails_and_overflow() {
        let mut random = 0x1234_5678_abcd_ef01u64;
        for len in 1..=129 {
            for _ in 0..16 {
                let entries: Vec<_> = (0..len)
                    .map(|_| {
                        random ^= random << 13;
                        random ^= random >> 7;
                        random ^= random << 17;
                        // Include subnormals and both signs, excluding NaN/Inf.
                        let bits = (random & 0x800f_ffff_ffff_ffff) | ((random % 2047) << 52);
                        (random as usize, f64::from_bits(bits))
                    })
                    .collect();
                assert_eq!(
                    fingerprint(entries.iter().copied()),
                    scalar_fingerprint(entries.iter().copied()),
                    "length {len}"
                );
            }
        }

        let half = 0.5f64;
        let cases = [
            vec![
                1e6,
                half.next_down(),
                half,
                half.next_up(),
                -half.next_down(),
                -half,
                -half.next_up(),
                1.5,
                -1.5,
                2.5,
                -2.5,
            ],
            vec![
                f64::MAX,
                -f64::MAX,
                f64::MIN_POSITIVE,
                f64::from_bits(1),
                -f64::from_bits(1),
                0.0,
                -0.0,
                1.0,
                -1.0,
            ],
            vec![
                f64::from_bits(1),
                -f64::from_bits(1),
                f64::from_bits(2),
                -f64::from_bits(2),
                0.0,
                -0.0,
                f64::from_bits(3),
                f64::from_bits(4),
            ],
            vec![0.0; 11],
            vec![
                f64::INFINITY,
                f64::NAN,
                1.0,
                -1.0,
                -f64::INFINITY,
                0.0,
                -0.0,
                2.0,
                3.0,
            ],
        ];
        for values in cases {
            for sign in [1.0, -1.0] {
                let entries: Vec<_> = values
                    .iter()
                    .enumerate()
                    .map(|(j, a)| (usize::MAX - j, a * sign))
                    .collect();
                assert_eq!(
                    fingerprint(entries.iter().copied()),
                    scalar_fingerprint(entries.iter().copied())
                );
            }
        }
    }

    #[test]
    fn fingerprints_match_on_fragmented_linked_rows_and_columns() {
        use crate::matrix::linked::LinkedMatrix;
        for len in [2, 7, 8, 9, 15, 16, 17, 64, 129] {
            for (rows, cols) in [(len, 3), (3, len)] {
                let mut matrix = LinkedMatrix::from_columns(rows, cols, |j| {
                    (0..rows).map(move |i| {
                        (
                            i,
                            (1 + i + 7 * j) as f64 * if i % 2 == 0 { 1.0 } else { -1.0 },
                        )
                    })
                });
                // Reusing the free list in deletion order reverses physical
                // slot order, while the logical rows and columns stay sorted.
                for i in (0..rows).step_by(2) {
                    for j in (0..cols).step_by(2) {
                        matrix.set(i, j, 0.0);
                    }
                }
                for i in (0..rows).step_by(2) {
                    for j in (0..cols).step_by(2) {
                        matrix.set(i, j, (1 + i + 7 * j) as f64);
                    }
                }
                for view in (0..rows)
                    .map(|i| matrix.row(i))
                    .chain((0..cols).map(|j| matrix.column(j)))
                {
                    let expected = scalar_fingerprint(view.iter());
                    assert_eq!(fingerprint(view.iter()), expected);
                    assert_eq!(fingerprint(view.to_vec().iter().copied()), expected);
                }
            }
        }
    }

    #[test]
    fn discovery_uses_requested_threads_and_preserves_filtered_group_order() {
        let key = |i: usize| (!i.is_multiple_of(5)).then_some(i % 17);
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
