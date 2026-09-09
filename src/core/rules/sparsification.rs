//! Bounded row cancellation on the native ranged model. Row/column storage,
//! bound proofs, cleanup queues, and postsolve are shared with other rules.
use crate::{
    core::model::{Model, RowDomain, shifted},
    matrix::sparse::Entries,
    postsolve::tape::Rule,
    problem::Bounds,
};
use std::time::Instant;

fn packed_sides(bounds: Bounds) -> usize {
    if bounds.equality() {
        1
    } else {
        usize::from(bounds.lower.is_finite()) + usize::from(bounds.upper.is_finite())
    }
}

fn worth_cancelling(length: usize, count: usize) -> bool {
    let required = if length < 20 { 10 } else { length / 2 };
    count.saturating_mul(2) >= length.saturating_add(required)
}

fn subtract<I: Iterator<Item = (usize, f64)> + Clone>(
    base: I,
    other: I,
    alpha: f64,
) -> Option<Entries> {
    let mut result = Vec::new();
    let mut max_new = 0.0_f64;
    let max_old = other.clone().fold(0.0_f64, |m, (_, a)| m.max(a.abs()));
    let mut base = base.peekable();
    let mut other = other.peekable();
    while base.peek().is_some() || other.peek().is_some() {
        let jb = base.peek().map_or(usize::MAX, |e| e.0);
        let jo = other.peek().map_or(usize::MAX, |e| e.0);
        let column = jb.min(jo);
        let b = if jb == column {
            base.next().unwrap().1
        } else {
            0.0
        };
        let a = if jo == column {
            other.next().unwrap().1
        } else {
            0.0
        };
        let v = a - alpha * b;
        if !v.is_finite() || (v != 0.0 && v.abs() < 1e-10 * a.abs().max((alpha * b).abs())) {
            return None;
        }
        if v != 0.0 {
            result.push((column, v));
            max_new = max_new.max(v.abs());
        }
    }
    (max_new <= 10.0 * max_old).then_some(result)
}

impl Model {
    pub fn sparsify_rows(&mut self, deadline: Instant) -> usize {
        let m = self.rows.len();
        let mut seen = vec![usize::MAX; m];
        let mut ratios = vec![0.0; m];
        let mut counts = vec![0usize; m];
        let mut candidates = Vec::new();
        // Bounds were explicit rows in the standalone pass. Include their
        // eventual packed entries when preserving its linear work allowance.
        let bound_entries: usize = self
            .bounds
            .iter()
            .enumerate()
            .filter(|(j, _)| self.alive[*j])
            .map(|(_, b)| packed_sides(*b))
            .sum();
        let limit = self
            .a
            .nnz()
            .saturating_add(bound_entries)
            .saturating_mul(8)
            .max(1024);
        let mut work = 0usize;
        let mut groups = 0;
        for reference in 0..m {
            if work >= limit || (reference % 64 == 0 && Instant::now() >= deadline) {
                break;
            }
            let RowDomain::Linear(bounds) = self.rows[reference] else {
                continue;
            };
            if self.a.row(reference).len() < 10 || packed_sides(bounds) == 0 {
                continue;
            }
            let base = self.a.row(reference);
            let auxiliary = (!bounds.equality()).then_some(self.bounds.len());
            let length = base.len() + usize::from(auxiliary.is_some());
            let minimum = base
                .iter()
                .map(|(_, a)| a.abs())
                .fold(f64::INFINITY, f64::min);
            if minimum < 1e-13 {
                continue;
            }
            let minimum_ratio = (1e-10 / minimum).max(1e-6);
            candidates.clear();
            'scan: for (j, value) in base {
                work += 1;
                for (i, other) in self.a.column(j) {
                    if work >= limit {
                        break 'scan;
                    }
                    work += 1;
                    if i == reference || i >= m || !matches!(self.rows[i], RowDomain::Linear(_)) {
                        continue;
                    }
                    if seen[i] != reference {
                        seen[i] = reference;
                        ratios[i] = 0.0;
                        counts[i] = 0;
                        candidates.push(i);
                    }
                    let ratio = other / value;
                    if !ratio.is_finite() {
                        continue;
                    }
                    if counts[i] == 0 {
                        ratios[i] = if ratio.abs() > minimum_ratio && ratio.abs() < 1e4 {
                            ratio
                        } else {
                            0.0
                        };
                        counts[i] = 1;
                    } else if ratio.abs() <= minimum_ratio {
                        ratios[i] = 0.0;
                        counts[i] = 0;
                    } else if (ratio - ratios[i]).abs() < 1e-10 {
                        counts[i] += 1;
                    } else if ratio.abs() < ratios[i].abs() {
                        ratios[i] = ratio;
                        counts[i] = 1;
                    }
                }
            }
            let mut updates = Vec::new();
            let mut targets = Vec::new();
            let mut saving = 0isize;
            for &i in &candidates {
                if !worth_cancelling(length, counts[i]) {
                    continue;
                }
                let cost = 3usize.saturating_mul(base.len() + self.a.row(i).len() + 1);
                if work.saturating_add(cost) > limit {
                    continue;
                }
                work += cost;
                let alpha = ratios[i];
                let Some(mut row) = subtract(base.iter(), self.a.row(i).iter(), alpha) else {
                    continue;
                };
                let RowDomain::Linear(target_bounds) = self.rows[i] else {
                    unreachable!();
                };
                let domain = if let Some(j) = auxiliary {
                    let max_old = self
                        .a
                        .row(i)
                        .iter()
                        .fold(0.0_f64, |v, (_, a)| v.max(a.abs()));
                    if alpha.abs() > 10.0 * max_old {
                        continue;
                    }
                    row.push((j, alpha));
                    self.rows[i]
                } else {
                    let Some(domain) = shifted(self.rows[i], alpha * bounds.lower) else {
                        continue;
                    };
                    domain
                };
                if row.len() >= self.a.row(i).len() {
                    continue;
                }
                saving +=
                    (packed_sides(target_bounds) * (self.a.row(i).len() - row.len())) as isize;
                updates.push((i, row, domain));
                targets.push((i, alpha));
            }
            if targets.is_empty() {
                continue;
            }
            if auxiliary.is_some() {
                saving += (packed_sides(bounds) * base.len()) as isize
                    - (base.len() + 1 + packed_sides(bounds)) as isize;
            }
            if saving <= 0 {
                continue;
            }
            let activity = if let Some(j) = auxiliary {
                let terms = base.to_vec();
                let mut equation = terms.clone();
                equation.push((j, -1.0));
                assert_eq!(self.add_variable(bounds), j);
                updates.push((reference, equation, RowDomain::Linear(Bounds::fixed(0.0))));
                Some((j, terms))
            } else {
                None
            };
            self.postsolve.rules.push(Rule::RowCombination {
                reference,
                targets,
                activity,
            });
            self.replace_rows_batch(updates);
            groups += 1;
        }
        groups
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn near_cancellations_and_overflow_are_not_rewritten() {
        let base: Entries = (0..20).map(|j| (j, 1.0)).collect();
        let mut other = base.clone();
        other[19].1 += 1e-12;
        assert!(subtract(base.iter().copied(), other.iter().copied(), 1.0).is_none());
        assert!(shifted(RowDomain::Linear(Bounds::fixed(-f64::MAX)), f64::MAX).is_none());
    }
}
