//! Exact LP folding through a verified equitable partition of the coefficient graph.
//! Equal coefficient multisets are a conservative alternative to rounded sum
//! signatures: they prove both block-sum identities without numerical tolerances.
use crate::{
    matrix::sparse::Entries,
    model::{Model, RowDomain, tape::Rule},
};
use std::{collections::HashMap, hash::Hash, time::Instant};

fn bits(x: f64) -> u64 {
    if x == 0.0 { 0 } else { x.to_bits() }
}

fn partition<K: Eq + Hash>(items: impl Iterator<Item = (usize, K)>, length: usize) -> Vec<usize> {
    let mut colors = vec![usize::MAX; length];
    let mut table = HashMap::new();
    for (i, key) in items {
        let next = table.len();
        colors[i] = *table.entry(key).or_insert(next);
    }
    colors
}

fn groups(colors: &[usize]) -> Vec<Vec<usize>> {
    let n = colors
        .iter()
        .filter(|&&c| c != usize::MAX)
        .max()
        .map_or(0, |c| c + 1);
    let mut groups = vec![Vec::new(); n];
    for (i, &c) in colors.iter().enumerate() {
        if c != usize::MAX {
            groups[c].push(i);
        }
    }
    groups
}

/// Accumulate an error-free expansion, then round its total on export.
/// The partition proof is exact; compressed coefficients have ordinary f64 rounding.
fn finite_sum(mut values: impl Iterator<Item = f64>) -> Option<f64> {
    let Some(first) = values.next() else {
        return Some(0.0);
    };
    let Some(second) = values.next() else {
        return Some(first);
    };
    let mut expansion = vec![first];
    for mut x in std::iter::once(second).chain(values) {
        let mut write = 0;
        for read in 0..expansion.len() {
            let y = expansion[read];
            let s = x + y;
            let b = s - x;
            let error = (x - (s - b)) + (y - b);
            if !s.is_finite() || !error.is_finite() {
                return None;
            }
            if error != 0.0 {
                expansion[write] = error;
                write += 1;
            }
            x = s;
        }
        expansion.truncate(write);
        if x != 0.0 {
            expansion.push(x);
        }
    }
    let total: f64 = expansion.into_iter().sum();
    total.is_finite().then_some(total)
}

fn spend(work: &mut usize, amount: usize) -> bool {
    if amount > *work {
        return false;
    }
    *work -= amount;
    true
}

impl Model {
    fn fold_refine<const ROW: bool>(
        &self,
        own: &[usize],
        opposite: &[usize],
        work: &mut usize,
        deadline: Instant,
    ) -> Option<Vec<usize>> {
        let mut result = vec![usize::MAX; own.len()];
        let mut signatures = HashMap::new();
        for (i, &color) in own.iter().enumerate() {
            if i % 256 == 0 && Instant::now() >= deadline {
                return None;
            }
            if color == usize::MAX {
                continue;
            }
            let length = if ROW {
                self.a.row(i).len()
            } else {
                self.a.column(i).len()
            };
            let cost = length
                .saturating_mul(length.max(1).ilog2() as usize + 1)
                .saturating_add(1);
            if !spend(work, cost) {
                return None;
            }
            let mut signature = Vec::with_capacity(length);
            if ROW {
                signature.extend(self.a.row(i).iter().map(|(j, a)| (opposite[j], bits(a))));
            } else {
                signature.extend(self.a.column(i).iter().map(|(j, a)| (opposite[j], bits(a))));
            }
            signature.sort_unstable();
            let next = signatures.len();
            result[i] = *signatures.entry((color, signature)).or_insert(next);
        }
        Some(result)
    }

    /// Restrict an LP to class-constant primal coordinates and average rows.
    /// The transpose partition identity proves that averaging any feasible
    /// point preserves feasibility and objective; no optimal point is lost.
    pub fn fold_lp(&mut self, deadline: Instant) -> usize {
        if !self.cones.is_empty()
            || self.objective.p.nnz() != 0
            || self.settings.folding.max_rounds == 0
        {
            return 0;
        }
        let mut work = self.settings.folding.work_limit.resolve(
            self.a
                .nnz()
                .saturating_mul(64)
                .saturating_add(self.rows.len())
                .saturating_add(self.bounds.len()),
        );
        if !spend(&mut work, self.rows.len().saturating_add(self.bounds.len()))
            || Instant::now() >= deadline
        {
            return 0;
        }
        let mut rows = partition(
            self.rows.iter().enumerate().filter_map(|(i, r)| match r {
                RowDomain::Linear(b) => Some((i, (bits(b.lower), bits(b.upper)))),
                _ => None,
            }),
            self.rows.len(),
        );
        let mut columns = partition(
            self.bounds.iter().enumerate().filter_map(|(j, b)| {
                self.alive[j]
                    .then_some((j, (bits(self.objective.c[j]), bits(b.lower), bits(b.upper))))
            }),
            self.bounds.len(),
        );
        let mut stable = false;
        for _ in 0..self.settings.folding.max_rounds {
            let Some(r) = self.fold_refine::<true>(&rows, &columns, &mut work, deadline) else {
                return 0;
            };
            let Some(c) = self.fold_refine::<false>(&columns, &r, &mut work, deadline) else {
                return 0;
            };
            stable = r == rows && c == columns;
            rows = r;
            columns = c;
            if stable {
                break;
            }
        }
        // A partial refinement is not necessarily equitable and cannot be used.
        if !stable {
            return 0;
        }
        let row_groups = groups(&rows);
        let column_groups = groups(&columns);
        let removed = row_groups
            .iter()
            .chain(&column_groups)
            .map(|g| g.len() - 1)
            .sum();
        if removed == 0 {
            return 0;
        }
        let mut costs = Vec::with_capacity(column_groups.len());
        for group in &column_groups {
            // The exact partition proves identical objective coefficients.
            // Scale once (one ordinary floating-point rounding), as in column
            // aggregation; demanding exact representability rejects e.g. 3*c.
            let cost = self.objective.c[group[0]] * group.len() as f64;
            if !cost.is_finite() {
                return 0;
            }
            costs.push(cost);
        }
        let mut updates = Vec::new();
        for group in &row_groups {
            if Instant::now() >= deadline {
                return 0;
            }
            let i = group[0];
            if !spend(
                &mut work,
                self.a
                    .row(i)
                    .len()
                    .saturating_mul(2)
                    .saturating_add(group.len()),
            ) {
                return 0;
            }
            let mut entries: Vec<_> = self.a.row(i).iter().map(|(j, a)| (columns[j], a)).collect();
            entries.sort_unstable_by_key(|&(c, _)| c);
            let mut merged = Entries::new();
            let mut k = 0;
            while k < entries.len() {
                let color = entries[k].0;
                let end = k + entries[k..].partition_point(|&(c, _)| c == color);
                let Some(value) = finite_sum(entries[k..end].iter().map(|&(_, a)| a)) else {
                    return 0;
                };
                if value != 0.0 {
                    merged.push((column_groups[color][0], value));
                }
                k = end;
            }
            merged.sort_unstable_by_key(|&(j, _)| j);
            updates.push((i, merged, self.rows[i]));
            updates.extend(
                group[1..]
                    .iter()
                    .map(|&i| (i, Entries::new(), RowDomain::Deleted)),
            );
        }
        if Instant::now() >= deadline {
            return 0;
        }
        self.replace_rows_batch(updates);
        for (group, cost) in column_groups.iter().zip(costs) {
            self.objective.c[group[0]] = cost;
            for &j in &group[1..] {
                debug_assert!(self.a.column(j).is_empty());
                self.objective.c[j] = 0.0;
                // Even an empty supplied Hessian has the original dimensions.
                // Invalidate reuse when folding changes the variable indexing.
                self.objective.p.remove_variable(j);
                self.alive[j] = false;
                self.revision += 1;
            }
        }
        self.postsolve.rules.push(Rule::LpFold {
            columns: column_groups.into_iter().filter(|g| g.len() > 1).collect(),
            rows: row_groups.into_iter().filter(|g| g.len() > 1).collect(),
        });
        removed
    }
}

#[cfg(test)]
mod tests {
    use super::finite_sum;
    #[test]
    fn compensated_sums_preserve_small_terms_and_reject_overflow() {
        assert_eq!(finite_sum([0.1, 0.1, 0.1, 0.1].into_iter()), Some(0.4));
        assert_eq!(finite_sum([0.1, 0.2].into_iter()), Some(0.1 + 0.2));
        assert_eq!(finite_sum([1e16, 1.0, -1e16].into_iter()), Some(1.0));
        assert_eq!(finite_sum([f64::MAX, f64::MAX].into_iter()), None);
    }
}
