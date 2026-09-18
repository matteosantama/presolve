//! Bounded exact-arithmetic equality dependencies, without changing A or P
//! until a complete relation has been established in scratch storage.
use crate::{
    matrix::sparse::Entries,
    model::tape::{Certificate, Rule},
    model::{Model, RowDomain},
};
use std::time::Instant;

struct BasisRow {
    pivot: usize,
    value: f64,
    entries: Entries,
    rhs: f64,
    // This scratch row is a linear combination of original model rows.
    proof: Entries,
}

fn exact_product(a: f64, b: f64) -> Option<f64> {
    let product = a * b;
    if a == 0.0 || b == 0.0 {
        return Some(0.0);
    }
    // Below this conservative threshold the product's rounding residual could
    // underflow, defeating the FMA exactness check. Skip such transformations.
    (product.is_finite()
        && product.abs() >= f64::MIN_POSITIVE * 18_014_398_509_481_984.0
        && a.mul_add(b, -product) == 0.0)
        .then_some(product)
}

fn exact_difference(a: f64, b: f64) -> Option<f64> {
    let value = a - b;
    let bv = a - value;
    let av = value + bv;
    let error = (a - av) + (bv - b);
    (value.is_finite() && bv.is_finite() && av.is_finite() && error == 0.0).then_some(value)
}

fn subtract(target: &Entries, base: &Entries, alpha: f64, out: &mut Entries) -> Option<()> {
    out.clear();
    let (mut i, mut j) = (0, 0);
    while i < target.len() || j < base.len() {
        let ti = target.get(i).map_or(usize::MAX, |e| e.0);
        let bj = base.get(j).map_or(usize::MAX, |e| e.0);
        let column = ti.min(bj);
        let a = if ti == column {
            let a = target[i].1;
            i += 1;
            a
        } else {
            0.0
        };
        let b = if bj == column {
            let b = base[j].1;
            j += 1;
            b
        } else {
            0.0
        };
        let value = exact_difference(a, exact_product(alpha, b)?)?;
        if value != 0.0 {
            out.push((column, value));
        }
    }
    Some(())
}

impl Model {
    /// A row with a private column cannot occur in a linear dependence. Peeling
    /// is confined to the same eligible equality subsystem as the exact search;
    /// inequalities and oversized equalities must not contribute to degrees.
    fn dependency_core(&self, deadline: Instant, work: &mut usize) -> Vec<usize> {
        let mut rows = Vec::new();
        let mut entries = 0usize;
        for (i, domain) in self.rows.iter().enumerate() {
            if i % 1024 == 0 && Instant::now() >= deadline {
                return Vec::new();
            }
            if !matches!(domain, RowDomain::Linear(b) if b.equality()) {
                continue;
            }
            let length = self.a.row(i).len();
            if length != 0 && length <= self.settings.dependencies.max_row_length {
                rows.push(i);
                entries = entries.saturating_add(length);
                // The unpeeled search cannot visit beyond this prefix under
                // the allowance. Do not scan the rest of a large model just
                // to discover that preprocessing cannot fit either.
                if entries > *work {
                    return rows;
                }
            }
        }
        // For a tiny allowance, spend it on the original search instead of
        // initializing column scratch. Account for both initialization and
        // incidence visits, leaving at least half the allowance after counting.
        let count_cost = entries.saturating_add(self.bounds.len());
        if rows.len() < 2 || count_cost > *work / 2 {
            return rows;
        }
        // If the subsystem contains every matrix entry, the arena's cached
        // singleton counts already prove whether any peeling can start. Dense
        // equality cores avoid column scratch and a redundant incidence scan.
        if entries == self.a.nnz() && !rows.iter().any(|&i| self.a.row_singletons(i) != 0) {
            return rows;
        }
        *work -= count_cost;
        let mut degree = vec![0usize; self.bounds.len()];
        let mut owner = vec![0usize; self.bounds.len()];
        for (position, &i) in rows.iter().enumerate() {
            if position % 256 == 0 && Instant::now() >= deadline {
                return Vec::new();
            }
            for (j, _) in self.a.row(i) {
                degree[j] += 1;
                owner[j] ^= position;
            }
        }
        let mut leaves: Vec<_> = degree
            .iter()
            .enumerate()
            .filter_map(|(j, &d)| (d == 1).then_some(j))
            .collect();
        while let Some(j) = leaves.pop() {
            if degree[j] != 1 {
                continue;
            }
            let position = owner[j];
            let i = rows[position];
            let length = self.a.row(i).len();
            if length > *work || Instant::now() >= deadline {
                break;
            }
            *work -= length;
            // XOR identifies the sole remaining row without searching a linked
            // column. Decrementing every incidence also invalidates queued leaves
            // from this same row, so no separate removed-row flags are needed.
            rows[position] = usize::MAX;
            for (k, _) in self.a.row(i) {
                degree[k] -= 1;
                owner[k] ^= position;
                if degree[k] == 1 {
                    leaves.push(k);
                }
            }
        }
        rows.retain(|&i| i != usize::MAX);
        rows
    }

    pub fn equality_dependencies(&mut self, deadline: Instant) -> Result<usize, Certificate> {
        let options = self.settings.dependencies;
        let mut work = options.work_limit.resolve(self.a.nnz().saturating_mul(4));
        if work == 0 || options.max_basis_rows == 0 || options.max_row_length == 0 {
            return Ok(0);
        }
        let candidates = self.dependency_core(deadline, &mut work);
        let mut basis: Vec<BasisRow> = Vec::new();
        let mut row = Vec::new();
        let mut proof = Vec::new();
        let mut next_row = Vec::new();
        let mut next_proof = Vec::new();
        let mut removed = 0;
        'rows: for i in candidates {
            if work == 0 || Instant::now() >= deadline {
                break;
            }
            let RowDomain::Linear(bounds) = self.rows[i] else {
                continue;
            };
            let length = self.a.row(i).len();
            if !bounds.equality() || length == 0 || length > options.max_row_length {
                continue;
            }
            if length > work {
                break;
            }
            work -= length;
            row.clear();
            row.extend(self.a.row(i));
            proof.clear();
            proof.push((i, 1.0));
            let mut rhs = bounds.lower;
            for reference in &basis {
                if work == 0 || Instant::now() >= deadline {
                    break 'rows;
                }
                work = work.saturating_sub(row.len().ilog2() as usize + 1);
                let Ok(at) = row.binary_search_by_key(&reference.pivot, |e| e.0) else {
                    continue;
                };
                let alpha = row[at].1 / reference.value;
                if !alpha.is_finite() || alpha == 0.0 {
                    continue 'rows;
                }
                let cost = row
                    .len()
                    .saturating_add(reference.entries.len())
                    .saturating_add(proof.len())
                    .saturating_add(reference.proof.len());
                if cost > work {
                    break 'rows;
                }
                work -= cost;
                let Some(new_rhs) =
                    exact_product(alpha, reference.rhs).and_then(|v| exact_difference(rhs, v))
                else {
                    continue 'rows;
                };
                if subtract(&row, &reference.entries, alpha, &mut next_row).is_none()
                    || next_row.len() > options.max_row_length
                    || subtract(&proof, &reference.proof, alpha, &mut next_proof).is_none()
                {
                    continue 'rows;
                }
                // Division must actually cancel the chosen pivot exactly.
                if next_row
                    .binary_search_by_key(&reference.pivot, |e| e.0)
                    .is_ok()
                {
                    continue 'rows;
                }
                std::mem::swap(&mut row, &mut next_row);
                std::mem::swap(&mut proof, &mut next_proof);
                rhs = new_rhs;
                if row.is_empty() {
                    break;
                }
            }
            if row.is_empty() {
                if rhs == 0.0 {
                    // Forward warm starts transfer the deleted row's multiplier
                    // to its surviving proof rows. Reverse recovery sets it to zero.
                    let coefficients = proof
                        .iter()
                        .filter(|e| e.0 != i)
                        .map(|&(j, a)| (j, -a))
                        .collect();
                    self.clear_row(i);
                    self.postsolve.rules.push(Rule::DependentRow {
                        row: i,
                        coefficients,
                    });
                    removed += 1;
                } else {
                    let scale: f64 = proof
                        .iter()
                        .map(|&(j, a)| {
                            let RowDomain::Linear(b) = self.rows[j] else {
                                unreachable!()
                            };
                            (a * b.lower).abs()
                        })
                        .sum();
                    if rhs.abs() > self.settings.numerics.feasibility * (1.0 + scale) {
                        return Err(self.primal_certificate(
                            proof.iter().map(|&(j, a)| (j, rhs.signum() * a)),
                            [],
                        ));
                    }
                }
            } else if basis.len() < options.max_basis_rows {
                let &(pivot, value) = row
                    .iter()
                    .max_by(|a, b| a.1.abs().total_cmp(&b.1.abs()).then_with(|| b.0.cmp(&a.0)))
                    .unwrap();
                basis.push(BasisRow {
                    pivot,
                    value,
                    entries: std::mem::take(&mut row),
                    rhs,
                    proof: std::mem::take(&mut proof),
                });
            }
        }
        Ok(removed)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn arithmetic_checks_reject_rounding_and_underflow() {
        assert_eq!(exact_product(0.5, 3.), Some(1.5));
        assert_eq!(exact_product(0.1, 3.), None);
        assert_eq!(exact_product(1e-300, 1e-100), None);
        assert_eq!(exact_product(f64::MAX, 2.), None);
        assert_eq!(exact_difference(2., 1.), Some(1.));
        assert_eq!(exact_difference(1., 1e-20), None);
    }
}
