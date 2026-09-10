//! Bounded exact-arithmetic equality dependencies, without changing A or P
//! until a complete relation has been established in scratch storage.
use crate::{
    core::model::{Model, RowDomain},
    matrix::sparse::Entries,
    postsolve::tape::{Certificate, Point, Recovery, Rule},
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
    pub fn equality_dependencies(&mut self, deadline: Instant) -> Result<usize, Certificate> {
        let options = self.dependencies;
        let mut work = options.work_limit.resolve(self.a.nnz().saturating_mul(4));
        if work == 0 || options.max_basis_rows == 0 || options.max_row_length == 0 {
            return Ok(0);
        }
        let mut basis: Vec<BasisRow> = Vec::new();
        let mut row = Vec::new();
        let mut proof = Vec::new();
        let mut next_row = Vec::new();
        let mut next_proof = Vec::new();
        let mut removed = 0;
        'rows: for i in 0..self.rows.len() {
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
                    self.replace_row(i, &[], RowDomain::Deleted);
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
                    if rhs.abs() > self.numerics.feasibility * (1.0 + scale) {
                        let mut point = Point::zeros(self.bounds.len(), self.rows.len());
                        for &(j, a) in &proof {
                            point.y[j] = rhs.signum() * a;
                        }
                        return Err(Certificate {
                            mode: Recovery::PrimalInfeasibility,
                            point,
                        });
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
