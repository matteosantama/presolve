//! Replace an implied one-sided variable bound by a doubleton inequality.
//! The invertible coordinate change keeps both variables, but removes the row.
use crate::result::RuleId;
use crate::{
    model::{Model, RowDomain, shifted, tape::Rule},
    problem::Bounds,
};
use std::time::Instant;

// Reject multiplication that silently loses a nonzero matrix/objective term.
fn product(a: f64, b: f64) -> Option<f64> {
    let value = a * b;
    (value.is_finite() && (value != 0.0 || a == 0.0 || b == 0.0)).then_some(value)
}

impl Model {
    pub fn bound_shift(&mut self, deadline: Instant) -> usize {
        self.enter(RuleId::BoundShift);
        let options = self.settings.bound_shift;
        if !options.max_ratio.is_finite() || options.max_ratio < 1.0 {
            return 0;
        }
        let mut remaining = options.work_limit.resolve(self.a.nnz().saturating_mul(64));
        let mut applied = 0;
        for row in 0..self.rows.len() {
            if remaining == 0 || Instant::now() >= deadline {
                break;
            }
            remaining -= 1;
            let RowDomain::Linear(domain) = self.rows[row] else {
                continue;
            };
            if domain.lower.is_finite() == domain.upper.is_finite() || self.a.row(row).len() != 2 {
                continue;
            }
            let entries = self.a.row(row).to_vec();
            let mut pivots = [0, 1];
            pivots.sort_unstable_by_key(|&index| self.a.column(entries[index].0).len());
            for index in pivots {
                let (column, pivot) = entries[index];
                let (other, coefficient) = entries[1 - index];
                if !self.objective.p.column(column).is_empty()
                    || self.a.column(column).len() > options.max_column_length
                {
                    continue;
                }
                let lower = (pivot > 0.0) == domain.lower.is_finite();
                let old = self.bounds[column];
                if (lower && old.upper.is_finite()) || (!lower && old.lower.is_finite()) {
                    continue;
                }
                let rhs = if domain.lower.is_finite() {
                    domain.lower
                } else {
                    domain.upper
                };
                let offset = rhs / pivot;
                let slope = -coefficient / pivot;
                if !offset.is_finite()
                    || !slope.is_finite()
                    || slope.abs() > options.max_ratio
                    || slope.abs() < 1.0 / options.max_ratio
                    || slope * pivot != -coefficient
                    || offset * pivot != rhs
                {
                    continue;
                }
                let bounds = self.bounds[other];
                let endpoint = if lower == (slope > 0.0) {
                    bounds.lower
                } else {
                    bounds.upper
                };
                let term = if endpoint.is_finite() {
                    let Some(term) = product(slope, endpoint) else {
                        continue;
                    };
                    term
                } else {
                    slope * endpoint
                };
                let extreme = term + offset;
                if (lower && extreme < old.lower)
                    || (!lower && extreme > old.upper)
                    || extreme.is_nan()
                {
                    continue;
                }
                let cost = self.a.column(column).iter().fold(0usize, |sum, (i, _)| {
                    sum.saturating_add(self.a.row(i).len().saturating_add(1))
                });
                if cost > remaining {
                    continue;
                }
                remaining -= cost;
                let Some(objective_shift) = product(slope, self.objective.c[column]) else {
                    continue;
                };
                let Some(constant_shift) = product(offset, self.objective.c[column]) else {
                    continue;
                };
                let objective = self.objective.c[other] + objective_shift;
                let constant = self.objective.constant + constant_shift;
                if !objective.is_finite() || !constant.is_finite() {
                    continue;
                }
                let mut updates = Vec::with_capacity(self.a.column(column).len());
                let mut fill = 0usize;
                let mut valid = true;
                for (i, value) in self.a.column(column) {
                    if Instant::now() >= deadline {
                        valid = false;
                        break;
                    }
                    if i == row {
                        continue;
                    }
                    if !matches!(self.rows[i], RowDomain::Linear(_)) {
                        valid = false;
                        break;
                    }
                    let Some(shift) = product(offset, value) else {
                        valid = false;
                        break;
                    };
                    let Some(domain) = shifted(self.rows[i], shift) else {
                        valid = false;
                        break;
                    };
                    let mut values = self.a.row(i).to_vec();
                    let position = values.binary_search_by_key(&other, |&(j, _)| j);
                    let previous = position.map_or(0.0, |k| values[k].1);
                    let Some(addition) = product(slope, value) else {
                        valid = false;
                        break;
                    };
                    let next = previous + addition;
                    if !next.is_finite()
                        || (next != 0.0 && next.abs() < 1e-12 * previous.abs().max(addition.abs()))
                    {
                        valid = false;
                        break;
                    }
                    match position {
                        Ok(k) if next == 0.0 => {
                            values.remove(k);
                        }
                        Ok(k) => values[k].1 = next,
                        Err(k) if next != 0.0 => {
                            values.insert(k, (other, next));
                            fill += 1;
                        }
                        Err(_) => {}
                    }
                    if fill > options.max_fill {
                        valid = false;
                        break;
                    }
                    updates.push((i, values, domain));
                }
                if !valid || Instant::now() >= deadline {
                    continue;
                }
                self.record(Rule::BoundShift {
                    column,
                    other,
                    row,
                    pivot,
                    slope,
                    offset,
                });
                self.objective.c[other] = objective;
                self.objective.constant = constant;
                self.replace_rows_batch(updates);
                self.clear_row(row);
                self.set_bounds(
                    column,
                    if lower {
                        Bounds {
                            lower: 0.0,
                            upper: f64::INFINITY,
                        }
                    } else {
                        Bounds {
                            lower: f64::NEG_INFINITY,
                            upper: 0.0,
                        }
                    },
                );
                applied += 1;
                break;
            }
        }
        applied
    }
}
