//! One-dimensional convex dominance with exact coefficient matching and
//! outward-rounded cost comparisons. No coefficient or objective is rewritten.
use crate::{model::Model, result::RuleId, settings::Progress};
use std::{
    hash::{DefaultHasher, Hasher},
    time::Instant,
};

#[derive(Clone, Copy)]
struct Point {
    a: f64,
    c: f64,
    j: usize,
}

fn difference(a: f64, b: f64) -> (f64, f64) {
    let d = a - b;
    (d.next_down(), d.next_up())
}
fn product(a: (f64, f64), b: (f64, f64)) -> (f64, f64) {
    let v = [a.0 * b.0, a.0 * b.1, a.1 * b.0, a.1 * b.1];
    if v.iter().any(|x| !x.is_finite()) {
        return (f64::NEG_INFINITY, f64::INFINITY);
    }
    (
        v.into_iter().fold(f64::INFINITY, f64::min).next_down(),
        v.into_iter().fold(f64::NEG_INFINITY, f64::max).next_up(),
    )
}
/// Certify that b lies strictly above the chord a--c in exact real arithmetic
/// on the stored f64 inputs. Uncertain, collinear, or overflowing cases stay.
fn above(a: Point, b: Point, c: Point) -> bool {
    let left = product(difference(b.a, a.a), difference(c.c, a.c));
    let right = product(difference(b.c, a.c), difference(c.a, a.a));
    left.1 < right.0
}

impl Model {
    /// Every eligible column has bounds [0,+inf] and no quadratic terms.
    /// A group must match exactly in every matrix entry except one row.
    /// Replacing an interior variable by a convex mixture of two survivors
    /// preserves Ax and cannot increase c'x. Fixed records recover its reduced
    /// cost as a nonnegative mixture of the survivors' costs plus the saving.
    pub fn convex_dominance(&mut self, deadline: Instant) -> usize {
        self.enter(RuleId::ConvexDominance);
        let mut candidates = Vec::new();
        for j in 0..self.bounds.len() {
            if j % 256 == 0 && Instant::now() >= deadline {
                return 0;
            }
            if !self.alive[j]
                || self.bounds[j].lower != 0.0
                || self.bounds[j].upper != f64::INFINITY
                || !self.linear_column(j)
                || self.a.column(j).is_empty()
            {
                continue;
            }
            let mut hash = DefaultHasher::new();
            for (i, _) in self.a.column(j) {
                hash.write_usize(i);
            }
            candidates.push((hash.finish(), j));
        }
        candidates.sort_unstable();
        let before = self.work_size();
        let mut removals = Vec::new();
        let mut removed_nonzeros = 0;
        let mut start = 0;
        let mut points = Vec::new();
        let mut hull: Vec<Point> = Vec::new();
        while start < candidates.len() {
            if Instant::now() >= deadline {
                return 0;
            }
            let mut end = start + 1;
            while end < candidates.len() && candidates[end].0 == candidates[start].0 {
                end += 1;
            }
            let group = &candidates[start..end];
            start = end;
            if group.len() < 3 {
                continue;
            }
            let reference: Vec<_> = self.a.column(group[0].1).iter().collect();
            let mut varying = None;
            let mut valid = true;
            for &(_, j) in group {
                if self.a.column(j).len() != reference.len() {
                    valid = false;
                    break;
                }
                for ((i, a), &(k, b)) in self.a.column(j).iter().zip(&reference) {
                    if i != k {
                        valid = false;
                        break;
                    }
                    if a != b {
                        if varying.is_some_and(|v| v != i) {
                            valid = false;
                            break;
                        }
                        varying = Some(i);
                    }
                }
                if !valid {
                    break;
                }
            }
            if !valid {
                continue;
            }
            let Some(varying) = varying else {
                continue;
            };
            points.clear();
            for &(_, j) in group {
                let a = self
                    .a
                    .column(j)
                    .iter()
                    .find(|&(i, _)| i == varying)
                    .unwrap()
                    .1;
                points.push(Point {
                    a,
                    c: self.objective.c[j],
                    j,
                });
            }
            points.sort_unstable_by(|a, b| {
                a.a.total_cmp(&b.a)
                    .then(a.c.total_cmp(&b.c))
                    .then(a.j.cmp(&b.j))
            });
            hull.clear();
            for &p in &points {
                if Instant::now() >= deadline {
                    return 0;
                }
                if hull.last().is_some_and(|last| last.a == p.a) {
                    removals.push(p.j);
                    removed_nonzeros += self.a.column(p.j).len();
                    continue;
                }
                while hull.len() >= 2 && above(hull[hull.len() - 2], hull[hull.len() - 1], p) {
                    let removed = hull.pop().unwrap();
                    removals.push(removed.j);
                    removed_nonzeros += self.a.column(removed.j).len();
                }
                hull.push(p);
            }
        }
        // Stage the whole pass before changing the model. Sparse local pruning
        // can alter the solver trajectory without reducing its matrix work
        // appreciably. Use the caller's convex-dominance reduction threshold;
        // the aggressive preset uses zero to accept every nonempty plan.
        // These columns have no Hessian entries, so the projected work reduction
        // is exactly their A entries, with no fill or coefficient changes.
        if Instant::now() >= deadline
            || !super::significant_progress(
                Progress::Nonzeros {
                    minimum_reduction: self.settings.convex_dominance.minimum_reduction,
                },
                before,
                before - removed_nonzeros,
                !removals.is_empty(),
            )
        {
            return 0;
        }
        // Commit the accepted batch as one soft-budget transaction. Partial
        // application would defeat the work-reduction gate. All fixes are at
        // the existing zero lower bound of a linear column.
        removals
            .into_iter()
            .map(|j| usize::from(self.fix(j, 0.0)))
            .sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn p(a: f64, c: f64) -> Point {
        Point { a, c, j: 0 }
    }
    #[test]
    fn conservative_orientation() {
        assert!(above(p(1., 1.), p(2., 3.), p(3., 3.)));
        assert!(!above(p(1., 1.), p(2., 2.), p(3., 3.)));
        assert!(!above(p(1., 1.), p(2., 1.), p(3., 3.)));
        assert!(!above(
            p(-f64::MAX, -f64::MAX),
            p(0., 1.),
            p(f64::MAX, f64::MAX)
        ));
    }
}
