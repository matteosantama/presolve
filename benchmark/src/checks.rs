//! Relative violations of optimality and certificate conditions for native
//! points of a native problem. Every measure is zero for an exact answer.
//!
//! With `A` the linear rows, `G` the conic rows, `u` the row multipliers
//! (`y` on linear rows, `-w` on conic rows) and `s = h - Gx` the conic slack:
//!
//! - `primal_rows`, `primal_bounds`: the largest bound violation of a linear
//!   row or variable, over `max(1, |bound|, Σⱼ|aᵢⱼxⱼ|)` (rows) or
//!   `max(1, |bound|)` (variables);
//! - `primal_cones`: the largest distance of a block of `s` from its cone,
//!   over `max(1, ‖s‖∞)` for zero, nonnegative and second-order cones;
//! - `stationarity`: the largest entry of `Px + c - Aᵀy - z + Gᵀw`, each over
//!   the largest of 1, `|cⱼ|`, `|(Px)ⱼ|`, `Σᵢ|aᵢⱼuᵢ|` and `|zⱼ|`;
//! - `dual_signs`: the largest multiplier on an infinite side, or distance of
//!   a block of `w` from its dual cone, over `max(1, ‖(y, z, w)‖∞)`;
//! - `dual_feasibility`: the largest bound multiplier `Px + c - Aᵀy + Gᵀw`,
//!   which absorbs the stationarity residual, on an infinite side, over the
//!   scale of `stationarity` without `|zⱼ|`;
//! - `complementarity`: `Σ|multiplier|·|distance to its side| + Σ|sᵀw|` over
//!   `max(1, |objective|)`;
//! - `duality_gap`: the objective minus the dual objective
//!   `c0 - ½xᵀPx + Σ multiplier·side - hᵀw` over finite sides, over
//!   `max(1, |objective|)`, with bound multipliers `Px + c - Aᵀy + Gᵀw` that
//!   absorb the stationarity residual; this bounds the objective error of a
//!   feasible `x` only when `dual_feasibility` is zero, and a tiny residual
//!   on a huge finite bound makes it large, so it only calibrates;
//! - `certificate_signs`: as `dual_signs`, for a Farkas certificate;
//! - `certificate_value`: `-v / Σ|terms of v|`, where `v`, the selected
//!   side-weighted sum minus `hᵀw`, is positive for a contradiction;
//! - `certificate_residual`: `‖Aᵀy + z - Gᵀw‖∞ / v`, infinite when `v ≤ 0`;
//! - `ray_objective`: `cᵀd / Σ|cⱼdⱼ|`, negative for a descent direction;
//! - `ray_hessian`, `ray_recession`: `‖Pd‖∞` and the largest violation of the
//!   recession cone of the rows, bounds and cones, over `-cᵀd` with
//!   `‖d‖∞ = 1`, infinite when `cᵀd ≥ 0`.

use presolve::Problem;
use presolve::postsolve::{PrimalCertificate, Solution};
use presolve::problem::{Bounds, Cone, Constraint};
use std::collections::BTreeMap;

pub type Measures = BTreeMap<&'static str, f64>;

/// A problem with its rows split into linear rows and cone blocks.
pub struct Model<'a> {
    problem: &'a Problem,
    /// Row index of each linear row, and of each conic row, in order.
    linear: Vec<usize>,
    conic: Vec<usize>,
    /// Each cone with the conic positions of its block.
    blocks: Vec<(Cone, Vec<usize>)>,
}

fn max_abs(v: &[f64]) -> f64 {
    v.iter().fold(0.0, |m, x| m.max(x.abs()))
}

/// Violation of `bounds` by `value`, with the magnitude of the violated side.
fn violation(bounds: Bounds, value: f64) -> (f64, f64) {
    if value < bounds.lower {
        (bounds.lower - value, bounds.lower.abs())
    } else if value > bounds.upper {
        (value - bounds.upper, bounds.upper.abs())
    } else {
        (0.0, 0.0)
    }
}

/// The side a multiplier selects, `None` when it is zero, and whether it
/// lies on an infinite side.
fn side(bounds: Bounds, m: f64) -> Option<(f64, bool)> {
    let side = if m > 0.0 {
        bounds.lower
    } else if m < 0.0 {
        bounds.upper
    } else {
        return None;
    };
    Some((side, side.is_infinite()))
}

/// Distance of `v` from a cone, or from its dual when `dual`. `None` when
/// the cone is not checked.
fn cone_distance(cone: Cone, v: &[f64], dual: bool) -> Option<f64> {
    match cone {
        Cone::Zero(_) if dual => Some(0.0),
        Cone::Zero(_) => Some(max_abs(v)),
        Cone::Nonnegative(_) => Some(v.iter().fold(0.0, |m, &x| m.max(-x))),
        Cone::SecondOrder(_) if v.is_empty() => Some(0.0),
        Cone::SecondOrder(_) => {
            let norm = v[1..].iter().fold(0.0_f64, |n, &x| n.hypot(x));
            Some((norm - v[0]).max(0.0))
        }
        _ => None,
    }
}

impl<'a> Model<'a> {
    pub fn new(problem: &'a Problem) -> Self {
        let mut linear = Vec::new();
        let mut conic = Vec::new();
        let mut blocks: Vec<(Cone, Vec<usize>)> =
            problem.cones.iter().map(|&c| (c, Vec::new())).collect();
        for (i, row) in problem.rows.iter().enumerate() {
            match *row {
                Constraint::Linear(_) => linear.push(i),
                Constraint::Cone { block, .. } => {
                    blocks[block].1.push(conic.len());
                    conic.push(i);
                }
            }
        }
        Self {
            problem,
            linear,
            conic,
            blocks,
        }
    }

    fn row_bounds(&self, k: usize) -> Bounds {
        match self.problem.rows[self.linear[k]] {
            Constraint::Linear(b) => b,
            Constraint::Cone { .. } => unreachable!(),
        }
    }

    fn rhs(&self, k: usize) -> f64 {
        match self.problem.rows[self.conic[k]] {
            Constraint::Cone { rhs, .. } => rhs,
            Constraint::Linear(_) => unreachable!(),
        }
    }

    /// `½xᵀPx + cᵀx + c0`.
    pub fn objective(&self, x: &[f64]) -> f64 {
        let px = self.hessian(x);
        let quadratic: f64 = px.iter().zip(x).map(|(a, b)| a * b).sum();
        let linear: f64 = self.problem.c.iter().zip(x).map(|(a, b)| a * b).sum();
        0.5 * quadratic + linear + self.problem.c0
    }

    /// `Ax` and `Σⱼ|aᵢⱼxⱼ|` for every row.
    fn activity(&self, x: &[f64]) -> (Vec<f64>, Vec<f64>) {
        let a = self.problem.a.as_ref();
        let mut value = vec![0.0; a.rows()];
        let mut size = vec![0.0; a.rows()];
        for (j, &xj) in x.iter().enumerate() {
            for (i, v) in a.column(j) {
                value[i] += v * xj;
                size[i] += (v * xj).abs();
            }
        }
        (value, size)
    }

    /// `Px` from the upper triangle.
    fn hessian(&self, x: &[f64]) -> Vec<f64> {
        let mut out = vec![0.0; x.len()];
        if let Some(p) = &self.problem.p {
            let p = p.as_ref();
            for (j, &xj) in x.iter().enumerate() {
                for (i, v) in p.column(j) {
                    out[i] += v * xj;
                    if i != j {
                        out[j] += v * x[i];
                    }
                }
            }
        }
        out
    }

    /// Row multipliers: `y` on linear rows and `-w` on conic rows.
    fn row_multipliers(&self, y: &[f64], w: &[f64]) -> Vec<f64> {
        let mut u = vec![0.0; self.problem.row_count()];
        for (&i, &v) in self.linear.iter().zip(y) {
            u[i] = v;
        }
        for (&i, &v) in self.conic.iter().zip(w) {
            u[i] = -v;
        }
        u
    }

    /// `Aᵀu` and `Σᵢ|aᵢⱼuᵢ|` for every column.
    fn transpose(&self, u: &[f64]) -> (Vec<f64>, Vec<f64>) {
        let a = self.problem.a.as_ref();
        (0..a.columns())
            .map(|j| {
                let (mut sum, mut abs) = (0.0, 0.0);
                for (i, v) in a.column(j) {
                    sum += v * u[i];
                    abs += (v * u[i]).abs();
                }
                (sum, abs)
            })
            .unzip()
    }

    fn block(&self, positions: &[usize], values: &[f64], scratch: &mut Vec<f64>) {
        scratch.clear();
        scratch.extend(positions.iter().map(|&k| values[k]));
    }

    /// Row and bound violations of `x`, and the conic slack `h - Gx`.
    fn primal_into(&self, x: &[f64], out: &mut Measures) -> Vec<f64> {
        let (activity, size) = self.activity(x);
        let rows = (0..self.linear.len())
            .map(|k| {
                let i = self.linear[k];
                let (v, side) = violation(self.row_bounds(k), activity[i]);
                v / side.max(size[i]).max(1.0)
            })
            .fold(0.0, f64::max);
        let bounds = x
            .iter()
            .enumerate()
            .map(|(j, &xj)| {
                let (v, side) = violation(self.problem.variable_bounds(j), xj);
                v / side.max(1.0)
            })
            .fold(0.0, f64::max);
        let slack: Vec<f64> = (0..self.conic.len())
            .map(|k| self.rhs(k) - activity[self.conic[k]])
            .collect();
        out.insert("primal_rows", rows);
        out.insert("primal_bounds", bounds);
        if !self.conic.is_empty() {
            let mut scratch = Vec::new();
            let mut worst = 0.0_f64;
            for (cone, positions) in &self.blocks {
                self.block(positions, &slack, &mut scratch);
                if let Some(d) = cone_distance(*cone, &scratch, false) {
                    worst = worst.max(d / max_abs(&scratch).max(1.0));
                }
            }
            out.insert("primal_cones", worst);
        }
        slack
    }

    /// Feasibility of `x` alone.
    pub fn primal(&self, x: &[f64]) -> Measures {
        let mut out = Measures::new();
        self.primal_into(x, &mut out);
        out
    }

    /// Largest multiplier on an infinite side or outside the dual cone.
    fn signs(&self, y: &[f64], z: &[f64], w: &[f64]) -> f64 {
        let scale = max_abs(y).max(max_abs(z)).max(max_abs(w)).max(1.0);
        self.wrong_signs(y, z, w) / scale
    }

    fn wrong_signs(&self, y: &[f64], z: &[f64], w: &[f64]) -> f64 {
        let mut worst = 0.0_f64;
        for (k, &m) in y.iter().enumerate() {
            if let Some((_, true)) = side(self.row_bounds(k), m) {
                worst = worst.max(m.abs());
            }
        }
        for (j, &m) in z.iter().enumerate() {
            if let Some((_, true)) = side(self.problem.variable_bounds(j), m) {
                worst = worst.max(m.abs());
            }
        }
        let mut scratch = Vec::new();
        for (cone, positions) in &self.blocks {
            self.block(positions, w, &mut scratch);
            worst = worst.max(cone_distance(*cone, &scratch, true).unwrap_or(0.0));
        }
        worst
    }

    /// All optimality measures of a primal-dual point.
    pub fn optimality(&self, s: &Solution) -> Measures {
        let mut out = Measures::new();
        let slack = self.primal_into(&s.x, &mut out);
        let px = self.hessian(&s.x);
        let u = self.row_multipliers(&s.y, &s.conic_dual);
        let (atu, atu_size) = self.transpose(&u);
        // Bound multipliers that absorb the stationarity residual.
        let absorbed: Vec<f64> = (0..s.x.len())
            .map(|j| px[j] + self.problem.c[j] - atu[j])
            .collect();
        let (mut residual, mut infeasibility) = (0.0_f64, 0.0_f64);
        for j in 0..s.x.len() {
            let scale = self.problem.c[j]
                .abs()
                .max(px[j].abs())
                .max(atu_size[j])
                .max(1.0);
            residual = residual.max((absorbed[j] - s.z[j]).abs() / scale.max(s.z[j].abs()));
            if let Some((_, true)) = side(self.problem.variable_bounds(j), absorbed[j]) {
                infeasibility = infeasibility.max(absorbed[j].abs() / scale);
            }
        }
        out.insert("stationarity", residual);
        out.insert("dual_signs", self.signs(&s.y, &s.z, &s.conic_dual));
        out.insert("dual_feasibility", infeasibility);
        let (activity, _) = self.activity(&s.x);
        let mut gap = 0.0;
        for (k, &m) in s.y.iter().enumerate() {
            if let Some((side, false)) = side(self.row_bounds(k), m) {
                gap += (m * (activity[self.linear[k]] - side)).abs();
            }
        }
        for (j, &m) in s.z.iter().enumerate() {
            if let Some((side, false)) = side(self.problem.variable_bounds(j), m) {
                gap += (m * (s.x[j] - side)).abs();
            }
        }
        let (mut a, mut b) = (Vec::new(), Vec::new());
        for (_, positions) in &self.blocks {
            self.block(positions, &slack, &mut a);
            self.block(positions, &s.conic_dual, &mut b);
            gap += a.iter().zip(&b).map(|(p, q)| p * q).sum::<f64>().abs();
        }
        let objective = self.objective(&s.x);
        out.insert("complementarity", gap / objective.abs().max(1.0));
        let quadratic: f64 = px.iter().zip(&s.x).map(|(a, b)| a * b).sum();
        let (sides, _) = self.side_sum(&s.y, &absorbed, &s.conic_dual);
        let dual = self.problem.c0 - 0.5 * quadratic + sides;
        out.insert(
            "duality_gap",
            (objective - dual).abs() / objective.abs().max(1.0),
        );
        out
    }

    /// The sum over finite sides of multiplier times side, minus `hᵀw`,
    /// and the sum of the absolute terms.
    fn side_sum(&self, y: &[f64], z: &[f64], w: &[f64]) -> (f64, f64) {
        let (mut value, mut size) = (0.0, 0.0);
        let mut add = |term: f64| {
            value += term;
            size += term.abs();
        };
        for (k, &m) in y.iter().enumerate() {
            if let Some((side, false)) = side(self.row_bounds(k), m) {
                add(m * side);
            }
        }
        for (j, &m) in z.iter().enumerate() {
            if let Some((side, false)) = side(self.problem.variable_bounds(j), m) {
                add(m * side);
            }
        }
        for (k, &w) in w.iter().enumerate() {
            add(-self.rhs(k) * w);
        }
        (value, size)
    }

    /// Measures of a Farkas certificate.
    pub fn farkas(&self, c: &PrimalCertificate) -> Measures {
        let mut out = Measures::new();
        out.insert("certificate_signs", self.signs(&c.y, &c.z, &c.conic_dual));
        let (value, size) = self.side_sum(&c.y, &c.z, &c.conic_dual);
        out.insert(
            "certificate_value",
            if size > 0.0 { -value / size } else { 0.0 },
        );
        let u = self.row_multipliers(&c.y, &c.conic_dual);
        let (atu, _) = self.transpose(&u);
        let residual = atu
            .iter()
            .zip(&c.z)
            .fold(0.0_f64, |m, (a, z)| m.max((a + z).abs()));
        out.insert(
            "certificate_residual",
            if value > 0.0 {
                residual / value
            } else {
                f64::INFINITY
            },
        );
        out
    }

    /// Measures of a recession direction.
    pub fn ray(&self, d: &[f64]) -> Measures {
        let mut out = Measures::new();
        let norm = max_abs(d);
        let d: Vec<f64> = d.iter().map(|v| v / norm.max(f64::MIN_POSITIVE)).collect();
        let (descent, size) = self
            .problem
            .c
            .iter()
            .zip(&d)
            .fold((0.0, 0.0), |(s, a), (c, d)| (s + c * d, a + (c * d).abs()));
        out.insert(
            "ray_objective",
            if size > 0.0 { descent / size } else { 0.0 },
        );
        let per = |v: f64| {
            if descent < 0.0 {
                v / -descent
            } else {
                f64::INFINITY
            }
        };
        out.insert("ray_hessian", per(max_abs(&self.hessian(&d))));
        let (activity, _) = self.activity(&d);
        let mut worst = 0.0_f64;
        for k in 0..self.linear.len() {
            let (v, _) = violation(self.row_bounds(k).recession(), activity[self.linear[k]]);
            worst = worst.max(v);
        }
        for (j, &dj) in d.iter().enumerate() {
            let (v, _) = violation(self.problem.variable_bounds(j).recession(), dj);
            worst = worst.max(v);
        }
        let direction: Vec<f64> = self.conic.iter().map(|&i| -activity[i]).collect();
        let mut scratch = Vec::new();
        for (cone, positions) in &self.blocks {
            self.block(positions, &direction, &mut scratch);
            worst = worst.max(cone_distance(*cone, &scratch, false).unwrap_or(0.0));
        }
        out.insert("ray_recession", per(worst));
        out
    }
}
