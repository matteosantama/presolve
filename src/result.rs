//! Presolve outcomes, model sizes, and execution statistics.
use crate::postsolve::{Postsolve, PrimalCertificate, Solution};
use crate::problem::Problem;
use std::time::Duration;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Size {
    pub variables: usize,
    pub linear_rows: usize,
    pub conic_rows: usize,
    /// Nonzero values, excluding explicitly stored zeros.
    pub a_nonzeros: usize,
    pub g_nonzeros: usize,
    /// Full symmetric storage: each off-diagonal pair counts twice.
    pub p_nonzeros: usize,
}
#[derive(Clone, Debug)]
pub struct Stats {
    pub elapsed: Duration,
    pub time_limit_reached: bool,
    pub before: Size,
    /// None when presolve stops on a certificate before producing a problem.
    pub after: Option<Size>,
    pub parallel_comparisons: usize,
    /// Whether Hessian coefficients or variable coordinates changed. A caller
    /// retaining its original matrix can reuse it when this is false.
    pub quadratic_changed: bool,
}

#[derive(Debug)]
pub struct PresolveResult {
    pub outcome: Outcome,
    pub stats: Stats,
}
// Returning unchanged data inline avoids an allocation on that path.
#[allow(clippy::large_enum_variant)]
#[derive(Debug)]
pub enum Outcome {
    /// Return the input allocations when no rule was applied.
    Unchanged(Problem),
    Reduced(Box<ReducedProblem>),
    /// Original-coordinate solution. Never emitted with remaining opaque cone blocks:
    /// the caller must check their feasibility even if all variables are fixed.
    Solved(Solution),
    Infeasible(PrimalCertificate),
    /// A feasible point and an objective-decreasing recession direction.
    Unbounded(UnboundednessCertificate),
}
#[derive(Debug)]
pub struct ReducedProblem {
    pub problem: Problem,
    pub postsolve: Postsolve,
}

#[derive(Clone, Debug)]
pub struct UnboundednessCertificate {
    pub point: Vec<f64>,
    pub ray: Vec<f64>,
}
