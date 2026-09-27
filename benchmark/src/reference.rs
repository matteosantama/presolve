//! Clarabel as a reference solver for native problems.
//!
//! A problem is exported with `Problem::into_conic` and solved as
//! `min ½xᵀPx + cᵀx  s.t.  Ax + s = b, s ∈ K`. Clarabel's dual `z` satisfies
//! `Px + c + Aᵀz = 0`, which is exactly what `ConicMap::dual` expects, so the
//! map turns it into native multipliers, and a primal-infeasibility `z` into a
//! native Farkas certificate.

use clarabel::algebra::CscMatrix as ClarabelMatrix;
use clarabel::solver::{
    DefaultSettingsBuilder, DefaultSolver, IPSolver, SolverStatus, SupportedConeT,
};
use presolve::Problem;
use presolve::matrix::CscMatrix;
use presolve::postsolve::{PrimalCertificate, Solution};
use presolve::problem::Cone;
use std::panic::{AssertUnwindSafe, catch_unwind};

/// A reference solve in native coordinates.
#[derive(Debug)]
pub struct Reference {
    /// Clarabel's status, as its `Debug` name.
    pub status: String,
    pub result: Answer,
    pub iterations: u32,
    pub seconds: f64,
}

#[derive(Debug)]
pub enum Answer {
    /// Solved to full accuracy: a native primal-dual point.
    Optimal(Solution),
    /// A native Farkas certificate.
    Infeasible(PrimalCertificate),
    /// A recession direction `d` with `cᵀd < 0`.
    Unbounded(Vec<f64>),
    /// Reduced accuracy, a limit, a numerical failure, or unsupported data.
    Unknown,
}

/// Solve `problem` with Clarabel, with `tolerance` on gap, feasibility and
/// infeasibility (Clarabel's default is 1e-8), at most 400 iterations, and a
/// time limit in seconds. Errors are data Clarabel rejects or a panic inside
/// it.
pub fn solve(problem: &Problem, tolerance: f64, time_limit: f64) -> Result<Reference, String> {
    let export = problem.clone().into_conic();
    let (data, map) = (export.problem, export.map);
    let n = data.c.len();
    let mut cones = Vec::with_capacity(data.cones.len());
    for &cone in &data.cones {
        cones.push(match cone {
            _ if cone.dimension() == 0 => continue,
            Cone::Zero(k) => SupportedConeT::ZeroConeT(k),
            Cone::Nonnegative(k) => SupportedConeT::NonnegativeConeT(k),
            Cone::SecondOrder(k) => SupportedConeT::SecondOrderConeT(k),
            Cone::Exponential => SupportedConeT::ExponentialConeT(),
            Cone::Power { alpha } => SupportedConeT::PowerConeT(alpha),
            other => return Err(format!("unsupported cone {other:?}")),
        });
    }
    let p = data
        .p
        .as_ref()
        .map_or_else(|| ClarabelMatrix::zeros((n, n)), convert);
    let a = convert(&data.a);
    let settings = DefaultSettingsBuilder::default()
        .verbose(false)
        .time_limit(time_limit)
        .tol_gap_abs(tolerance)
        .tol_gap_rel(tolerance)
        .tol_feas(tolerance)
        .tol_ktratio(tolerance * 100.0)
        .tol_infeas_abs(tolerance)
        .tol_infeas_rel(tolerance)
        .max_iter(400)
        .build()
        .map_err(|e| e.to_string())?;
    let solved = catch_unwind(AssertUnwindSafe(|| {
        let mut solver = DefaultSolver::new(&p, &data.c, &a, &data.b, &cones, settings)
            .map_err(|e| e.to_string())?;
        solver.solve();
        Ok::<_, String>(solver.solution)
    }))
    .map_err(|_| "Clarabel panicked".to_string())??;
    let result = match solved.status {
        SolverStatus::Solved => {
            let native = map.dual(&solved.z);
            Answer::Optimal(Solution {
                x: solved.x,
                y: native.y,
                z: native.z,
                conic_dual: native.conic_dual,
                conic_slack: solved.s[map.conic_offset()..].to_vec(),
            })
        }
        SolverStatus::PrimalInfeasible => Answer::Infeasible(map.dual(&solved.z)),
        SolverStatus::DualInfeasible => Answer::Unbounded(solved.x),
        _ => Answer::Unknown,
    };
    Ok(Reference {
        status: format!("{:?}", solved.status),
        result,
        iterations: solved.iterations,
        seconds: solved.solve_time,
    })
}

fn convert(m: &CscMatrix) -> ClarabelMatrix<f64> {
    ClarabelMatrix::new(
        m.rows(),
        m.columns(),
        m.column_pointers().to_vec(),
        m.row_indices().to_vec(),
        m.values().to_vec(),
    )
}
