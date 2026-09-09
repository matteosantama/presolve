//! Sparse presolve for convex quadratic objectives, ranged linear constraints,
//! variable bounds, and general convex cone blocks.
//!
//! Minimize `0.5 xᵀ P x + cᵀ x + c0`, subject to `row_bounds.lower ≤ A x ≤
//! row_bounds.upper`, `variable_bounds.lower ≤ x ≤ variable_bounds.upper`,
//! and `G x + s = h, s ∈ K`. Supply the upper triangle of the positive semidefinite
//! Hessian. All input validation, including dimensions, sparse structure,
//! finite values, convexity, and cone geometry, is the caller's responsibility.
//! Linear and conic rows share one sparse matrix, with explicit row domains.
//! Known cone blocks may be simplified; postsolve restores original coordinates.
//!
//! Native multipliers satisfy `P x + c - Aᵀ y - z + Gᵀ w = 0`.
//! Linear row and variable-bound multipliers are positive on lower sides and
//! negative on upper sides; conic multipliers use the usual dual-cone sign.
//! Postsolve maps coordinates without clipping or interiority adjustments.

//!
//! ```
//! use ::presolve::{presolve, Outcome, Settings};
//! use ::presolve::matrix::CscMatrix;
//! use ::presolve::problem::{Bounds, ProblemData};
//! let problem = ProblemData {
//!     p: None, c: vec![1.0], objective_constant: 0.0,
//!     a: CscMatrix::zeros(0, 1)?, rows: vec![],
//!     variable_bounds: vec![Bounds { lower: 2.0, upper: 4.0 }], cones: vec![],
//! };
//! let result = presolve(problem, &Settings::default());
//! if let Outcome::Solved(solution) = result.outcome {
//!     assert_eq!(solution.x, [2.0]);
//!     assert_eq!(solution.z, [1.0]);
//! }
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

mod core;
pub mod matrix;
pub mod postsolve;
mod presolve;
pub mod problem;
pub mod result;
pub mod settings;

pub use presolve::presolve;
pub use problem::Problem;
pub use result::Outcome;
pub use settings::Settings;
