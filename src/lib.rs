//! Sparse presolve for convex quadratic objectives, ranged linear constraints,
//! variable bounds, and general convex cone blocks.
//!
//! Minimize `0.5 xᵀ P x + cᵀ x + c0` subject to variable bounds `l ≤ x ≤ u`
//! and one constraint per row `aᵢ` of a shared sparse matrix `A`, whose kind
//! is given by the row's tag in `Problem::rows`:
//!
//! - `Constraint::Linear(bounds)`: a ranged row `bounds.lower ≤ aᵢ x ≤ bounds.upper`;
//! - `Constraint::Cone { rhs, block }`: one coordinate of a cone block,
//!   `aᵢ x + sᵢ = rhs`, where the slack vector `s` of the block's consecutive
//!   rows lies in `Problem::cones[block]`.
//!
//! Supply the upper triangle of the positive semidefinite Hessian `P`. All
//! input validation, including dimensions, sparse structure, finite values,
//! convexity, and cone geometry, is the caller's responsibility. Known cone
//! blocks may be simplified; postsolve restores original coordinates.
//!
//! Writing `A` for the linear rows and `G` for the conic rows, native
//! multipliers satisfy `P x + c - Aᵀ y - z + Gᵀ w = 0`.
//! Linear row and variable-bound multipliers are positive on lower sides and
//! negative on upper sides; conic multipliers use the usual dual-cone sign.
//! Postsolve maps coordinates without clipping or interiority adjustments.

//!
//! ```
//! use ::presolve::{Outcome, Presolver, Problem};
//! use ::presolve::matrix::CscMatrix;
//! use ::presolve::problem::Bounds;
//! let problem = Problem {
//!     p: None, c: vec![1.0], c0: 0.0,
//!     a: CscMatrix::zeros(0, 1)?.into(), rows: vec![],
//!     variable_bounds: vec![Bounds { lower: 2.0, upper: 4.0 }], cones: vec![],
//! };
//! let result = Presolver::default().presolve(problem);
//! if let Outcome::Solved(solution) = result.outcome {
//!     assert_eq!(solution.x, [2.0]);
//!     assert_eq!(solution.z, [1.0]);
//! }
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

mod executor;
pub mod matrix;
mod model;
pub mod postsolve;
mod presolve;
pub mod problem;
pub mod result;
mod rules;
pub mod settings;

pub use presolve::{InitError, Presolver};
pub use problem::Problem;
pub use result::Outcome;
pub use settings::Settings;
