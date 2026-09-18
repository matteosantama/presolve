//! Owned native problem data with explicit, consuming CSC export.
mod bounds;
mod cone;
mod conic;

pub use bounds::Bounds;
pub use cone::Cone;
pub(crate) use cone::Membership;
pub use conic::{ConicData, ConicExport, ConicMap};

use crate::{
    matrix::{CscMatrix, sparse::SymmetricMatrix},
    model::objective::Objective,
};

/// One row of the shared constraint matrix. Cone coordinates in a block are
/// consecutive; `block` indexes the cone list, independently of opaque IDs.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Constraint {
    Linear(Bounds),
    Cone { rhs: f64, block: usize },
}

/// Minimize `0.5 xᵀ P x + cᵀ x + c0` over `x` with `variable_bounds[j].lower ≤
/// x[j] ≤ variable_bounds[j].upper`, subject to one constraint per row of `a`:
/// row `i` is a ranged linear row when `rows[i]` is `Constraint::Linear`, and
/// one coordinate of a cone block when it is `Constraint::Cone`. Consecutive
/// cone rows sharing a `block` form that block, whose slacks lie in
/// `cones[block]` and whose right-hand sides are the tags' `rhs` values.
///
/// `p` holds the upper triangle of the symmetric positive semidefinite `P`,
/// including its diagonal; `None` means a linear objective. A reduced
/// problem's `p` is again an upper triangle. Empty `variable_bounds` means
/// every variable is free. Presolve works on its own copies of both matrices;
/// an unchanged problem hands the caller's buffers back untouched.
///
/// Sparse columns must have sorted unique row indices and finite
/// coefficients, cone dimensions must agree with their row blocks, and
/// dimensions and nonzero counts must fit below `u32::MAX`. No validation is
/// performed by `presolve`.
#[derive(Clone, Debug)]
pub struct Problem {
    pub p: Option<CscMatrix>,
    pub c: Vec<f64>,
    pub c0: f64,
    pub a: CscMatrix,
    pub rows: Vec<Constraint>,
    pub variable_bounds: Vec<Bounds>,
    pub cones: Vec<Cone>,
}
pub(crate) fn row_indices(rows: &[Constraint]) -> (Vec<usize>, Vec<usize>) {
    let mut linear = Vec::new();
    let mut conic = Vec::new();
    for (i, row) in rows.iter().enumerate() {
        match row {
            Constraint::Linear(_) => linear.push(i),
            Constraint::Cone { .. } => conic.push(i),
        }
    }
    (linear, conic)
}
impl Problem {
    pub fn variable_count(&self) -> usize {
        self.c.len()
    }
    pub fn row_count(&self) -> usize {
        self.rows.len()
    }
    /// Bounds of variable `j`, treating empty bounds as free.
    pub fn variable_bounds(&self, j: usize) -> Bounds {
        if self.variable_bounds.is_empty() {
            Bounds::FREE
        } else {
            self.variable_bounds[j]
        }
    }
    /// The working objective: a symmetric copy of the Hessian and the moved
    /// linear coefficients.
    pub(crate) fn take_objective(&mut self) -> Objective {
        Objective {
            p: self.p.as_ref().map_or_else(
                || SymmetricMatrix::zeros(self.c.len()),
                |p| SymmetricMatrix::from_upper_columns(p.columns(), |j| p.as_ref().column(j)),
            ),
            c: std::mem::take(&mut self.c),
            constant: self.c0,
            scratch: Default::default(),
        }
    }
}
