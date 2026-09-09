//! Sparse interchange matrices and views. Editable storage is internal.
mod csc;
pub(crate) mod linked;
pub(crate) mod quadratic;
pub(crate) mod sparse;

#[cfg(test)]
pub(crate) use csc::test_matrix;
pub use csc::{CscMatrix, CscMatrixRef, MatrixError};
pub use quadratic::QuadraticRef;
