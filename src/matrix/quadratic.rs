//! Lazy quadratic export. The working Hessian remains owned until requested.
use crate::matrix::sparse::SparseMatrix;
use crate::matrix::{CscMatrix, CscMatrixRef};
use std::sync::Arc;

#[derive(Clone, Debug)]
pub(crate) enum Quadratic {
    Csc(CscMatrix),
    Sparse(Box<SparseQuadratic>),
}
#[derive(Clone, Debug)]
pub(crate) struct SparseQuadratic {
    matrix: SparseMatrix,
    compact_to_stable_columns: Arc<Vec<usize>>,
    stable_to_compact_columns: Arc<Vec<usize>>,
}
/// Upper-triangular access independent of the reduced Hessian's storage.
#[derive(Clone, Copy)]
pub struct QuadraticRef<'a> {
    inner: &'a Quadratic,
}
impl<'a> QuadraticRef<'a> {
    pub fn columns(self) -> usize {
        match self.inner {
            Quadratic::Csc(p) => p.columns(),
            Quadratic::Sparse(p) => p.compact_to_stable_columns.len(),
        }
    }
    /// Borrow existing CSC buffers if available; this never triggers packing.
    pub fn as_csc(self) -> Option<CscMatrixRef<'a>> {
        match self.inner {
            Quadratic::Csc(p) => Some(p.as_ref()),
            Quadratic::Sparse(_) => None,
        }
    }
    pub fn column(self, j: usize) -> impl Iterator<Item = (usize, f64)> + 'a {
        let (csc, sparse) = match self.inner {
            Quadratic::Csc(p) => (Some(p.as_ref()), None),
            Quadratic::Sparse(p) => (None, Some(p.as_ref())),
        };
        csc.into_iter()
            .flat_map(move |p| p.column(j))
            .chain(sparse.into_iter().flat_map(move |p| {
                p.matrix
                    .column(p.compact_to_stable_columns[j])
                    .iter()
                    .filter_map(move |&(i, v)| {
                        let row = p.stable_to_compact_columns[i];
                        (row <= j).then_some((row, v))
                    })
            }))
    }
}
impl Quadratic {
    pub fn as_ref(&self) -> QuadraticRef<'_> {
        QuadraticRef { inner: self }
    }
    pub fn from_sparse(
        matrix: SparseMatrix,
        compact_to_stable_columns: Arc<Vec<usize>>,
        stable_to_compact_columns: Arc<Vec<usize>>,
    ) -> Self {
        Self::Sparse(Box::new(SparseQuadratic {
            matrix,
            compact_to_stable_columns,
            stable_to_compact_columns,
        }))
    }
    pub fn working(&self) -> SparseMatrix {
        match self {
            Self::Csc(p) => SparseMatrix::from_upper_columns(p.columns(), |j| p.as_ref().column(j)),
            Self::Sparse(_) => {
                let p = self.as_ref();
                SparseMatrix::from_upper_columns(p.columns(), |j| p.column(j))
            }
        }
    }
    pub fn into_csc(self) -> CscMatrix {
        match self {
            Self::Csc(p) => p,
            Self::Sparse(p) => {
                let n = p.compact_to_stable_columns.len();
                let nnz = p.matrix.nnz();
                let capacity = nnz.min(nnz / 2 + n);
                let mut pointers = Vec::with_capacity(n + 1);
                let mut rows = Vec::with_capacity(capacity);
                let mut values = Vec::with_capacity(capacity);
                pointers.push(0);
                for (j, &old) in p.compact_to_stable_columns.iter().enumerate() {
                    for &(i, v) in p.matrix.column(old) {
                        let row = p.stable_to_compact_columns[i];
                        if row <= j {
                            rows.push(row);
                            values.push(v);
                        }
                    }
                    pointers.push(values.len());
                }
                CscMatrix::from_parts(n, n, pointers, rows, values)
            }
        }
    }
}
