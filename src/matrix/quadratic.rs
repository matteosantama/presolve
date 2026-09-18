//! Lazy quadratic export. The working Hessian remains owned until requested.
use crate::matrix::CscMatrix;
use crate::matrix::sparse::SymmetricMatrix;
use std::sync::Arc;

/// Upper triangle of a symmetric Hessian, in either the caller's CSC form or a
/// reduced model's working storage. Packing to CSC happens only on request.
#[derive(Clone, Debug)]
pub struct QuadraticMatrix(pub(crate) Quadratic);

#[derive(Clone, Debug)]
pub(crate) enum Quadratic {
    Csc(CscMatrix),
    Sparse(Box<SparseQuadratic>),
}
#[derive(Clone, Debug)]
pub(crate) struct SparseQuadratic {
    matrix: SymmetricMatrix,
    compact_to_stable_columns: Arc<Vec<usize>>,
    stable_to_compact_columns: Arc<Vec<usize>>,
}
impl From<CscMatrix> for QuadraticMatrix {
    fn from(matrix: CscMatrix) -> Self {
        Self(Quadratic::Csc(matrix))
    }
}
impl QuadraticMatrix {
    pub fn columns(&self) -> usize {
        match &self.0 {
            Quadratic::Csc(p) => p.columns(),
            Quadratic::Sparse(p) => p.compact_to_stable_columns.len(),
        }
    }
    /// Borrow existing CSC buffers if available; this never triggers packing.
    pub fn as_csc(&self) -> Option<&CscMatrix> {
        match &self.0 {
            Quadratic::Csc(p) => Some(p),
            Quadratic::Sparse(_) => None,
        }
    }
    /// Upper-triangular column `j`, independent of the storage form.
    pub fn column(&self, j: usize) -> impl Iterator<Item = (usize, f64)> + '_ {
        self.0.column(j)
    }
    /// Move existing CSC buffers or pack the working storage once.
    pub fn into_csc(self) -> CscMatrix {
        self.0.into_csc()
    }
}
impl Quadratic {
    fn column(&self, j: usize) -> impl Iterator<Item = (usize, f64)> + '_ {
        let (csc, sparse) = match self {
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
    pub fn from_sparse(
        matrix: SymmetricMatrix,
        compact_to_stable_columns: Arc<Vec<usize>>,
        stable_to_compact_columns: Arc<Vec<usize>>,
    ) -> Self {
        Self::Sparse(Box::new(SparseQuadratic {
            matrix,
            compact_to_stable_columns,
            stable_to_compact_columns,
        }))
    }
    pub fn working(&self) -> SymmetricMatrix {
        match self {
            Self::Csc(p) => {
                SymmetricMatrix::from_upper_columns(p.columns(), |j| p.as_ref().column(j))
            }
            Self::Sparse(p) => {
                SymmetricMatrix::from_upper_columns(p.compact_to_stable_columns.len(), |j| {
                    self.column(j)
                })
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
