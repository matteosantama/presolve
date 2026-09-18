//! Owned native problem data with explicit, consuming CSC export.
mod bounds;
mod cone;
mod conic;

pub use bounds::Bounds;
pub use cone::Cone;
pub(crate) use cone::Membership;
pub use conic::{ConicData, ConicExport, ConicMap};

use crate::matrix::{CscMatrix, QuadraticMatrix, linked::LinkedMatrix};
use std::sync::Arc;

/// One row of the shared constraint matrix. Cone coordinates in a block are
/// consecutive; `block` indexes the cone list, independently of opaque IDs.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Constraint {
    Linear(Bounds),
    Cone { rhs: f64, block: usize },
}

/// A native problem. The two matrices are opaque: a caller's CSC input is
/// kept as is, and a reduced problem keeps the working model's editable
/// storage until CSC is requested. Empty variable bounds mean free variables.
/// `p` stores the upper triangle of the symmetric PSD Hessian, including its
/// diagonal. `a` contains all rows; cone dimensions must agree with
/// consecutive row blocks. Sparse columns must have sorted unique row indices
/// and finite coefficients, and dimensions and nonzero counts must fit below
/// u32::MAX. No validation is performed by `presolve`.
#[derive(Clone, Debug)]
pub struct Problem {
    pub p: Option<QuadraticMatrix>,
    pub c: Vec<f64>,
    pub objective_constant: f64,
    pub a: ConstraintMatrix,
    pub rows: Vec<Constraint>,
    pub variable_bounds: Vec<Bounds>,
    pub cones: Vec<Cone>,
}

/// Constraint storage: the caller's CSC matrix, sparse columns copied once
/// into the editable representation, or a reduced model's rows addressed
/// through survivor maps.
#[derive(Clone, Debug)]
pub struct ConstraintMatrix(pub(crate) Matrix);
#[derive(Clone, Debug)]
pub(crate) enum Matrix {
    Csc(CscMatrix),
    Linked {
        /// `None` while the working model owns the storage; see `working_matrix`.
        matrix: Option<LinkedMatrix>,
        columns: usize,
        // Compact row -> stable row, stable column -> compact column.
        compact_to_stable_rows: Vec<usize>,
        stable_to_compact_columns: Arc<Vec<usize>>,
    },
}
impl From<CscMatrix> for ConstraintMatrix {
    fn from(matrix: CscMatrix) -> Self {
        Self(Matrix::Csc(matrix))
    }
}
impl ConstraintMatrix {
    /// Build editable storage directly from a caller's sparse columns. This
    /// copies coefficients once and avoids an intermediate CSC copy when a
    /// solver retains its original data.
    pub fn from_columns<I: Iterator<Item = (usize, f64)>>(
        rows: usize,
        columns: usize,
        column: impl Fn(usize) -> I,
    ) -> Self {
        Self(Matrix::Linked {
            matrix: Some(LinkedMatrix::from_columns(rows, columns, column)),
            columns,
            compact_to_stable_rows: (0..rows).collect(),
            stable_to_compact_columns: Arc::new((0..columns).collect()),
        })
    }
    pub(crate) fn linked(
        matrix: LinkedMatrix,
        columns: usize,
        compact_to_stable_rows: Vec<usize>,
        stable_to_compact_columns: Arc<Vec<usize>>,
    ) -> Self {
        Self(Matrix::Linked {
            matrix: Some(matrix),
            columns,
            compact_to_stable_rows,
            stable_to_compact_columns,
        })
    }
    pub fn rows(&self) -> usize {
        match &self.0 {
            Matrix::Csc(a) => a.rows(),
            Matrix::Linked {
                compact_to_stable_rows,
                ..
            } => compact_to_stable_rows.len(),
        }
    }
    pub fn columns(&self) -> usize {
        match &self.0 {
            Matrix::Csc(a) => a.columns(),
            Matrix::Linked { columns, .. } => *columns,
        }
    }
    /// Borrow existing CSC buffers if available; this never triggers packing.
    pub fn as_csc(&self) -> Option<&CscMatrix> {
        match &self.0 {
            Matrix::Csc(a) => Some(a),
            Matrix::Linked { .. } => None,
        }
    }
    /// Visit a row. Reduced rows are visited in O(nnz(row)); CSC storage
    /// requires column lookups, so export need not use this slow input path.
    pub fn row(&self, i: usize) -> impl Iterator<Item = (usize, f64)> + '_ {
        let (csc, linked, map) = match &self.0 {
            Matrix::Csc(a) => (Some(a), None, &[][..]),
            Matrix::Linked {
                matrix,
                compact_to_stable_rows,
                stable_to_compact_columns,
                ..
            } => (
                None,
                Some(
                    matrix
                        .as_ref()
                        .expect("editable input belongs to the working model")
                        .row(compact_to_stable_rows[i])
                        .iter(),
                ),
                stable_to_compact_columns.as_slice(),
            ),
        };
        csc.into_iter()
            .flat_map(move |a| {
                (0..a.columns())
                    .filter_map(move |j| a.get(i, j).filter(|&v| v != 0.0).map(|v| (j, v)))
            })
            .chain(linked.into_iter().flatten().map(move |(j, v)| (map[j], v)))
    }
    /// Move existing CSC buffers or pack linked rows once. No solver-specific
    /// expansion of ranged constraints or bounds takes place here.
    pub fn into_csc(self) -> CscMatrix {
        match self.0 {
            Matrix::Csc(a) => a,
            Matrix::Linked {
                matrix,
                columns,
                compact_to_stable_rows,
                stable_to_compact_columns,
            } => pack_rows(compact_to_stable_rows.len(), columns, |i| {
                matrix
                    .as_ref()
                    .expect("working storage is returned before export")
                    .row(compact_to_stable_rows[i])
                    .iter()
                    .map(|(j, v)| (stable_to_compact_columns[j], v))
            }),
        }
    }
    pub(crate) fn restore_working_matrix(&mut self, matrix: LinkedMatrix) {
        if let Matrix::Linked {
            matrix: slot @ None,
            ..
        } = &mut self.0
        {
            *slot = Some(matrix);
        }
    }
    /// Editable storage for the working model. Input already in stable linked
    /// form is lent out and returned by `restore_working_matrix`; anything else
    /// is copied once.
    pub(crate) fn working_matrix(&mut self) -> LinkedMatrix {
        let (m, n) = (self.rows(), self.columns());
        let identity = |map: &[usize], len: usize| map.iter().copied().eq(0..len);
        match &mut self.0 {
            Matrix::Linked {
                matrix,
                compact_to_stable_rows,
                stable_to_compact_columns,
                ..
            } if identity(compact_to_stable_rows, m) && identity(stable_to_compact_columns, n) => {
                matrix.take().expect("working storage is lent at most once")
            }
            Matrix::Csc(a) => {
                LinkedMatrix::from_columns(a.rows(), a.columns(), |j| a.as_ref().column(j))
            }
            Matrix::Linked {
                matrix,
                compact_to_stable_rows,
                stable_to_compact_columns,
                ..
            } => {
                let matrix = matrix
                    .as_ref()
                    .expect("working storage is lent at most once");
                let a = pack_rows(compact_to_stable_rows.len(), n, |i| {
                    matrix
                        .row(compact_to_stable_rows[i])
                        .iter()
                        .map(|(j, v)| (stable_to_compact_columns[j], v))
                });
                LinkedMatrix::from_columns(a.rows(), a.columns(), |j| a.as_ref().column(j))
            }
        }
    }
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
    /// Pack both matrices into CSC storage, moving existing buffers. After
    /// this, `as_csc` on either matrix is `Some`.
    pub fn into_csc(self) -> Problem {
        Problem {
            p: self.p.map(|p| p.into_csc().into()),
            a: self.a.into_csc().into(),
            ..self
        }
    }
    pub(crate) fn take_objective(&mut self) -> crate::model::objective::Objective {
        crate::model::objective::Objective {
            p: self.p.as_ref().map_or_else(
                || crate::matrix::sparse::SymmetricMatrix::zeros(self.c.len()),
                |p| p.0.working(),
            ),
            c: std::mem::take(&mut self.c),
            constant: self.objective_constant,
            scratch: Default::default(),
        }
    }
}
pub(crate) fn pack_rows<I: Iterator<Item = (usize, f64)>>(
    m: usize,
    n: usize,
    visit: impl Fn(usize) -> I,
) -> CscMatrix {
    let mut pointers = vec![0; n + 1];
    for i in 0..m {
        for (j, _) in visit(i) {
            pointers[j + 1] += 1;
        }
    }
    for j in 0..n {
        pointers[j + 1] += pointers[j];
    }
    let mut next = pointers[..n].to_vec();
    let mut ri = vec![0; pointers[n]];
    let mut values = vec![0.; pointers[n]];
    for i in 0..m {
        for (j, v) in visit(i) {
            let at = next[j];
            ri[at] = i;
            values[at] = v;
            next[j] += 1;
        }
    }
    CscMatrix::from_parts(m, n, pointers, ri, values)
}
