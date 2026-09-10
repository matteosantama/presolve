//! Owned native data and explicit, consuming CSC export.
mod bounds;
mod cone;
mod conic;

pub use bounds::Bounds;
pub use cone::Cone;
pub(crate) use cone::Membership;
pub use conic::{ConicData, ConicExport, ConicMap};

use crate::matrix::{CscMatrix, QuadraticRef, linked::LinkedMatrix, quadratic::Quadratic};
use std::sync::Arc;

/// One row of the shared constraint matrix. Cone coordinates in a block are
/// consecutive; `block` indexes the cone list, independently of opaque IDs.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Constraint {
    Linear(Bounds),
    Cone { rhs: f64, block: usize },
}

/// Caller-validated interchange data. Constraint storage defaults to CSC;
/// `ConstraintMatrix` also accepts sparse columns directly. `p` stores the upper triangle of the
/// symmetric PSD Hessian, including its diagonal. `a` contains all rows.
/// Empty variable bounds mean free variables. Cone dimensions must agree with
/// consecutive row blocks. Sparse columns must have sorted unique row indices
/// and finite coefficients. Dimensions and internal nonzero counts must fit
/// below u32::MAX. No validation is performed by `From` or `presolve`.
#[derive(Clone, Debug)]
pub struct ProblemData<M = CscMatrix> {
    pub p: Option<CscMatrix>,
    pub c: Vec<f64>,
    pub objective_constant: f64,
    pub a: M,
    pub rows: Vec<Constraint>,
    pub variable_bounds: Vec<Bounds>,
    pub cones: Vec<Cone>,
}

/// Owned constraint storage, built directly from a caller's sparse columns.
/// This construction copies coefficients once into the editable representation,
/// avoiding an intermediate CSC copy when a solver retains its original data.
#[derive(Clone, Debug)]
pub struct ConstraintMatrix(Matrix);
impl From<CscMatrix> for ConstraintMatrix {
    fn from(matrix: CscMatrix) -> Self {
        Self(Matrix::Csc(matrix))
    }
}
impl ConstraintMatrix {
    pub fn from_columns<I: Iterator<Item = (usize, f64)>>(
        rows: usize,
        columns: usize,
        column: impl Fn(usize) -> I,
    ) -> Self {
        Self(Matrix::Linked {
            matrix: LinkedMatrix::from_columns(rows, columns, column),
            compact_to_stable_rows: (0..rows).collect(),
            stable_to_compact_columns: Arc::new((0..columns).collect()),
        })
    }
}

/// A native owned problem. Sparse layout is private: a reduced problem retains
/// editable row storage until an explicit CSC export is requested.
#[derive(Clone, Debug)]
pub struct Problem {
    pub(crate) p: Option<Quadratic>,
    pub(crate) c: Vec<f64>,
    pub(crate) objective_constant: f64,
    pub(crate) a: Matrix,
    pub(crate) rows: Vec<Constraint>,
    pub(crate) variable_bounds: Vec<Bounds>,
    pub(crate) cones: Vec<Cone>,
    pub(crate) linear_rows: Vec<usize>,
    pub(crate) conic_rows: Vec<usize>,
}
#[derive(Clone, Debug)]
pub(crate) enum Matrix {
    Vacant,
    // Temporary state while already-editable input is owned by Model.
    Moved {
        compact_to_stable_rows: Vec<usize>,
        stable_to_compact_columns: Arc<Vec<usize>>,
    },
    Csc(CscMatrix),
    Linked {
        matrix: LinkedMatrix,
        // Compact row -> stable row, stable column -> compact column.
        compact_to_stable_rows: Vec<usize>,
        stable_to_compact_columns: Arc<Vec<usize>>,
    },
}
impl<M: Into<ConstraintMatrix>> From<ProblemData<M>> for Problem {
    fn from(data: ProblemData<M>) -> Self {
        let (linear_rows, conic_rows) = row_indices(&data.rows);
        Self {
            p: data.p.map(Quadratic::Csc),
            c: data.c,
            objective_constant: data.objective_constant,
            a: data.a.into().0,
            rows: data.rows,
            variable_bounds: data.variable_bounds,
            cones: data.cones,
            linear_rows,
            conic_rows,
        }
    }
}
impl From<Problem> for ProblemData {
    fn from(p: Problem) -> Self {
        p.into_csc()
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
    pub fn linear_row_count(&self) -> usize {
        self.linear_rows.len()
    }
    pub fn conic_row_count(&self) -> usize {
        self.conic_rows.len()
    }
    pub fn p(&self) -> Option<QuadraticRef<'_>> {
        self.p.as_ref().map(Quadratic::as_ref)
    }
    pub fn c(&self) -> &[f64] {
        &self.c
    }
    pub fn objective_constant(&self) -> f64 {
        self.objective_constant
    }
    pub fn constraints(&self) -> &[Constraint] {
        &self.rows
    }
    pub fn cones(&self) -> &[Cone] {
        &self.cones
    }
    pub fn variable_bounds(&self, j: usize) -> Bounds {
        if self.variable_bounds.is_empty() {
            Bounds::FREE
        } else {
            self.variable_bounds[j]
        }
    }
    pub fn row_bounds(&self, i: usize) -> Bounds {
        match self.rows[self.linear_rows[i]] {
            Constraint::Linear(b) => b,
            _ => unreachable!(),
        }
    }
    pub fn conic_rhs(&self, i: usize) -> f64 {
        match self.rows[self.conic_rows[i]] {
            Constraint::Cone { rhs, .. } => rhs,
            _ => unreachable!(),
        }
    }
    /// Visit a row in the shared matrix. Reduced rows are visited in O(nnz(row));
    /// CSC inputs require column lookups. Export need not use this slow input path.
    pub fn matrix_row(&self, i: usize) -> impl Iterator<Item = (usize, f64)> + '_ {
        let (csc, linked, map) = match &self.a {
            Matrix::Csc(a) => (Some(a), None, &[][..]),
            Matrix::Linked {
                matrix,
                compact_to_stable_rows,
                stable_to_compact_columns,
            } => (
                None,
                Some(matrix.row(compact_to_stable_rows[i]).iter()),
                stable_to_compact_columns.as_slice(),
            ),
            Matrix::Moved { .. } | Matrix::Vacant => {
                unreachable!("editable input belongs to the working model")
            }
        };
        csc.into_iter()
            .flat_map(move |a| {
                (0..a.columns())
                    .filter_map(move |j| a.get(i, j).filter(|&v| v != 0.0).map(|v| (j, v)))
            })
            .chain(linked.into_iter().flatten().map(move |(j, v)| (map[j], v)))
    }
    pub fn row(&self, i: usize) -> impl Iterator<Item = (usize, f64)> + '_ {
        self.matrix_row(self.linear_rows[i])
    }
    pub fn conic_row(&self, i: usize) -> impl Iterator<Item = (usize, f64)> + '_ {
        self.matrix_row(self.conic_rows[i])
    }
    /// Consume the objective and discard constraints, for callers that already
    /// exported their own constraint representation using row iterators.
    pub fn into_objective(self) -> (Option<CscMatrix>, Vec<f64>, f64) {
        (
            self.p.map(Quadratic::into_csc),
            self.c,
            self.objective_constant,
        )
    }
    /// Move existing CSC buffers or pack linked rows once. No solver-specific
    /// expansion of ranged constraints or bounds takes place here.
    pub fn into_csc(self) -> ProblemData {
        let n = self.variable_count();
        let a = match self.a {
            Matrix::Moved { .. } | Matrix::Vacant => {
                unreachable!("working storage is returned before export")
            }
            Matrix::Csc(a) => a,
            Matrix::Linked {
                matrix,
                compact_to_stable_rows,
                stable_to_compact_columns,
            } => pack_rows(compact_to_stable_rows.len(), n, |i| {
                matrix
                    .row(compact_to_stable_rows[i])
                    .iter()
                    .map(|(j, v)| (stable_to_compact_columns[j], v))
            }),
        };
        ProblemData {
            p: self.p.map(Quadratic::into_csc),
            c: self.c,
            objective_constant: self.objective_constant,
            a,
            rows: self.rows,
            variable_bounds: self.variable_bounds,
            cones: self.cones,
        }
    }
    pub(crate) fn take_objective(&mut self) -> crate::core::objective::Objective {
        crate::core::objective::Objective {
            p: self.p.as_ref().map_or_else(
                || crate::matrix::sparse::SymmetricMatrix::zeros(self.c.len()),
                Quadratic::working,
            ),
            c: std::mem::take(&mut self.c),
            constant: self.objective_constant,
            scratch: Default::default(),
        }
    }
    pub(crate) fn restore_working_matrix(&mut self, matrix: LinkedMatrix) {
        if matches!(self.a, Matrix::Moved { .. }) {
            let Matrix::Moved {
                compact_to_stable_rows,
                stable_to_compact_columns,
                ..
            } = std::mem::replace(&mut self.a, Matrix::Vacant)
            else {
                unreachable!()
            };
            self.a = Matrix::Linked {
                matrix,
                compact_to_stable_rows,
                stable_to_compact_columns,
            };
        }
    }
    pub(crate) fn working_matrix(&mut self) -> LinkedMatrix {
        let direct = match &self.a {
            Matrix::Linked {
                compact_to_stable_rows,
                stable_to_compact_columns,
                ..
            } => {
                compact_to_stable_rows
                    .iter()
                    .copied()
                    .eq(0..self.row_count())
                    && stable_to_compact_columns
                        .iter()
                        .copied()
                        .eq(0..self.variable_count())
            }
            _ => false,
        };
        if direct {
            let Matrix::Linked {
                matrix,
                compact_to_stable_rows,
                stable_to_compact_columns,
            } = std::mem::replace(&mut self.a, Matrix::Vacant)
            else {
                unreachable!()
            };
            self.a = Matrix::Moved {
                compact_to_stable_rows,
                stable_to_compact_columns,
            };
            return matrix;
        }
        match &self.a {
            Matrix::Moved { .. } | Matrix::Vacant => unreachable!(),
            Matrix::Csc(a) => {
                LinkedMatrix::from_columns(a.rows(), a.columns(), |j| a.as_ref().column(j))
            }
            Matrix::Linked {
                matrix,
                compact_to_stable_rows,
                stable_to_compact_columns,
            } => {
                let a = pack_rows(compact_to_stable_rows.len(), self.variable_count(), |i| {
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
