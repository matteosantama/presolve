use thiserror::Error;

/// Errors produced while constructing a matrix.
#[derive(Clone, Debug, Error, PartialEq)]
#[non_exhaustive]
pub enum MatrixError {
    #[error("matrix dimensions overflow addressable memory")]
    DimensionOverflow,
    #[error("matrix values must all be finite")]
    NonFiniteValue,
    #[error("CSC column pointers must have length ncols + 1")]
    InvalidColumnPointerCount,
    #[error("CSC column pointers must begin at zero and be nondecreasing")]
    InvalidColumnPointers,
    #[error("CSC index and value arrays must have equal lengths")]
    MismatchedSparseArrays,
    #[error("CSC row index {row} is out of bounds for a matrix with {rows} rows")]
    RowIndexOutOfBounds { row: usize, rows: usize },
    #[error("CSC row indices must be strictly increasing within each column")]
    UnsortedRowIndices,
    #[error("triplet row, column, and value arrays must have equal lengths")]
    MismatchedTripletArrays,
    #[error("triplet column index {column} is out of bounds for {columns} columns")]
    ColumnIndexOutOfBounds { column: usize, columns: usize },
}

/// Owned CSC buffers. Consuming construction and export move the arrays.
#[derive(Clone, Debug, PartialEq)]
pub struct CscMatrix {
    rows: usize,
    columns: usize,
    values: Vec<f64>,
    column_pointers: Vec<usize>,
    row_indices: Vec<usize>,
}

impl CscMatrix {
    /// Adopt caller-validated CSC arrays without scanning them.
    pub fn from_parts(
        rows: usize,
        columns: usize,
        column_pointers: Vec<usize>,
        row_indices: Vec<usize>,
        values: Vec<f64>,
    ) -> Self {
        Self {
            rows,
            columns,
            column_pointers,
            row_indices,
            values,
        }
    }

    pub fn new(
        rows: usize,
        columns: usize,
        column_pointers: Vec<usize>,
        row_indices: Vec<usize>,
        values: Vec<f64>,
    ) -> Result<Self, MatrixError> {
        CscMatrixRef::new(rows, columns, &column_pointers, &row_indices, &values)?;
        Ok(Self::from_parts(
            rows,
            columns,
            column_pointers,
            row_indices,
            values,
        ))
    }

    pub fn from_triplets(
        rows: usize,
        columns: usize,
        row_indices: Vec<usize>,
        column_indices: Vec<usize>,
        values: Vec<f64>,
    ) -> Result<Self, MatrixError> {
        if row_indices.len() != column_indices.len() || row_indices.len() != values.len() {
            return Err(MatrixError::MismatchedTripletArrays);
        }

        let mut triplets = Vec::with_capacity(values.len());
        for ((row, column), value) in row_indices.into_iter().zip(column_indices).zip(values) {
            if row >= rows {
                return Err(MatrixError::RowIndexOutOfBounds { row, rows });
            }
            if column >= columns {
                return Err(MatrixError::ColumnIndexOutOfBounds { column, columns });
            }
            if !value.is_finite() {
                return Err(MatrixError::NonFiniteValue);
            }
            triplets.push((column, row, value));
        }
        triplets.sort_unstable_by_key(|&(column, row, _)| (column, row));

        let mut canonical: Vec<(usize, usize, f64)> = Vec::with_capacity(triplets.len());
        for (column, row, value) in triplets {
            if let Some((last_column, last_row, last_value)) = canonical.last_mut()
                && *last_column == column
                && *last_row == row
            {
                *last_value += value;
                if !last_value.is_finite() {
                    return Err(MatrixError::NonFiniteValue);
                }
                continue;
            }
            canonical.push((column, row, value));
        }

        let mut column_pointers = vec![
            0;
            columns
                .checked_add(1)
                .ok_or(MatrixError::DimensionOverflow)?
        ];
        let mut canonical_rows = Vec::new();
        let mut canonical_values = Vec::new();
        let mut next_column = 0;
        for (column, row, value) in canonical {
            while next_column < column {
                next_column += 1;
                column_pointers[next_column] = canonical_values.len();
            }
            if value != 0.0 {
                canonical_rows.push(row);
                canonical_values.push(value);
            }
        }
        while next_column < columns {
            next_column += 1;
            column_pointers[next_column] = canonical_values.len();
        }

        Self::new(
            rows,
            columns,
            column_pointers,
            canonical_rows,
            canonical_values,
        )
    }

    pub fn identity(size: usize) -> Result<Self, MatrixError> {
        size.checked_add(1).ok_or(MatrixError::DimensionOverflow)?;
        Self::new(
            size,
            size,
            (0..=size).collect(),
            (0..size).collect(),
            vec![1.0; size],
        )
    }

    pub fn zeros(rows: usize, columns: usize) -> Result<Self, MatrixError> {
        Self::new(
            rows,
            columns,
            vec![
                0;
                columns
                    .checked_add(1)
                    .ok_or(MatrixError::DimensionOverflow)?
            ],
            vec![],
            vec![],
        )
    }

    #[inline]
    pub fn rows(&self) -> usize {
        self.rows
    }

    #[inline]
    pub fn columns(&self) -> usize {
        self.columns
    }

    #[inline]
    pub fn column_pointers(&self) -> &[usize] {
        &self.column_pointers
    }

    #[inline]
    pub fn row_indices(&self) -> &[usize] {
        &self.row_indices
    }

    #[inline]
    pub fn values(&self) -> &[f64] {
        &self.values
    }

    /// Mutate coefficients without changing the sparsity pattern.
    pub fn values_mut(&mut self) -> &mut [f64] {
        &mut self.values
    }

    #[inline]
    pub fn get(&self, row: usize, column: usize) -> Option<f64> {
        if row >= self.rows || column >= self.columns {
            return None;
        }
        let start = self.column_pointers[column];
        let end = self.column_pointers[column + 1];
        let local = self.row_indices[start..end].binary_search(&row).ok()?;
        Some(self.values[start + local])
    }
}

#[cfg(test)]
pub(crate) fn matrix_mul(matrix: &CscMatrix, vector: &[f64]) -> Vec<f64> {
    let mut product = vec![0.0; matrix.rows()];
    for (column, value) in vector.iter().copied().enumerate() {
        for index in matrix.column_pointers()[column]..matrix.column_pointers()[column + 1] {
            product[matrix.row_indices()[index]] += matrix.values()[index] * value;
        }
    }
    product
}

#[cfg(test)]
pub(crate) fn matrix_transpose_mul(matrix: &CscMatrix, vector: &[f64]) -> Vec<f64> {
    let mut product = vec![0.0; matrix.columns()];
    for (column, result) in product.iter_mut().enumerate() {
        for index in matrix.column_pointers()[column]..matrix.column_pointers()[column + 1] {
            *result += matrix.values()[index] * vector[matrix.row_indices()[index]];
        }
    }
    product
}

/// Visit every stored entry in column order.
#[cfg(test)]
pub(crate) fn visit_entries(matrix: &CscMatrix, mut visit: impl FnMut(usize, usize, f64)) {
    for column in 0..matrix.columns() {
        for index in matrix.column_pointers()[column]..matrix.column_pointers()[column + 1] {
            visit(matrix.row_indices()[index], column, matrix.values()[index]);
        }
    }
}

/// Build a small sparse fixture from row-major values. Test code only.
#[cfg(test)]
pub(crate) fn test_matrix(
    rows: usize,
    columns: usize,
    values: Vec<f64>,
) -> Result<CscMatrix, MatrixError> {
    assert_eq!(values.len(), rows * columns);
    let mut ri = Vec::new();
    let mut ci = Vec::new();
    let mut data = Vec::new();
    for (i, v) in values.into_iter().enumerate() {
        if v != 0.0 {
            ri.push(i / columns);
            ci.push(i % columns);
            data.push(v);
        }
    }
    CscMatrix::from_triplets(rows, columns, ri, ci, data)
}

#[cfg(test)]
mod arithmetic_tests {
    use super::*;

    #[test]
    fn rectangular_products_handle_explicit_sparse_zeros() {
        let canonical: CscMatrix =
            crate::matrix::test_matrix(2, 3, vec![2.0, 0.0, -1.0, 0.0, 3.0, 4.0]).unwrap();
        let sparse: CscMatrix = CscMatrix::new(
            2,
            3,
            vec![0, 2, 3, 5],
            vec![0, 1, 1, 0, 1],
            vec![2.0, 0.0, 3.0, -1.0, 4.0],
        )
        .unwrap();
        for matrix in [&canonical, &sparse] {
            assert_eq!(matrix_mul(matrix, &[1.0, 2.0, 3.0]), vec![-1.0, 18.0]);
            assert_eq!(
                matrix_transpose_mul(matrix, &[2.0, -1.0]),
                vec![4.0, -3.0, -6.0]
            );
        }
        let mut entries = Vec::new();
        visit_entries(&sparse, |row, col, value| entries.push((row, col, value)));
        assert_eq!(
            entries,
            vec![
                (0, 0, 2.0),
                (1, 0, 0.0),
                (1, 1, 3.0),
                (0, 2, -1.0),
                (1, 2, 4.0)
            ]
        );
    }
}

/// A borrowed CSC matrix. Callers supply sorted, unique row indices per column.
#[derive(Clone, Copy, Debug)]
pub struct CscMatrixRef<'a> {
    rows: usize,
    columns: usize,
    column_pointers: &'a [usize],
    row_indices: &'a [usize],
    values: &'a [f64],
}
impl<'a> CscMatrixRef<'a> {
    /// Borrow caller-validated CSC arrays without scanning them.
    /// Invalid dimensions, indices, or values may panic or produce incorrect results.
    pub fn from_parts(
        rows: usize,
        columns: usize,
        column_pointers: &'a [usize],
        row_indices: &'a [usize],
        values: &'a [f64],
    ) -> Self {
        Self {
            rows,
            columns,
            column_pointers,
            row_indices,
            values,
        }
    }
    pub fn column(self, column: usize) -> impl ExactSizeIterator<Item = (usize, f64)> + 'a {
        let range = self.column_pointers[column]..self.column_pointers[column + 1];
        self.row_indices[range.clone()]
            .iter()
            .copied()
            .zip(self.values[range].iter().copied())
    }

    pub fn new(
        rows: usize,
        columns: usize,
        column_pointers: &'a [usize],
        row_indices: &'a [usize],
        values: &'a [f64],
    ) -> Result<Self, MatrixError> {
        if column_pointers.len() != columns.saturating_add(1) {
            return Err(MatrixError::InvalidColumnPointerCount);
        }
        if row_indices.len() != values.len() {
            return Err(MatrixError::MismatchedSparseArrays);
        }
        if column_pointers.first() != Some(&0)
            || column_pointers.last() != Some(&values.len())
            || column_pointers.windows(2).any(|pair| pair[0] > pair[1])
        {
            return Err(MatrixError::InvalidColumnPointers);
        }
        if !values.iter().all(|value| value.is_finite()) {
            return Err(MatrixError::NonFiniteValue);
        }
        for column in 0..columns {
            let start = column_pointers[column];
            let end = column_pointers[column + 1];
            let indices = &row_indices[start..end];
            for &row in indices {
                if row >= rows {
                    return Err(MatrixError::RowIndexOutOfBounds { row, rows });
                }
            }
            if indices.windows(2).any(|pair| pair[0] >= pair[1]) {
                return Err(MatrixError::UnsortedRowIndices);
            }
        }

        Ok(Self {
            rows,
            columns,
            column_pointers,
            row_indices,
            values,
        })
    }
    #[inline]
    pub fn rows(&self) -> usize {
        self.rows
    }
    #[inline]
    pub fn columns(&self) -> usize {
        self.columns
    }
    #[inline]
    pub fn column_pointers(&self) -> &'a [usize] {
        self.column_pointers
    }
    #[inline]
    pub fn row_indices(&self) -> &'a [usize] {
        self.row_indices
    }
    #[inline]
    pub fn values(&self) -> &'a [f64] {
        self.values
    }
    #[inline]
    pub fn get(&self, row: usize, column: usize) -> Option<f64> {
        if row >= self.rows || column >= self.columns {
            return None;
        }
        let start = self.column_pointers[column];
        let end = self.column_pointers[column + 1];
        let local = self.row_indices[start..end].binary_search(&row).ok()?;
        Some(self.values[start + local])
    }
    /// Check full-storage symmetry in O(nnz + n) work and O(n) scratch space.
    /// Pairs use the tolerance `1e-12 * (1 + max(abs(a), abs(b)))`.
    pub fn is_symmetric(&self) -> bool {
        if self.rows != self.columns {
            return false;
        }
        // Columns are visited in increasing order, so each transposed-column
        // cursor advances monotonically instead of binary-searching each pair.
        let mut next = self.column_pointers[..self.columns].to_vec();
        for j in 0..self.columns {
            for k in self.column_pointers[j]..self.column_pointers[j + 1] {
                let i = self.row_indices[k];
                let at = &mut next[i];
                let end = self.column_pointers[i + 1];
                while *at < end && self.row_indices[*at] < j {
                    *at += 1;
                }
                let other = if *at < end && self.row_indices[*at] == j {
                    self.values[*at]
                } else {
                    0.0
                };
                let value = self.values[k];
                if (value - other).abs() > 1e-12 * (1.0 + value.abs().max(other.abs())) {
                    return false;
                }
            }
        }
        true
    }
    pub fn to_owned(self) -> CscMatrix {
        CscMatrix::from_parts(
            self.rows,
            self.columns,
            self.column_pointers.to_vec(),
            self.row_indices.to_vec(),
            self.values.to_vec(),
        )
    }
}
impl CscMatrix {
    #[inline]
    pub fn as_ref(&self) -> CscMatrixRef<'_> {
        CscMatrixRef {
            rows: self.rows,
            columns: self.columns,
            column_pointers: &self.column_pointers,
            row_indices: &self.row_indices,
            values: &self.values,
        }
    }
    /// Transfer the arrays without copying their elements.
    pub fn into_parts(self) -> (Vec<usize>, Vec<usize>, Vec<f64>) {
        (self.column_pointers, self.row_indices, self.values)
    }
}
