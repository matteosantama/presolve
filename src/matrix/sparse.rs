// SPDX-License-Identifier: Apache-2.0
// Modified for this library; copyright and attribution notices are in NOTICE.
//! Editable sparse storage with sorted row and column vectors.
//! Indices remain stable until final packing.

#[cfg(test)]
use crate::matrix::CscMatrix;

pub(crate) type Entries = Vec<(usize, f64)>;

#[derive(Clone, Debug)]
pub(crate) struct SparseMatrix {
    rows: Vec<Entries>,
    columns: Vec<Entries>,
    nnz: usize,
    pub revision: usize,
}

impl SparseMatrix {
    pub fn zeros(rows: usize, columns: usize) -> Self {
        Self {
            rows: vec![Vec::new(); rows],
            columns: vec![Vec::new(); columns],
            nnz: 0,
            revision: 0,
        }
    }

    #[cfg(test)]
    pub fn from_matrix(matrix: &CscMatrix) -> Self {
        let mut out = Self::zeros(matrix.rows(), matrix.columns());
        for j in 0..matrix.columns() {
            for (i, v) in matrix.as_ref().column(j) {
                out.set(i, j, v);
            }
        }
        out.revision = 0;
        out
    }
    pub fn from_upper_columns<I: Iterator<Item = (usize, f64)>>(
        n: usize,
        column: impl Fn(usize) -> I,
    ) -> Self {
        let mut lengths = vec![0; n];
        for j in 0..n {
            for (i, v) in column(j) {
                if v != 0.0 {
                    lengths[i] += 1;
                    if i != j {
                        lengths[j] += 1;
                    }
                }
            }
        }
        let mut rows: Vec<Entries> = lengths.into_iter().map(Vec::with_capacity).collect();
        for j in 0..n {
            for (i, v) in column(j) {
                if v != 0.0 {
                    rows[i].push((j, v));
                    if i != j {
                        rows[j].push((i, v));
                    }
                }
            }
        }
        // Lower entries arrive before the diagonal and upper entries, in order.
        let nnz = rows.iter().map(Vec::len).sum();
        Self {
            columns: rows.clone(),
            rows,
            nnz,
            revision: 0,
        }
    }

    pub fn add_row(&mut self) -> usize {
        self.revision += 1;
        let i = self.rows.len();
        self.rows.push(Vec::new());
        i
    }

    pub fn add_column(&mut self) -> usize {
        self.revision += 1;
        let j = self.columns.len();
        self.columns.push(Vec::new());
        j
    }

    pub fn nnz(&self) -> usize {
        self.nnz
    }
    pub fn row(&self, row: usize) -> &[(usize, f64)] {
        &self.rows[row]
    }
    pub fn column(&self, col: usize) -> &[(usize, f64)] {
        &self.columns[col]
    }

    pub fn get(&self, row: usize, column: usize) -> f64 {
        self.rows[row]
            .binary_search_by_key(&column, |&(j, _)| j)
            .map_or(0.0, |index| self.rows[row][index].1)
    }

    pub fn set(&mut self, row: usize, column: usize, value: f64) {
        assert!(value.is_finite());
        let old = self.get(row, column);
        if old == value {
            return;
        }
        self.revision += 1;
        Self::set_entry(&mut self.rows[row], column, value);
        Self::set_entry(&mut self.columns[column], row, value);
        if old == 0.0 {
            self.nnz += 1;
        }
        if value == 0.0 {
            self.nnz -= 1;
        }
    }

    fn set_entry(entries: &mut Entries, index: usize, value: f64) {
        match entries.binary_search_by_key(&index, |&(i, _)| i) {
            Ok(at) if value == 0.0 => {
                entries.remove(at);
            }
            Ok(at) => entries[at].1 = value,
            Err(at) if value != 0.0 => entries.insert(at, (index, value)),
            Err(_) => {}
        }
    }

    pub fn remove_row(&mut self, row: usize) -> Entries {
        self.revision += 1;
        let entries = std::mem::take(&mut self.rows[row]);
        for &(column, _) in &entries {
            Self::set_entry(&mut self.columns[column], row, 0.0);
        }
        self.nnz -= entries.len();
        entries
    }

    pub fn remove_column(&mut self, column: usize) -> Entries {
        self.revision += 1;
        let entries = std::mem::take(&mut self.columns[column]);
        for &(row, _) in &entries {
            Self::set_entry(&mut self.rows[row], column, 0.0);
        }
        self.nnz -= entries.len();
        entries
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sparse_views_stay_consistent_through_fill_cancellation_and_deletion() {
        let mut a = SparseMatrix::zeros(2, 3);
        a.set(0, 0, 2.0);
        a.set(0, 2, 3.0);
        a.set(1, 0, -4.0);
        a.set(1, 1, 5.0);
        a.set(1, 0, 0.0);
        a.set(1, 2, 6.0);
        assert_eq!(a.row(1), [(1, 5.0), (2, 6.0)]);
        assert_eq!(a.column(0), [(0, 2.0)]);
        assert_eq!(a.column(2), [(0, 3.0), (1, 6.0)]);
        assert_eq!(a.remove_column(2), [(0, 3.0), (1, 6.0)]);
        assert_eq!(a.remove_row(0), [(0, 2.0)]);
        a.set(0, 2, 7.0);
        assert_eq!(a.column(2), [(0, 7.0)]);
        assert_eq!(a.nnz(), 2);
    }
}
