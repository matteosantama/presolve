// SPDX-License-Identifier: Apache-2.0
// Modified for this library; copyright and attribution notices are in NOTICE.
//! Editable symmetric Hessian storage with one sorted adjacency vector per variable.
//! Indices remain stable until final packing.

#[cfg(test)]
use crate::matrix::CscMatrix;

pub(crate) type Entries = Vec<(usize, f64)>;

#[derive(Clone, Debug)]
pub(crate) struct SymmetricMatrix {
    rows: Vec<Entries>,
    nnz: usize,
    pub revision: usize,
}

impl SymmetricMatrix {
    pub fn zeros(n: usize) -> Self {
        Self {
            rows: vec![Vec::new(); n],
            nnz: 0,
            revision: 0,
        }
    }

    #[cfg(test)]
    pub fn from_matrix(matrix: &CscMatrix) -> Self {
        assert_eq!(matrix.rows(), matrix.columns());
        Self::from_upper_columns(matrix.columns(), |j| {
            matrix.as_ref().column(j).filter(move |&(i, _)| i <= j)
        })
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
            rows,
            nnz,
            revision: 0,
        }
    }

    pub fn add_variable(&mut self) -> usize {
        self.revision += 1;
        let j = self.rows.len();
        self.rows.push(Vec::new());
        j
    }

    pub fn nnz(&self) -> usize {
        self.nnz
    }
    pub fn row(&self, row: usize) -> &[(usize, f64)] {
        &self.rows[row]
    }
    pub fn column(&self, col: usize) -> &[(usize, f64)] {
        &self.rows[col]
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
        if row != column {
            Self::set_entry(&mut self.rows[column], row, value);
        }
        let count = if row == column { 1 } else { 2 };
        if old == 0.0 {
            self.nnz += count;
        }
        if value == 0.0 {
            self.nnz -= count;
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

    /// Remove both symmetric views while keeping all other variable IDs stable.
    pub fn remove_variable(&mut self, column: usize) -> Entries {
        self.revision += 1;
        let entries = std::mem::take(&mut self.rows[column]);
        for &(row, _) in &entries {
            if row != column {
                Self::set_entry(&mut self.rows[row], column, 0.0);
                self.nnz -= 1;
            }
        }
        self.nnz -= entries.len();
        entries
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn symmetric_views_survive_fill_cancellation_removal_and_growth() {
        let mut p = SymmetricMatrix::zeros(3);
        p.set(0, 0, 2.0);
        p.set(0, 2, 3.0);
        p.set(1, 0, -4.0);
        p.set(1, 1, 5.0);
        p.set(0, 1, 0.0);
        p.set(1, 2, 6.0);
        assert_eq!(p.nnz(), 6);
        assert_eq!(p.row(0), [(0, 2.0), (2, 3.0)]);
        assert_eq!(p.column(2), [(0, 3.0), (1, 6.0)]);
        assert!(std::ptr::eq(p.row(2), p.column(2)));
        let revision = p.revision;
        p.set(2, 0, 3.0);
        assert_eq!(p.revision, revision);
        assert_eq!(p.remove_variable(2), [(0, 3.0), (1, 6.0)]);
        assert_eq!(p.nnz(), 2);
        assert_eq!(p.remove_variable(0), [(0, 2.0)]);
        assert_eq!(p.nnz(), 1);
        assert_eq!(p.add_variable(), 3);
        p.set(0, 3, 7.0);
        assert_eq!(p.row(3), [(0, 7.0)]);
        assert_eq!(p.nnz(), 3);
        let cloned = p.clone();
        p.remove_variable(3);
        assert_eq!(cloned.get(3, 0), 7.0);
        assert_eq!(p.nnz(), 1);
    }

    #[test]
    fn upper_triangle_builds_sorted_full_adjacency_without_duplicate_payloads() {
        let columns = [
            vec![(0, 2.0)],
            vec![(0, -1.0)],
            vec![(0, 0.0), (1, 3.0), (2, 4.0)],
        ];
        let p = SymmetricMatrix::from_upper_columns(3, |j| columns[j].iter().copied());
        assert_eq!(p.row(0), [(0, 2.0), (1, -1.0)]);
        assert_eq!(p.row(1), [(0, -1.0), (2, 3.0)]);
        assert_eq!(p.row(2), [(1, 3.0), (2, 4.0)]);
        assert_eq!(p.nnz(), 6);
        assert_eq!(p.revision, 0);
    }
}
