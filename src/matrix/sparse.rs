// SPDX-License-Identifier: Apache-2.0
// Modified for this library; copyright and attribution notices are in NOTICE.
//! Editable symmetric Hessian. Off-diagonal coefficients are stored once;
//! sorted adjacency lists contain compact references to the shared coefficients.

use crate::matrix::CscMatrix;

pub(crate) type Entries = Vec<(usize, f64)>;
const NONE: u32 = u32::MAX;

// Row degrees fit the public u32 dimension limit. Arena offsets remain usize
// because spare capacity and released slots can exceed the live nonzero count.
#[derive(Clone, Copy, Debug, Default)]
struct Slot {
    start: usize,
    len: u32,
    capacity: u32,
}

/// Sorted support stays contiguous; only coefficient reads are indirect.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Link {
    column: u32,
    coefficient: u32,
}
const EMPTY_LINK: Link = Link {
    column: NONE,
    coefficient: NONE,
};

#[derive(Clone, Debug)]
pub(crate) struct SymmetricMatrix {
    values: Vec<f64>,
    /// Released value slots store the next free ID in their raw bits. They
    /// are unreachable from all live adjacency lists.
    free: u32,
    free_count: usize,
    entries: Vec<Link>,
    slots: Vec<Slot>,
    /// Allocated only when the matrix first acquires a nonzero diagonal.
    diagonal: Vec<f64>,
    dead: usize,
    nnz: usize,
    pub revision: usize,
}

#[derive(Clone, Copy)]
pub(crate) struct View<'a> {
    matrix: &'a SymmetricMatrix,
    row: usize,
}
impl<'a> View<'a> {
    #[inline]
    pub fn len(self) -> usize {
        self.matrix.slots[self.row].len as usize
            + usize::from(self.matrix.diagonal(self.row) != 0.0)
    }
    #[inline]
    pub fn is_empty(self) -> bool {
        self.len() == 0
    }
    #[inline]
    pub fn iter(self) -> Iter<'a> {
        let slot = self.matrix.slots[self.row];
        Iter {
            values: &self.matrix.values,
            entries: self.matrix.entries[slot.start..slot.start + slot.len as usize].iter(),
            row: self.row,
            diagonal: self.matrix.diagonal(self.row),
        }
    }
    pub fn to_vec(self) -> Entries {
        self.iter().collect()
    }
}
impl<'a> IntoIterator for View<'a> {
    type Item = (usize, f64);
    type IntoIter = Iter<'a>;
    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}

#[derive(Clone)]
pub(crate) struct Iter<'a> {
    values: &'a [f64],
    entries: std::slice::Iter<'a, Link>,
    row: usize,
    diagonal: f64,
}
impl Iter<'_> {
    /// Look up nondecreasing columns, skipping support indices without fetching
    /// their coefficients. Recreate the cursor before a backwards request.
    #[inline]
    pub fn advance_to(&mut self, column: usize) -> f64 {
        while self
            .entries
            .as_slice()
            .first()
            .is_some_and(|link| (link.column as usize) < column)
        {
            self.entries.next();
        }
        if column == self.row {
            return self.diagonal;
        }
        if column > self.row {
            self.diagonal = 0.0;
        }
        self.entries
            .as_slice()
            .first()
            .filter(|link| link.column as usize == column)
            .map_or(0.0, |link| self.values[link.coefficient as usize])
    }
}
impl Iterator for Iter<'_> {
    type Item = (usize, f64);
    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        let link = self.entries.as_slice().first();
        if self.diagonal != 0.0 && link.is_none_or(|e| e.column as usize > self.row) {
            return Some((self.row, std::mem::replace(&mut self.diagonal, 0.0)));
        }
        let link = self.entries.next()?;
        Some((link.column as usize, self.values[link.coefficient as usize]))
    }
    fn size_hint(&self) -> (usize, Option<usize>) {
        let len = self.entries.len() + usize::from(self.diagonal != 0.0);
        (len, Some(len))
    }
}
impl ExactSizeIterator for Iter<'_> {}

impl SymmetricMatrix {
    pub fn zeros(n: usize) -> Self {
        Self {
            values: Vec::new(),
            free: NONE,
            free_count: 0,
            entries: Vec::new(),
            slots: vec![Slot::default(); n],
            diagonal: Vec::new(),
            dead: 0,
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
        let mut matrix = Self::zeros(n);
        for j in 0..n {
            for (i, v) in column(j) {
                if v == 0.0 {
                    continue;
                }
                if i == j {
                    if matrix.diagonal.is_empty() {
                        matrix.diagonal.resize(n, 0.0);
                    }
                    matrix.diagonal[j] = v;
                    matrix.nnz += 1;
                } else {
                    matrix.slots[i].len += 1;
                    matrix.slots[j].len += 1;
                    matrix.nnz += 2;
                }
            }
        }
        let mut count = 0;
        for slot in &mut matrix.slots {
            slot.start = count;
            slot.capacity = slot.len;
            count += slot.len as usize;
        }
        matrix.entries.resize(count, EMPTY_LINK);
        matrix.values.reserve_exact(count / 2);
        let mut next: Vec<_> = matrix.slots.iter().map(|s| s.start).collect();
        for j in 0..n {
            for (i, v) in column(j) {
                if v != 0.0 && i != j {
                    let id = matrix.values.len() as u32;
                    matrix.values.push(v);
                    matrix.entries[next[i]] = Link {
                        column: j as u32,
                        coefficient: id,
                    };
                    matrix.entries[next[j]] = Link {
                        column: i as u32,
                        coefficient: id,
                    };
                    next[i] += 1;
                    next[j] += 1;
                }
            }
        }
        matrix
    }
    #[inline]
    fn diagonal(&self, row: usize) -> f64 {
        self.diagonal.get(row).copied().unwrap_or(0.0)
    }
    pub fn add_variable(&mut self) -> usize {
        self.revision += 1;
        let j = self.slots.len();
        self.slots.push(Slot::default());
        if !self.diagonal.is_empty() {
            self.diagonal.push(0.0);
        }
        j
    }
    pub fn pack_upper(
        &self,
        compact_to_stable: &[usize],
        stable_to_compact: &[usize],
    ) -> CscMatrix {
        let n = compact_to_stable.len();
        let capacity = self.nnz.min(self.nnz / 2 + n);
        let mut pointers = Vec::with_capacity(n + 1);
        let mut rows = Vec::with_capacity(capacity);
        let mut values = Vec::with_capacity(capacity);
        pointers.push(0);
        // Presolve retains stable variable order. Stop at the diagonal in
        // that common case; arbitrary output orders still use the full view.
        let ordered = compact_to_stable.windows(2).all(|pair| pair[0] < pair[1]);
        for (j, &old) in compact_to_stable.iter().enumerate() {
            for (i, v) in self.column(old) {
                if ordered && i > old {
                    break;
                }
                let row = stable_to_compact[i];
                if row <= j {
                    rows.push(row);
                    values.push(v);
                }
            }
            pointers.push(values.len());
        }
        CscMatrix::from_parts(n, n, pointers, rows, values)
    }
    pub fn nnz(&self) -> usize {
        self.nnz
    }
    #[inline]
    pub fn row(&self, row: usize) -> View<'_> {
        View { matrix: self, row }
    }
    #[inline]
    pub fn column(&self, col: usize) -> View<'_> {
        self.row(col)
    }
    /// The substitution kernel visits only pairs with column >= row.
    pub fn upper_row(&self, row: usize) -> Iter<'_> {
        let mut iter = self.row(row).iter();
        let entries = iter.entries.as_slice();
        let at = entries.partition_point(|link| (link.column as usize) < row);
        iter.entries = entries[at..].iter();
        iter
    }
    fn search(&self, row: usize, column: usize) -> Result<usize, usize> {
        let slot = self.slots[row];
        self.entries[slot.start..slot.start + slot.len as usize]
            .binary_search_by_key(&(column as u32), |link| link.column)
    }
    pub fn get(&self, row: usize, column: usize) -> f64 {
        if row == column {
            return self.diagonal(row);
        }
        self.search(row, column).map_or(0.0, |at| {
            self.values[self.entries[self.slots[row].start + at].coefficient as usize]
        })
    }
    pub fn set(&mut self, row: usize, column: usize, value: f64) {
        assert!(value.is_finite());
        if row == column {
            let old = self.diagonal(row);
            if old == value {
                return;
            }
            if self.diagonal.is_empty() {
                self.diagonal.resize(self.slots.len(), 0.0);
            }
            self.diagonal[row] = value;
            self.nnz += usize::from(old == 0.0);
            self.nnz -= usize::from(value == 0.0);
        } else {
            match self.search(row, column) {
                Ok(at) => {
                    let id = self.entries[self.slots[row].start + at].coefficient;
                    if self.values[id as usize] == value {
                        return;
                    }
                    if value == 0.0 {
                        self.remove_entry(row, at);
                        let other = self.search(column, row).unwrap();
                        self.remove_entry(column, other);
                        self.release(id);
                        self.nnz -= 2;
                    } else {
                        self.values[id as usize] = value;
                    }
                }
                Err(at) => {
                    if value == 0.0 {
                        return;
                    }
                    let id = if self.free != NONE {
                        let id = self.free;
                        self.free = self.values[id as usize].to_bits() as u32;
                        self.free_count -= 1;
                        self.values[id as usize] = value;
                        id
                    } else {
                        let id =
                            u32::try_from(self.values.len()).expect("Hessian exceeds u32 storage");
                        assert_ne!(id, NONE);
                        self.values.push(value);
                        id
                    };
                    self.insert_entry(
                        row,
                        at,
                        Link {
                            column: column as u32,
                            coefficient: id,
                        },
                    );
                    let other = self.search(column, row).unwrap_err();
                    self.insert_entry(
                        column,
                        other,
                        Link {
                            column: row as u32,
                            coefficient: id,
                        },
                    );
                    self.nnz += 2;
                }
            }
        }
        self.revision += 1;
        self.compact_values_if_fragmented();
    }
    fn release(&mut self, id: u32) {
        self.values[id as usize] = f64::from_bits(u64::from(self.free));
        self.free = id;
        self.free_count += 1;
    }
    /// A coefficient ID stays stable during adjacency relocation. When many
    /// coefficients have disappeared, rebuild their arena and remap only live
    /// adjacency slots; released slots may still contain obsolete IDs.
    #[inline]
    fn compact_values_if_fragmented(&mut self) {
        if self.values.len() < 1024 || 2 * self.free_count < self.values.len() {
            return;
        }
        self.compact_values();
    }
    #[cold]
    #[inline(never)]
    fn compact_values(&mut self) {
        let mut map = vec![NONE; self.values.len()];
        let mut values = Vec::with_capacity(self.values.len() - self.free_count);
        for slot in &self.slots {
            for link in &mut self.entries[slot.start..slot.start + slot.len as usize] {
                let mapped = &mut map[link.coefficient as usize];
                if *mapped == NONE {
                    *mapped = values.len() as u32;
                    values.push(self.values[link.coefficient as usize]);
                }
                link.coefficient = *mapped;
            }
        }
        self.values = values;
        self.free = NONE;
        self.free_count = 0;
    }
    fn remove_entry(&mut self, row: usize, at: usize) {
        let slot = self.slots[row];
        self.entries.copy_within(
            slot.start + at + 1..slot.start + slot.len as usize,
            slot.start + at,
        );
        self.slots[row].len -= 1;
    }
    fn insert_entry(&mut self, row: usize, at: usize, entry: Link) {
        let slot = self.slots[row];
        if slot.len < slot.capacity {
            self.entries.copy_within(
                slot.start + at..slot.start + slot.len as usize,
                slot.start + at + 1,
            );
            self.entries[slot.start + at] = entry;
            self.slots[row].len += 1;
            return;
        }
        let capacity = slot.capacity.saturating_mul(2).max(4);
        let start = self.entries.len();
        self.entries.reserve(capacity as usize);
        self.entries.extend_from_within(slot.start..slot.start + at);
        self.entries.push(entry);
        self.entries
            .extend_from_within(slot.start + at..slot.start + slot.len as usize);
        self.entries.resize(start + capacity as usize, EMPTY_LINK);
        self.slots[row] = Slot {
            start,
            len: slot.len + 1,
            capacity,
        };
        self.dead += slot.capacity as usize;
        self.compact_if_fragmented();
    }
    fn compact_if_fragmented(&mut self) {
        if self.entries.len() < 1024 || 2 * self.dead < self.entries.len() {
            return;
        }
        let mut entries = Vec::with_capacity(self.entries.len() - self.dead);
        for slot in &mut self.slots {
            let start = entries.len();
            entries.extend_from_slice(&self.entries[slot.start..slot.start + slot.len as usize]);
            *slot = Slot {
                start,
                len: slot.len,
                capacity: slot.len,
            };
        }
        self.entries = entries;
        self.dead = 0;
    }
    pub fn remove_variable(&mut self, column: usize) {
        self.revision += 1;
        if self.diagonal(column) != 0.0 {
            self.diagonal[column] = 0.0;
            self.nnz -= 1;
        }
        let slot = std::mem::take(&mut self.slots[column]);
        self.dead += slot.capacity as usize;
        for at in slot.start..slot.start + slot.len as usize {
            let link = self.entries[at];
            let row = link.column as usize;
            let other = self.search(row, column).unwrap();
            self.remove_entry(row, other);
            self.release(link.coefficient);
        }
        self.nnz -= 2 * slot.len as usize;
        self.compact_values_if_fragmented();
    }
}
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shared_coefficients_survive_recycling_compaction_and_packing() {
        let n = 80;
        let mut p = SymmetricMatrix::from_upper_columns(n, |j| {
            (0..=j).map(move |i| (i, (1 + i + j) as f64))
        });
        assert_eq!(p.values.len(), n * (n - 1) / 2);
        assert_eq!(p.entries.len(), n * (n - 1));
        for i in 0..60 {
            p.remove_variable(i);
        }
        assert!(p.values.len() < n * (n - 1) / 4);
        // Reuse stable variable IDs after coefficients and adjacency have moved.
        for j in 0..n {
            p.set(0, j, -(j as f64) - 1.0);
            p.set(j, 0, -(j as f64) - 2.0);
        }
        let active: Vec<_> = std::iter::once(0).chain(60..n).collect();
        let mut map = vec![usize::MAX; n];
        for (new, &old) in active.iter().enumerate() {
            map[old] = new;
        }
        // Remove the reintroduced intermediate variables before packing.
        for j in 1..60 {
            p.remove_variable(j);
        }
        let packed = p.pack_upper(&active, &map);
        let mut count = 0;
        for (new, &old) in active.iter().enumerate() {
            let full = p.row(old).to_vec();
            let mut cursor = p.upper_row(old);
            for j in old..n {
                assert_eq!(cursor.advance_to(j), p.get(old, j));
            }
            assert!(full.windows(2).all(|pair| pair[0].0 < pair[1].0));
            assert_eq!(p.row(old).iter().len(), full.len());
            assert_eq!(
                p.upper_row(old).collect::<Entries>(),
                full.iter()
                    .copied()
                    .filter(|&(j, _)| j >= old)
                    .collect::<Entries>()
            );
            let expected: Entries = full
                .iter()
                .filter_map(|&(i, v)| (map[i] <= new).then_some((map[i], v)))
                .collect();
            assert_eq!(packed.as_ref().column(new).collect::<Entries>(), expected);
            for &(j, v) in &full {
                assert_eq!(p.get(j, old), v);
                let expected = if j == 0 || old == 0 {
                    -(j.max(old) as f64) - 2.0
                } else {
                    (1 + j + old) as f64
                };
                assert_eq!(v, expected);
            }
            count += full.len();
        }
        assert_eq!(p.nnz(), count);
        // Updating an existing coefficient changes no adjacency or allocation.
        let ids = p.entries.clone();
        let coefficients = p.values.len();
        p.set(0, 79, 42.0);
        assert_eq!(p.get(79, 0), 42.0);
        assert_eq!(p.entries, ids);
        assert_eq!(p.values.len(), coefficients);
    }

    #[test]
    fn diagonal_only_storage_and_iterator_tails() {
        let mut p = SymmetricMatrix::zeros(4);
        assert!(p.diagonal.is_empty());
        for j in 0..4 {
            p.set(j, j, (j + 1) as f64);
        }
        assert!(p.values.is_empty() && p.entries.is_empty());
        p.set(1, 0, 3.0);
        p.set(1, 3, 5.0);
        let mut row = p.row(1).iter();
        assert_eq!(row.len(), 3);
        assert_eq!(row.next(), Some((0, 3.0)));
        assert_eq!(row.len(), 2);
        assert_eq!(row.clone().collect::<Entries>(), [(1, 2.0), (3, 5.0)]);
        assert_eq!(row.next(), Some((1, 2.0)));
        assert_eq!(row.next(), Some((3, 5.0)));
        assert_eq!(row.len(), 0);
        assert_eq!(row.next(), None);
        assert_eq!(p.add_variable(), 4);
        p.set(4, 4, 7.0);
        assert_eq!(p.row(4).to_vec(), [(4, 7.0)]);
    }

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
        assert_eq!(p.row(0).to_vec(), [(0, 2.0), (2, 3.0)]);
        assert_eq!(p.column(2).to_vec(), [(0, 3.0), (1, 6.0)]);
        assert_eq!(p.row(2).to_vec(), p.column(2).to_vec());
        let revision = p.revision;
        p.set(2, 0, 3.0);
        assert_eq!(p.revision, revision);
        assert_eq!(p.row(2).to_vec(), [(0, 3.0), (1, 6.0)]);
        p.remove_variable(2);
        assert_eq!(p.nnz(), 2);
        assert_eq!(p.row(0).to_vec(), [(0, 2.0)]);
        p.remove_variable(0);
        assert_eq!(p.nnz(), 1);
        assert_eq!(p.add_variable(), 3);
        p.set(0, 3, 7.0);
        assert_eq!(p.row(3).to_vec(), [(0, 7.0)]);
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
        assert_eq!(p.row(0).to_vec(), [(0, 2.0), (1, -1.0)]);
        assert_eq!(p.row(1).to_vec(), [(0, -1.0), (2, 3.0)]);
        assert_eq!(p.row(2).to_vec(), [(1, 3.0), (2, 4.0)]);
        assert_eq!(p.nnz(), 6);
        assert_eq!(p.revision, 0);
    }

    #[test]
    fn arena_growth_and_compaction_match_a_dense_reference() {
        let n = 40;
        let mut p = SymmetricMatrix::zeros(n);
        let mut dense = vec![vec![0.0; n]; n];
        let mut seed = 77u64;
        for step in 0..6000 {
            seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let i = (seed >> 33) as usize % n;
            let j = (seed >> 13) as usize % n;
            let value = ((seed >> 40) % 5) as f64 - 2.0;
            if step % 97 == 0 {
                p.remove_variable(i);
                dense[i].fill(0.0);
                for row in &mut dense {
                    row[i] = 0.0;
                }
            } else {
                p.set(i, j, value);
                dense[i][j] = value;
                dense[j][i] = value;
            }
            let mut count = 0;
            for (i, row) in dense.iter().enumerate() {
                let expected: Entries = row
                    .iter()
                    .enumerate()
                    .filter_map(|(j, &v)| (v != 0.0).then_some((j, v)))
                    .collect();
                assert_eq!(p.row(i).to_vec(), expected);
                count += expected.len();
            }
            assert_eq!(p.nnz(), count);
            // Relocation doubles a full slot and compaction reclaims dead
            // space once it covers half the arena, so growth stays bounded.
            assert!(p.entries.len() <= 8 * (count + n) + 2048);
        }
    }
}
