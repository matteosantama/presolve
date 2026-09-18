// SPDX-License-Identifier: Apache-2.0
// Modified for this library; copyright and attribution notices are in NOTICE.
//! Editable symmetric Hessian storage with one sorted adjacency slice per
//! variable, all held in a single arena. Indices remain stable until final
//! packing.

use crate::matrix::CscMatrix;

pub(crate) type Entries = Vec<(usize, f64)>;

/// Where a variable's sorted adjacency lives in the arena. A row that
/// outgrows its slot moves to the end of the arena; the old slot is dead
/// space until the next compaction.
#[derive(Clone, Copy, Debug, Default)]
struct Slot {
    start: usize,
    len: usize,
    capacity: usize,
}

/// One allocation for every row keeps construction to two counting passes
/// instead of one allocation per variable, and lets scans over consecutive
/// variables read consecutive memory.
#[derive(Clone, Debug)]
pub(crate) struct SymmetricMatrix {
    entries: Entries,
    slots: Vec<Slot>,
    /// Arena positions no live slot covers.
    dead: usize,
    nnz: usize,
    pub revision: usize,
}

impl SymmetricMatrix {
    pub fn zeros(n: usize) -> Self {
        Self {
            entries: Vec::new(),
            slots: vec![Slot::default(); n],
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
        let mut slots = vec![Slot::default(); n];
        for j in 0..n {
            for (i, v) in column(j) {
                if v != 0.0 {
                    slots[i].len += 1;
                    if i != j {
                        slots[j].len += 1;
                    }
                }
            }
        }
        let mut nnz = 0;
        for slot in &mut slots {
            slot.start = nnz;
            slot.capacity = slot.len;
            nnz += slot.len;
        }
        let mut entries = vec![(0, 0.0); nnz];
        let mut next: Vec<_> = slots.iter().map(|slot| slot.start).collect();
        for j in 0..n {
            for (i, v) in column(j) {
                if v != 0.0 {
                    entries[next[i]] = (j, v);
                    next[i] += 1;
                    if i != j {
                        entries[next[j]] = (i, v);
                        next[j] += 1;
                    }
                }
            }
        }
        // Lower entries arrive before the diagonal and upper entries, in order.
        Self {
            entries,
            slots,
            dead: 0,
            nnz,
            revision: 0,
        }
    }

    pub fn add_variable(&mut self) -> usize {
        self.revision += 1;
        let j = self.slots.len();
        self.slots.push(Slot::default());
        j
    }

    /// Upper triangle in CSC over compact columns. `compact_to_stable` lists
    /// the surviving stable columns in output order and `stable_to_compact`
    /// is its inverse, with `usize::MAX` for removed columns.
    pub fn pack_upper(
        &self,
        compact_to_stable: &[usize],
        stable_to_compact: &[usize],
    ) -> CscMatrix {
        let n = compact_to_stable.len();
        let nnz = self.nnz();
        let capacity = nnz.min(nnz / 2 + n);
        let mut pointers = Vec::with_capacity(n + 1);
        let mut rows = Vec::with_capacity(capacity);
        let mut values = Vec::with_capacity(capacity);
        pointers.push(0);
        for (j, &old) in compact_to_stable.iter().enumerate() {
            for &(i, v) in self.column(old) {
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
    pub fn row(&self, row: usize) -> &[(usize, f64)] {
        let slot = self.slots[row];
        &self.entries[slot.start..slot.start + slot.len]
    }
    pub fn column(&self, col: usize) -> &[(usize, f64)] {
        self.row(col)
    }

    pub fn get(&self, row: usize, column: usize) -> f64 {
        let entries = self.row(row);
        entries
            .binary_search_by_key(&column, |&(j, _)| j)
            .map_or(0.0, |index| entries[index].1)
    }

    pub fn set(&mut self, row: usize, column: usize, value: f64) {
        assert!(value.is_finite());
        let old = self.get(row, column);
        if old == value {
            return;
        }
        self.revision += 1;
        self.set_entry(row, column, value);
        if row != column {
            self.set_entry(column, row, value);
        }
        let count = if row == column { 1 } else { 2 };
        if old == 0.0 {
            self.nnz += count;
        }
        if value == 0.0 {
            self.nnz -= count;
        }
    }

    fn set_entry(&mut self, row: usize, index: usize, value: f64) {
        let slot = self.slots[row];
        let end = slot.start + slot.len;
        match self.entries[slot.start..end].binary_search_by_key(&index, |&(i, _)| i) {
            Ok(at) if value == 0.0 => {
                self.entries
                    .copy_within(slot.start + at + 1..end, slot.start + at);
                self.slots[row].len -= 1;
            }
            Ok(at) => self.entries[slot.start + at].1 = value,
            Err(at) if value != 0.0 => self.insert_entry(row, at, (index, value)),
            Err(_) => {}
        }
    }

    /// Insert into a slot with spare room in place; otherwise move the row to
    /// the end of the arena with room to grow, as a `Vec` would reallocate.
    fn insert_entry(&mut self, row: usize, at: usize, entry: (usize, f64)) {
        let slot = self.slots[row];
        if slot.len < slot.capacity {
            let end = slot.start + slot.len;
            self.entries
                .copy_within(slot.start + at..end, slot.start + at + 1);
            self.entries[slot.start + at] = entry;
            self.slots[row].len += 1;
            return;
        }
        let capacity = (2 * slot.capacity).max(4);
        let start = self.entries.len();
        self.entries.reserve(capacity);
        self.entries.extend_from_within(slot.start..slot.start + at);
        self.entries.push(entry);
        self.entries
            .extend_from_within(slot.start + at..slot.start + slot.len);
        self.entries.resize(start + capacity, (0, 0.0));
        self.slots[row] = Slot {
            start,
            len: slot.len + 1,
            capacity,
        };
        self.dead += slot.capacity;
        self.compact_if_fragmented();
    }

    /// Rebuild the arena once dead space dominates it. Row contents and their
    /// order are unchanged, so nothing observable depends on when this runs.
    fn compact_if_fragmented(&mut self) {
        const MIN_ARENA: usize = 1024;
        if self.entries.len() < MIN_ARENA || 2 * self.dead < self.entries.len() {
            return;
        }
        let mut entries = Vec::with_capacity(self.entries.len() - self.dead);
        for slot in &mut self.slots {
            let start = entries.len();
            entries.extend_from_slice(&self.entries[slot.start..slot.start + slot.len]);
            *slot = Slot {
                start,
                len: slot.len,
                capacity: slot.len,
            };
        }
        self.entries = entries;
        self.dead = 0;
    }

    /// Remove both symmetric views while keeping all other variable IDs stable.
    pub fn remove_variable(&mut self, column: usize) {
        self.revision += 1;
        let slot = std::mem::take(&mut self.slots[column]);
        self.dead += slot.capacity;
        // The released region is dead space that no other row's removal
        // touches (a removal only shifts inside that row's own slot), so it
        // can be read in place while the mirrored entries are cleared.
        for at in slot.start..slot.start + slot.len {
            let row = self.entries[at].0;
            if row != column {
                self.set_entry(row, column, 0.0);
                self.nnz -= 1;
            }
        }
        self.nnz -= slot.len;
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
        assert_eq!(p.row(2), [(0, 3.0), (1, 6.0)]);
        p.remove_variable(2);
        assert_eq!(p.nnz(), 2);
        assert_eq!(p.row(0), [(0, 2.0)]);
        p.remove_variable(0);
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
                assert_eq!(p.row(i), expected);
                count += expected.len();
            }
            assert_eq!(p.nnz(), count);
            // Relocation doubles a full slot and compaction reclaims dead
            // space once it covers half the arena, so growth stays bounded.
            assert!(p.entries.len() <= 8 * (count + n) + 2048);
        }
    }
}
