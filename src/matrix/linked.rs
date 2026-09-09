//! Experimental arena-backed orthogonal lists for the constraint matrix.
use crate::matrix::sparse::Entries;
const NONE: u32 = u32::MAX;

#[derive(Clone, Copy, Debug)]
struct Node {
    value: f64,
    row: u32,
    col: u32,
    prev: [u32; 2],
    next: [u32; 2],
}
#[derive(Clone, Copy, Debug)]
struct List {
    head: u32,
    tail: u32,
    len: usize,
}
impl Default for List {
    fn default() -> Self {
        Self {
            head: NONE,
            tail: NONE,
            len: 0,
        }
    }
}
#[derive(Clone, Debug)]
pub(crate) struct LinkedMatrix {
    nodes: Vec<Node>,
    lists: [Vec<List>; 2],
    free: u32,
    nnz: usize,
    column_cursors: Vec<u32>,
}
#[derive(Clone, Copy, Debug)]
pub(crate) struct View<'a> {
    matrix: &'a LinkedMatrix,
    list: List,
    axis: usize,
}
#[derive(Clone)]
pub(crate) struct Iter<'a> {
    matrix: &'a LinkedMatrix,
    cursor: Cursor,
    remaining: usize,
}
#[derive(Clone, Copy)]
pub(crate) struct Cursor {
    next: u32,
    axis: usize,
}
impl Cursor {
    pub fn next(&mut self, matrix: &LinkedMatrix) -> Option<(usize, f64)> {
        if self.next == NONE {
            return None;
        }
        let node = matrix.nodes[self.next as usize];
        self.next = node.next[self.axis];
        Some((
            if self.axis == 0 { node.col } else { node.row } as usize,
            node.value,
        ))
    }
}
impl Iterator for Iter<'_> {
    type Item = (usize, f64);
    fn next(&mut self) -> Option<Self::Item> {
        let entry = self.cursor.next(self.matrix)?;
        self.remaining -= 1;
        Some(entry)
    }
    fn size_hint(&self) -> (usize, Option<usize>) {
        (self.remaining, Some(self.remaining))
    }
}
impl ExactSizeIterator for Iter<'_> {}
impl<'a> View<'a> {
    pub fn len(self) -> usize {
        self.list.len
    }
    pub fn is_empty(self) -> bool {
        self.len() == 0
    }
    pub fn cursor(self) -> Cursor {
        Cursor {
            next: self.list.head,
            axis: self.axis,
        }
    }
    pub fn iter(self) -> Iter<'a> {
        Iter {
            matrix: self.matrix,
            cursor: self.cursor(),
            remaining: self.len(),
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
#[cfg(test)]
impl PartialEq for View<'_> {
    fn eq(&self, rhs: &Self) -> bool {
        self.iter().eq(rhs.iter())
    }
}
#[cfg(test)]
impl<const N: usize> PartialEq<[(usize, f64); N]> for View<'_> {
    fn eq(&self, rhs: &[(usize, f64); N]) -> bool {
        self.iter().eq(rhs.iter().copied())
    }
}
impl LinkedMatrix {
    pub fn zeros(rows: usize, columns: usize) -> Self {
        assert!(rows < NONE as usize && columns < NONE as usize);
        Self {
            nodes: Vec::new(),
            lists: [vec![List::default(); rows], vec![List::default(); columns]],
            free: NONE,
            nnz: 0,
            column_cursors: vec![NONE; columns],
        }
    }
    pub fn from_columns<I: Iterator<Item = (usize, f64)>>(
        rows: usize,
        columns: usize,
        column: impl Fn(usize) -> I,
    ) -> Self {
        let visit = |f: &mut dyn FnMut(usize, usize, f64)| {
            for j in 0..columns {
                for (i, v) in column(j) {
                    f(i, j, v);
                }
            }
        };
        let mut out = Self::zeros(rows, columns);
        // Presolve repeatedly scans complete rows. Lay the arena out in row
        // order while retaining sorted links in both dimensions.
        visit(&mut |row, _, value| {
            if value != 0.0 {
                out.lists[0][row].len += 1;
            }
        });
        let mut offsets = Vec::with_capacity(rows);
        let mut nnz = 0;
        for list in &mut out.lists[0] {
            offsets.push(nnz);
            if list.len != 0 {
                list.head = u32::try_from(nnz).unwrap();
                list.tail = u32::try_from(nnz + list.len - 1).unwrap();
            }
            nnz += list.len;
        }
        assert!(nnz < NONE as usize);
        out.nodes = vec![
            Node {
                value: 0.0,
                row: 0,
                col: 0,
                prev: [NONE; 2],
                next: [NONE; 2]
            };
            nnz
        ];
        visit(&mut |row, col, value| {
            if value == 0.0 {
                return;
            }
            let id = offsets[row] as u32;
            offsets[row] += 1;
            let row_list = out.lists[0][row];
            let col_list = &mut out.lists[1][col];
            out.nodes[id as usize] = Node {
                value,
                row: row as u32,
                col: col as u32,
                prev: [
                    if id == row_list.head { NONE } else { id - 1 },
                    col_list.tail,
                ],
                next: [if id == row_list.tail { NONE } else { id + 1 }, NONE],
            };
            if col_list.tail == NONE {
                col_list.head = id;
            } else {
                out.nodes[col_list.tail as usize].next[1] = id;
            }
            col_list.tail = id;
            col_list.len += 1;
            out.column_cursors[col] = id;
        });
        out.nnz = nnz;
        out
    }
    pub fn add_column(&mut self) -> usize {
        let i = self.lists[1].len();
        self.lists[1].push(List::default());
        self.column_cursors.push(NONE);
        i
    }
    pub fn nnz(&self) -> usize {
        self.nnz
    }
    pub fn row(&self, row: usize) -> View<'_> {
        View {
            matrix: self,
            list: self.lists[0][row],
            axis: 0,
        }
    }
    pub fn column(&self, col: usize) -> View<'_> {
        View {
            matrix: self,
            list: self.lists[1][col],
            axis: 1,
        }
    }
    fn key(&self, node: u32, axis: usize) -> usize {
        let n = &self.nodes[node as usize];
        if axis == 0 {
            n.col as usize
        } else {
            n.row as usize
        }
    }
    fn position(&self, axis: usize, index: usize, key: usize) -> u32 {
        let list = self.lists[axis][index];
        if list.tail == NONE || self.key(list.tail, axis) < key {
            return NONE;
        }
        let mut at = list.head;
        while at != NONE && self.key(at, axis) < key {
            at = self.nodes[at as usize].next[axis];
        }
        at
    }
    // Reuse a live node near the previous insertion. Batched rows are visited
    // in increasing order, so each column cursor advances instead of restarting.
    fn column_position(&self, col: usize, row: usize) -> u32 {
        let list = self.lists[1][col];
        if list.tail == NONE || self.key(list.tail, 1) < row {
            return NONE;
        }
        if self.key(list.head, 1) >= row {
            return list.head;
        }
        let mut at = self.column_cursors[col];
        if at == NONE {
            at = list.head;
        }
        // Prefer an endpoint when it is closer in row-index space.
        for candidate in [list.head, list.tail] {
            if self.key(candidate, 1).abs_diff(row) < self.key(at, 1).abs_diff(row) {
                at = candidate;
            }
        }
        if self.key(at, 1) < row {
            while at != NONE && self.key(at, 1) < row {
                at = self.nodes[at as usize].next[1];
            }
        } else {
            while self.nodes[at as usize].prev[1] != NONE {
                let prev = self.nodes[at as usize].prev[1];
                if self.key(prev, 1) < row {
                    break;
                }
                at = prev;
            }
        }
        at
    }
    pub fn get(&self, row: usize, col: usize) -> f64 {
        let (axis, index, key) = if self.lists[0][row].len <= self.lists[1][col].len {
            (0, row, col)
        } else {
            (1, col, row)
        };
        let at = self.position(axis, index, key);
        if at != NONE && self.key(at, axis) == key {
            self.nodes[at as usize].value
        } else {
            0.0
        }
    }
    fn insert(&mut self, row: usize, col: usize, value: f64, next: [u32; 2]) {
        let indices = [row, col];
        let prev = std::array::from_fn(|axis| {
            if next[axis] == NONE {
                self.lists[axis][indices[axis]].tail
            } else {
                self.nodes[next[axis] as usize].prev[axis]
            }
        });
        let node = Node {
            value,
            row: row.try_into().unwrap(),
            col: col.try_into().unwrap(),
            prev,
            next,
        };
        let id = if self.free == NONE {
            let id = u32::try_from(self.nodes.len()).unwrap();
            assert_ne!(id, NONE);
            self.nodes.push(node);
            id
        } else {
            let id = self.free;
            self.free = self.nodes[id as usize].next[0];
            self.nodes[id as usize] = node;
            id
        };
        for axis in 0..2 {
            let list = &mut self.lists[axis][indices[axis]];
            if prev[axis] == NONE {
                list.head = id;
            } else {
                self.nodes[prev[axis] as usize].next[axis] = id;
            }
            if next[axis] == NONE {
                list.tail = id;
            } else {
                self.nodes[next[axis] as usize].prev[axis] = id;
            }
            list.len += 1;
        }
        self.column_cursors[col] = id;
        self.nnz += 1;
    }
    fn remove(&mut self, id: u32) {
        let node = self.nodes[id as usize];
        if self.column_cursors[node.col as usize] == id {
            self.column_cursors[node.col as usize] = if node.next[1] != NONE {
                node.next[1]
            } else {
                node.prev[1]
            };
        }
        for (axis, index) in [node.row as usize, node.col as usize]
            .into_iter()
            .enumerate()
        {
            let list = &mut self.lists[axis][index];
            if node.prev[axis] == NONE {
                list.head = node.next[axis];
            } else {
                self.nodes[node.prev[axis] as usize].next[axis] = node.next[axis];
            }
            if node.next[axis] == NONE {
                list.tail = node.prev[axis];
            } else {
                self.nodes[node.next[axis] as usize].prev[axis] = node.prev[axis];
            }
            list.len -= 1;
        }
        self.nodes[id as usize].next[0] = self.free;
        self.free = id;
        self.nnz -= 1;
    }
    #[cfg(test)]
    pub fn set(&mut self, row: usize, col: usize, value: f64) {
        assert!(value.is_finite());
        let at = self.position(0, row, col);
        if at != NONE && self.key(at, 0) == col {
            if value == 0.0 {
                self.remove(at);
            } else {
                self.nodes[at as usize].value = value;
            }
        } else if value != 0.0 {
            self.insert(row, col, value, [at, self.column_position(col, row)]);
        }
    }
    #[cfg(test)]
    pub fn remove_row(&mut self, row: usize) -> Entries {
        let mut entries = Vec::with_capacity(self.lists[0][row].len);
        while self.lists[0][row].head != NONE {
            let id = self.lists[0][row].head;
            let node = self.nodes[id as usize];
            entries.push((node.col as usize, node.value));
            self.remove(id);
        }
        entries
    }
    pub fn remove_column(&mut self, col: usize) -> Entries {
        let mut entries = Vec::with_capacity(self.lists[1][col].len);
        while self.lists[1][col].head != NONE {
            let id = self.lists[1][col].head;
            let node = self.nodes[id as usize];
            entries.push((node.row as usize, node.value));
            self.remove(id);
        }
        entries
    }
    /// Keep existing nodes and their column links while merging a sorted row.
    /// Save the old coefficients only when the caller needs them for lock updates.
    fn update_row<const SAVE: bool, const TRACK: bool>(
        &mut self,
        row: usize,
        entries: &[(usize, f64)],
        touched: &mut Vec<usize>,
    ) -> Entries {
        let mut old = if SAVE {
            Vec::with_capacity(self.lists[0][row].len)
        } else {
            Vec::new()
        };
        let mut at = self.lists[0][row].head;
        for &(col, value) in entries {
            debug_assert!(value.is_finite() && value != 0.0);
            while at != NONE && self.key(at, 0) < col {
                let node = self.nodes[at as usize];
                if SAVE {
                    old.push((node.col as usize, node.value));
                }
                if TRACK {
                    touched.push(node.col as usize);
                }
                self.remove(at);
                at = node.next[0];
            }
            if at != NONE && self.key(at, 0) == col {
                let node = self.nodes[at as usize];
                if SAVE {
                    old.push((col, node.value));
                }
                if node.value != value {
                    if TRACK {
                        touched.push(col);
                    }
                    self.nodes[at as usize].value = value;
                }
                at = node.next[0];
            } else {
                if TRACK {
                    touched.push(col);
                }
                self.insert(row, col, value, [at, self.column_position(col, row)]);
            }
        }
        while at != NONE {
            let node = self.nodes[at as usize];
            if SAVE {
                old.push((node.col as usize, node.value));
            }
            if TRACK {
                touched.push(node.col as usize);
            }
            self.remove(at);
            at = node.next[0];
        }
        old
    }

    pub fn replace_row(&mut self, row: usize, entries: &[(usize, f64)]) -> Entries {
        self.update_row::<true, false>(row, entries, &mut Vec::new())
    }

    pub fn replace_rows(&mut self, mut updates: Vec<(usize, Entries)>) -> Vec<usize> {
        updates.sort_unstable_by_key(|(row, _)| *row);
        let mut touched = Vec::new();
        for (row, entries) in updates {
            self.update_row::<false, true>(row, &entries, &mut touched);
        }
        touched.sort_unstable();
        touched.dedup();
        touched
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn row_replacement_retains_live_nodes_and_returns_original_coefficients() {
        let mut a = LinkedMatrix::zeros(3, 5);
        for i in 0..3 {
            for j in [0, 2, 4] {
                a.set(i, j, 1.0);
            }
        }
        let keep = a.position(0, 1, 2);
        assert_eq!(
            a.replace_row(1, &[(1, 3.0), (2, 4.0), (3, 5.0)]),
            [(0, 1.0), (2, 1.0), (4, 1.0)]
        );
        assert_eq!(a.position(0, 1, 2), keep);
        assert_eq!(a.column(2).to_vec(), [(0, 1.0), (1, 4.0), (2, 1.0)]);
        assert_eq!(a.remove_row(1), [(1, 3.0), (2, 4.0), (3, 5.0)]);
        a.replace_row(1, &[(0, 7.0), (2, 8.0), (4, 9.0)]);
        assert_eq!(a.column(2).to_vec(), [(0, 1.0), (1, 8.0), (2, 1.0)]);
    }

    #[test]
    fn linked_mutations_match_dense_storage_and_recycle_nodes() {
        assert_eq!(std::mem::size_of::<Node>(), 32);
        let mut a = LinkedMatrix::zeros(13, 17);
        let mut dense = [[0.0; 17]; 13];
        let mut seed = 3123u64;
        for step in 0..1000 {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
            let row = (seed >> 32) as usize % 13;
            let col = (seed >> 16) as usize % 17;
            let value = (seed % 7) as f64 - 3.0;
            match step % 11 {
                0 => {
                    a.remove_row(row);
                    dense[row].fill(0.0);
                }
                1 => {
                    a.remove_column(col);
                    for row in &mut dense {
                        row[col] = 0.0;
                    }
                }
                2 => {
                    let next = (row + 1) % 13;
                    let entries = vec![(col, 2.0)];
                    let touched = a.replace_rows(vec![(row, entries.clone()), (next, entries)]);
                    assert!(touched.windows(2).all(|w| w[0] < w[1]));
                    dense[row].fill(0.0);
                    dense[row][col] = 2.0;
                    dense[next].fill(0.0);
                    dense[next][col] = 2.0;
                }
                3 => {
                    let old = a.row(row).to_vec();
                    assert_eq!(a.replace_row(row, &[(col, 2.0)]), old);
                    dense[row].fill(0.0);
                    dense[row][col] = 2.0;
                }
                _ => {
                    a.set(row, col, value);
                    dense[row][col] = value;
                }
            }
            let mut count = 0;
            for (i, row) in dense.iter().enumerate() {
                let expected: Entries = row
                    .iter()
                    .enumerate()
                    .filter_map(|(j, &v)| (v != 0.0).then_some((j, v)))
                    .collect();
                assert_eq!(a.row(i).to_vec(), expected);
                count += expected.len();
                for (j, &v) in row.iter().enumerate() {
                    assert_eq!(a.get(i, j), v);
                }
            }
            for j in 0..17 {
                let expected: Entries = dense
                    .iter()
                    .enumerate()
                    .filter_map(|(i, row)| (row[j] != 0.0).then_some((i, row[j])))
                    .collect();
                assert_eq!(a.column(j).to_vec(), expected);
                let cursor = a.column_cursors[j];
                if cursor != NONE {
                    let node = a.nodes[cursor as usize];
                    assert_eq!(node.col as usize, j);
                    assert_eq!(a.position(1, j, node.row as usize), cursor);
                }
            }
            assert_eq!(a.nnz(), count);
            assert!(a.nodes.len() <= 13 * 17);
            for axis in 0..2 {
                for list in &a.lists[axis] {
                    let mut at = list.head;
                    let mut previous = NONE;
                    let mut visited = 0;
                    while at != NONE {
                        assert_eq!(a.nodes[at as usize].prev[axis], previous);
                        previous = at;
                        at = a.nodes[at as usize].next[axis];
                        visited += 1;
                        assert!(visited <= a.nnz());
                    }
                    assert_eq!(previous, list.tail);
                    assert_eq!(visited, list.len);
                }
            }
        }
    }
}
