// SPDX-License-Identifier: Apache-2.0
// Modified for this library; copyright and attribution notices are in NOTICE.
//! Deduplicated work queues keyed by stable problem indices.
//! A drained batch is distinct from the next scheduling round.

#[derive(Clone, Debug, Default)]
pub(crate) struct Worklist {
    entries: Vec<usize>,
    queued: Vec<bool>,
}

impl Worklist {
    pub fn new(size: usize) -> Self {
        Self {
            entries: Vec::new(),
            queued: vec![false; size],
        }
    }

    pub fn push(&mut self, index: usize) {
        if index >= self.queued.len() {
            self.queued.resize(index + 1, false);
        }
        if !std::mem::replace(&mut self.queued[index], true) {
            self.entries.push(index);
        }
    }

    pub fn pop(&mut self) -> Option<usize> {
        let index = self.entries.pop()?;
        self.queued[index] = false;
        Some(index)
    }

    pub fn take_round(&mut self) -> Vec<usize> {
        // Rounds tend to repeat in size, so size the next buffer up front
        // instead of regrowing it from empty through every doubling step.
        let next = Vec::with_capacity(self.entries.len());
        let entries = std::mem::replace(&mut self.entries, next);
        for &index in &entries {
            self.queued[index] = false;
        }
        entries
    }

    /// Move the current round into `round`, whose previous allocation then
    /// receives the next round. Round-based rules stop allocating per round.
    pub fn swap_round(&mut self, round: &mut Vec<usize>) {
        round.clear();
        std::mem::swap(&mut self.entries, round);
        for &index in round.iter() {
            self.queued[index] = false;
        }
    }

    /// Reuse the allocation for a different index space, discarding entries.
    pub fn reset(&mut self, size: usize) {
        self.clear();
        self.queued.clear();
        self.queued.resize(size, false);
    }

    pub fn clear(&mut self) {
        for &index in &self.entries {
            self.queued[index] = false;
        }
        self.entries.clear();
    }

    #[cfg(test)]
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
}

/// Rows whose activity changed, queued for bound propagation and for the
/// singleton-column scan. Every change is pushed to both consumers, so one
/// flag byte per row serves both lists: a push over a column touches one
/// scattered location per incident row instead of two. Each list keeps its
/// own first-push order and is drained independently.
#[derive(Clone, Debug, Default)]
pub(crate) struct ActivityRows {
    flags: Vec<u8>,
    propagation: Vec<usize>,
    singleton: Vec<usize>,
}

const PROPAGATION: u8 = 1;
const SINGLETON: u8 = 2;

impl ActivityRows {
    pub fn new(rows: usize) -> Self {
        Self {
            flags: vec![0; rows],
            propagation: Vec::new(),
            singleton: Vec::new(),
        }
    }

    pub fn push(&mut self, row: usize) {
        let flags = self.flags[row];
        if flags & PROPAGATION == 0 {
            self.propagation.push(row);
        }
        if flags & SINGLETON == 0 {
            self.singleton.push(row);
        }
        self.flags[row] = PROPAGATION | SINGLETON;
    }

    fn drain(entries: &mut Vec<usize>, flags: &mut [u8], bit: u8) -> Vec<usize> {
        let next = Vec::with_capacity(entries.len());
        let entries = std::mem::replace(entries, next);
        for &row in &entries {
            flags[row] &= !bit;
        }
        entries
    }

    /// Rows pending for bound propagation.
    pub fn take_round(&mut self) -> Vec<usize> {
        Self::drain(&mut self.propagation, &mut self.flags, PROPAGATION)
    }

    /// Rows pending for the singleton-column scan.
    pub fn take_singleton_round(&mut self) -> Vec<usize> {
        Self::drain(&mut self.singleton, &mut self.flags, SINGLETON)
    }

    pub fn iter(&self) -> impl Iterator<Item = usize> + '_ {
        self.propagation.iter().copied()
    }
}

pub(crate) struct Queues {
    pub empty_rows: Worklist,
    pub singleton_rows: Worklist,
    pub doubleton_rows: Worklist,
    pub short_equalities: Worklist,
    pub empty_columns: Worklist,
    pub singleton_columns: Worklist,
    pub changed_activities: ActivityRows,
    pub fixed_columns: Worklist,
    pub unlocked_columns: Worklist,
}

impl Queues {
    pub fn new(rows: usize, columns: usize) -> Self {
        Self {
            empty_rows: Worklist::new(rows),
            singleton_rows: Worklist::new(rows),
            doubleton_rows: Worklist::new(rows),
            short_equalities: Worklist::new(rows),
            empty_columns: Worklist::new(columns),
            singleton_columns: Worklist::new(columns),
            changed_activities: ActivityRows::new(rows),
            fixed_columns: Worklist::new(columns),
            unlocked_columns: Worklist::new(columns),
        }
    }

    pub fn row_changed(&mut self, row: usize, size: usize, equality: bool) {
        match size {
            0 => self.empty_rows.push(row),
            1 => self.singleton_rows.push(row),
            2 if equality => self.doubleton_rows.push(row),
            3.. if equality => self.short_equalities.push(row),
            _ => {}
        }
        self.changed_activities.push(row);
    }

    pub fn column_changed(&mut self, column: usize, size: usize) {
        self.unlocked_columns.push(column);
        match size {
            0 => self.empty_columns.push(column),
            1 => self.singleton_columns.push(column),
            _ => {}
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn changes_during_propagation_are_scheduled_for_the_next_round() {
        let mut queue = Worklist::new(3);
        queue.push(1);
        queue.push(1);
        queue.push(2);
        let round = queue.take_round();
        queue.push(1);
        queue.push(1);
        assert_eq!(round, [1, 2]);
        assert_eq!(queue.pop(), Some(1));
        assert!(queue.is_empty());
    }

    #[test]
    fn swapped_rounds_reuse_the_previous_round_and_unmark_entries() {
        let mut queue = Worklist::new(3);
        let mut round = vec![7, 8, 9];
        queue.push(2);
        queue.push(0);
        queue.swap_round(&mut round);
        assert_eq!(round, [2, 0]);
        assert!(queue.is_empty());
        queue.push(2);
        queue.swap_round(&mut round);
        assert_eq!(round, [2]);
        queue.swap_round(&mut round);
        assert!(round.is_empty());
    }

    #[test]
    fn shared_activity_flags_keep_independent_rounds_for_both_consumers() {
        let mut rows = ActivityRows::new(4);
        rows.push(2);
        rows.push(0);
        rows.push(2);
        assert_eq!(rows.iter().collect::<Vec<_>>(), [2, 0]);
        assert_eq!(rows.take_round(), [2, 0]);
        // The singleton consumer still holds both rows; a repeated push
        // reaches only the drained propagation list.
        rows.push(0);
        rows.push(3);
        assert_eq!(rows.take_singleton_round(), [2, 0, 3]);
        assert_eq!(rows.take_round(), [0, 3]);
        assert!(rows.take_round().is_empty());
        assert!(rows.take_singleton_round().is_empty());
    }
}
