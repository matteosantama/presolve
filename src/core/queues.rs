// SPDX-License-Identifier: Apache-2.0
// Modified for this library; copyright and attribution notices are in NOTICE.
//! Deduplicated work queues keyed by stable problem indices.
//! A drained batch is distinct from the next scheduling round.

#[derive(Clone, Debug)]
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
        let entries = std::mem::take(&mut self.entries);
        for &index in &entries {
            self.queued[index] = false;
        }
        entries
    }

    pub fn iter(&self) -> impl Iterator<Item = usize> + '_ {
        self.entries.iter().copied()
    }

    #[cfg(test)]
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
}

pub(crate) struct Queues {
    pub empty_rows: Worklist,
    pub singleton_rows: Worklist,
    pub doubleton_rows: Worklist,
    pub short_equalities: Worklist,
    pub empty_columns: Worklist,
    pub singleton_columns: Worklist,
    pub changed_activities: Worklist,
    pub fixed_columns: Worklist,
    pub unlocked_columns: Worklist,
    pub singleton_activity_rows: Worklist,
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
            changed_activities: Worklist::new(rows),
            fixed_columns: Worklist::new(columns),
            unlocked_columns: Worklist::new(columns),
            singleton_activity_rows: Worklist::new(rows),
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
        self.singleton_activity_rows.push(row);
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
}
