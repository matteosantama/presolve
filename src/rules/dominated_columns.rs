// SPDX-License-Identifier: Apache-2.0
// Modified for this library; copyright and attribution notices are in NOTICE.
//! Fix a column dominated by another. Moving along `(Δx_j, Δx_k) = (δ, -δ)`
//! keeps every row feasible and never increases a linear objective, so a
//! feasible point can slide until the dominated variable reaches its bound.

use crate::result::RuleId;
use crate::{
    model::tape::{Certificate, Side},
    model::{Model, RowDomain},
    problem::Bounds,
};

/// Candidates come from a column's shortest row. Nearly every dominated pair
/// in the benchmark corpora shares a row this short, and longer rows would
/// spend the allowance on visits that rarely pay off.
const MAX_CANDIDATE_ROW: usize = 256;

/// Members of one support-hash run are compared with this many followers,
/// which keeps a large run of identical supports linear instead of quadratic.
const MAX_GROUP_NEIGHBOURS: usize = 32;

const UPPER_OPEN: u8 = 1;
const LOWER_OPEN: u8 = 2;
const KNOWN: u8 = 4;
const UNKNOWN_FINGERPRINT: u64 = u64::MAX;

/// Per-pass column data, computed on first use so a pass touches only the
/// columns it visits: equality-row fingerprint and which bound sides are
/// absent or implied by a row. A candidate pair is then screened in O(1).
#[derive(Debug, Default)]
pub(crate) struct DominatedScratch {
    fingerprints: Vec<u64>,
    open: Vec<u8>,
}

impl DominatedScratch {
    fn fingerprint(&mut self, model: &Model, j: usize) -> u64 {
        if self.fingerprints[j] == UNKNOWN_FINGERPRINT {
            self.fingerprints[j] = model.equality_fingerprint(j);
        }
        self.fingerprints[j]
    }

    fn open(&mut self, model: &mut Model, j: usize) -> u8 {
        if self.open[j] & KNOWN == 0 {
            self.open[j] = model.open_sides(j) | KNOWN;
        }
        self.open[j]
    }
}

impl Model {
    /// Entries in equality and ranged rows must agree exactly between the two
    /// columns of a dominated pair, so their hash is a cheap necessary condition.
    fn equality_fingerprint(&self, j: usize) -> u64 {
        let mut hash: u64 = 0xcbf2_9ce4_8422_2325;
        for (i, a) in self.a.column(j) {
            if matches!(self.rows[i], RowDomain::Linear(b) if b.lower.is_finite() && b.upper.is_finite())
            {
                for word in [i as u64, a.to_bits()] {
                    hash ^= word;
                    hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
                }
            }
        }
        hash
    }

    /// A side is open when it is infinite or a row already implies it: the
    /// shifted point satisfies every row and every other bound, hence that
    /// side too. Cached activities only; a cancellation-prone row is skipped
    /// rather than recomputed, which is conservative.
    fn open_sides(&mut self, j: usize) -> u8 {
        let mut open = 0;
        for (side, flag) in [(Side::Upper, UPPER_OPEN), (Side::Lower, LOWER_OPEN)] {
            if !side.value(self.bounds[j]).is_finite() || self.bound_implied(j, side, false) {
                open |= flag;
            }
        }
        open
    }

    /// Whether `c_j <= c_k` and every row allows shifting weight from `k` to `j`.
    /// Returns the merge steps performed, since most pairs fail early.
    fn dominates(&self, j: usize, k: usize) -> (bool, usize) {
        if self.objective.c[j] > self.objective.c[k] {
            return (false, 0);
        }
        let mut steps = 0;
        let mut left = self.a.column(j).iter().peekable();
        let mut right = self.a.column(k).iter().peekable();
        while left.peek().is_some() || right.peek().is_some() {
            steps += 1;
            let (i, d) = match (left.peek().copied(), right.peek().copied()) {
                (Some((i, a)), Some((l, b))) if i == l => {
                    left.next();
                    right.next();
                    (i, a - b)
                }
                (Some((i, a)), Some((l, _))) if i < l => {
                    left.next();
                    (i, a)
                }
                (Some((i, a)), None) => {
                    left.next();
                    (i, a)
                }
                (_, Some((l, b))) => {
                    right.next();
                    (l, -b)
                }
                (None, None) => unreachable!(),
            };
            let RowDomain::Linear(bounds) = self.rows[i] else {
                return (false, steps);
            };
            let allowed = match (bounds.lower.is_finite(), bounds.upper.is_finite()) {
                (true, true) => d == 0.0,
                (false, true) => d <= 0.0,
                (true, false) => d >= 0.0,
                (false, false) => true,
            };
            if !allowed || !d.is_finite() {
                return (false, steps);
            }
        }
        (true, steps)
    }

    /// Whether a proved `dominant` over `dominated` relation could fix a column.
    fn shift_target(
        &mut self,
        scratch: &mut DominatedScratch,
        dominant: usize,
        dominated: usize,
    ) -> bool {
        let up = self.bounds[dominant].upper;
        let down = self.bounds[dominated].lower;
        (down.is_finite() && scratch.open(self, dominant) & UPPER_OPEN != 0)
            || (up.is_finite() && scratch.open(self, dominated) & LOWER_OPEN != 0)
            || (!up.is_finite() && !down.is_finite())
    }

    /// Examine one unordered pair. Returns the fixed column, if any, and
    /// charges the merge work performed. Both columns must be alive and linear.
    fn dominated_pair(
        &mut self,
        scratch: &mut DominatedScratch,
        j: usize,
        k: usize,
        work: &mut usize,
    ) -> Result<Option<usize>, Certificate> {
        if scratch.fingerprint(self, j) != scratch.fingerprint(self, k) {
            return Ok(None);
        }
        // Screen on actual bounds first: a direction whose sides are both
        // finite can only fix a column through an implied side, so compute the
        // open sides before paying for a merge on such pairs.
        let bj = self.bounds[j];
        let bk = self.bounds[k];
        let quick = |dominant: Bounds, dominated: Bounds| {
            !dominant.upper.is_finite() || !dominated.lower.is_finite()
        };
        if !quick(bj, bk)
            && !quick(bk, bj)
            && !self.shift_target(scratch, j, k)
            && !self.shift_target(scratch, k, j)
        {
            return Ok(None);
        }
        let (forward, steps) = self.dominates(j, k);
        let (backward, more) = if forward {
            (false, 0)
        } else {
            self.dominates(k, j)
        };
        *work = work.saturating_sub(steps + more);
        let (dominant, dominated) = if forward {
            (j, k)
        } else if backward {
            (k, j)
        } else {
            return Ok(None);
        };
        if !self.shift_target(scratch, dominant, dominated) {
            return Ok(None);
        }
        let up = self.bounds[dominant].upper;
        let down = self.bounds[dominated].lower;
        let (column, value) = if down.is_finite() && scratch.open(self, dominant) & UPPER_OPEN != 0
        {
            (dominated, down)
        } else if up.is_finite() && scratch.open(self, dominated) & LOWER_OPEN != 0 {
            (dominant, up)
        } else {
            if self.objective.c[dominant] == self.objective.c[dominated] {
                return Ok(None);
            }
            return Err(self.dual_certificate([(dominant, 1.0), (dominated, -1.0)]));
        };
        if self.fix(column, value) {
            scratch.fingerprints[column] = UNKNOWN_FINGERPRINT;
            Ok(Some(column))
        } else {
            Ok(None)
        }
    }

    fn begin_pass(&mut self) -> (DominatedScratch, usize) {
        let n = self.bounds.len();
        let mut scratch = std::mem::take(&mut self.dominated_scratch);
        scratch.fingerprints.clear();
        scratch.fingerprints.resize(n, UNKNOWN_FINGERPRINT);
        scratch.open.clear();
        scratch.open.resize(n, 0);
        let work = self
            .settings
            .dominated_columns
            .work_limit
            .resolve(self.a.nnz().saturating_mul(2));
        (scratch, work)
    }

    /// Test the pairs inside each run of columns sharing a constraint support
    /// hash. The parallel-column scan already sorted these, so this costs only
    /// the merges within collision runs. Returns the number of comparisons.
    pub(super) fn dominated_support_groups(
        &mut self,
        entries: &[(super::parallel::ColumnKey, usize)],
    ) -> Result<usize, Certificate> {
        let (mut scratch, mut work) = self.begin_pass();
        let mut comparisons = 0;
        let mut start = 0;
        let result = (|| {
            while start < entries.len() {
                let support = super::parallel::support_hash(entries[start].0.0);
                let mut end = start + 1;
                while end < entries.len()
                    && super::parallel::support_hash(entries[end].0.0) == support
                {
                    end += 1;
                }
                let run = &entries[start..end];
                start = end;
                if run.len() < 2
                    || run
                        .iter()
                        .any(|&(_, j)| !self.objective.p.column(j).is_empty())
                {
                    continue;
                }
                for (at, &(_, j)) in run.iter().enumerate() {
                    for &(_, k) in run[at + 1..].iter().take(MAX_GROUP_NEIGHBOURS) {
                        if work == 0 {
                            return Ok(());
                        }
                        if !self.linear_column(j) {
                            break;
                        }
                        if !self.linear_column(k) {
                            continue;
                        }
                        comparisons += 1;
                        work -= 1;
                        if self.dominated_pair(&mut scratch, j, k, &mut work)? == Some(j) {
                            break;
                        }
                    }
                }
            }
            Ok(())
        })();
        self.dominated_scratch = scratch;
        result.map(|()| comparisons)
    }

    /// General search: candidates for `k` are the columns of its shortest row.
    /// The fixed variable's recovered reduced cost has the sign of the
    /// dominating column's, so no dual transformation is recorded beyond `Fixed`.
    pub fn dominated_columns(&mut self) -> Result<usize, Certificate> {
        self.enter(RuleId::DominatedColumns);
        let (mut scratch, mut work) = self.begin_pass();
        let n = self.bounds.len();
        let mut fixed = 0;
        let result = (|| {
            for k in 0..n {
                if work == 0 {
                    break;
                }
                if !self.linear_column(k) || self.bounds[k].equality() {
                    continue;
                }
                let Some((row, _)) = self
                    .a
                    .column(k)
                    .iter()
                    .min_by_key(|&(i, _)| (self.a.row(i).len(), i))
                else {
                    continue;
                };
                if self.a.row(row).len() > MAX_CANDIDATE_ROW {
                    continue;
                }
                let mut cursor = self.a.row(row).cursor();
                while let Some((j, _)) = cursor.next(&self.a) {
                    // Every visit is charged so a pass stays proportional to nnz.
                    if work == 0 {
                        break;
                    }
                    work -= 1;
                    // A pair with nested supports is found from its smaller
                    // column, whose every row contains the other. Skip the
                    // mirror visit.
                    if j == k
                        || !self.alive[k]
                        || (self.a.column(j).len(), j) < (self.a.column(k).len(), k)
                        || !self.linear_column(j)
                    {
                        continue;
                    }
                    if let Some(column) = self.dominated_pair(&mut scratch, j, k, &mut work)? {
                        fixed += 1;
                        if column == k {
                            break;
                        }
                    }
                }
            }
            Ok(())
        })();
        self.dominated_scratch = scratch;
        result.map(|()| fixed)
    }
}
