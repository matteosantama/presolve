# Efficiency review: substitution, objective/Hessian, sparsification, dependencies

Tree: `new-presolve-rules` at `2528782`. Scope: `src/model/objective.rs`,
`src/matrix/sparse.rs`, `src/rules/substitution.rs`, `src/rules/sparsification.rs`,
`src/rules/dependencies.rs`, `src/rules/variables.rs`, and the substitution,
row-replacement and tape-record paths of `src/model/mod.rs` / `src/model/tape.rs`.

No file under `src/`, `tests/` or `benchmark/` was modified and no timing or
size benchmark was run. Numbers come from three sources:

1. **Counters.** A copy of the tree in the scratchpad was instrumented with
   static counters at the hot spots (calls, entries merged, arena elements
   shifted, snapshot hits/misses, tape record census, ...) and run once per
   problem over the full corpus for the default and the aggressive preset
   (one presolve per problem, no trials, `nice -n 19`). Corpus totals are in
   the appendix; per-problem columns are quoted inline.
2. **Stored timings.** Per-problem times are the minimum over the four stored
   runs in `benchmark/results/plain-csc-20260918/time/new-{1..4}.jsonl`
   (sum of per-problem medians 527 ms; BOYD1 78.3, BOYD2 46.7, CONT-300 35.2,
   CONT-200 16.1, MAROS-R7 15.0, FIT2P 12.1, STOCFOR3 9.1 ms).
3. **Two micro-measurements** (scratch, `rustc -O`): `copy_within` of
   16-byte elements in runs of 22 to 500 costs 0.6 to 1.5 ns per element;
   a `binary_search_by_key` on a 111-entry row costs 6.5 to 9 ns per lookup,
   and a monotone cursor walk over the same row costs about 2 ns per pair.

Percentages below are of the 527 ms corpus total unless a problem is named.
Estimates are honest orders of magnitude; the "settles it" column says what
measurement would confirm each.

## Summary

| rank | finding | est. impact (default preset) | risk to identical reductions | effort |
|---|---|---|---|---|
| 1 | `short_equalities`: column degree loaded and row maximum computed before the cheap structural rejections (1.26M wasted scattered loads, 305k extra row walks) | 3–6 ms (0.6–1.1%); CONT-300 4–7%, CONT-200/201 5–8% | none | trivial |
| 2 | Equation snapshots: 188k `Arc<Equation>` (376k mallocs) for the tape; store snapshots in one arena with `(start, len)` handles | 5–9 ms (1–1.7%); CONT-200/201/100 5–12%, STOCFOR3 5–9% | none (storage only) | moderate |
| 3 | `try_substitute` pair scan: one binary search per candidate pair (996k `p.get`); walk each affected row with a cursor instead | 3–4 ms (0.6–0.8%); DUAL1–4 lose about half their 5.5 ms; CVXQP2_L, AUG3DC, Q25FV47 | none | small |
| 4 | `SymmetricMatrix::remove_variable` builds a `Vec` that both callers discard (84k mallocs, 170k copied entries) | 2–3 ms (0.4–0.6%); WOODW, MAROS-R7, CVXQP2_L, STCQP1/2 | none | trivial |
| 5 | `replace_row` calls `column_changed` for every retained column even when its coefficient and lock contribution are unchanged (587k of 1.04M) | 2–4 ms (0.4–0.8%); UBH1, DEGEN3, GREENBEA/B about 10% | medium: changes first-push order inside a round; needs `run size` on both presets | small |
| 6 | Singleton round: `rows[i]` (24 B, scattered) read for 937k rows; `row_singletons(i) == 0` already implies the row is not deleted | 1–2 ms (0.2–0.4%); BOYD2 about 0.8 ms, CONT-300 0.4 ms | none | trivial |
| 7 | Tape `Rule` is 80 B while 85% of records are `TightenedBound` (32 B payload); pack the large variants to reach 56–64 B without boxing, and pre-reserve | 1–3 ms (0.2–0.5%); BOYD1 1–2 ms | none | moderate, mechanical |
| 8 | `sparsify_rows` allocates and zeroes a 24-B `Scan` per row up front (17.8 MB memset over the corpus, wasted on problems with no 10-entry row) | 1–2 ms (0.2–0.3%); BOYD2, CONT-* | none | trivial |
| 9 | `fix` collects the shifted domains in a temporary `Vec` (51k mallocs) | about 1 ms (0.2%) | none | trivial |
| 10 | `singleton_range` recomputes the full residual activity (1.29M entries) and 63% of the results are two infinite extremes | 2–4 ms possible, but every cheaper formulation changes floating-point order | not identical; listed for completeness | – |
| A1 | Aggressive only: Hessian arena applies staged updates one at a time (13.5M `set`, 3 binary searches each, 248M elements shifted = 4 GB memmove, 744k relocations) | large on AUG2DC/AUG2D/CVXQP1_L/CONT-*; nil by default (10k sets) | none | moderate |
| A2 | Aggressive only: `equality_dependencies` tests each row against every basis row by binary search (1.5M searches for 12.5k hits) | CONT-300/201, QAP* under aggressive | none | small |
| A3 | Aggressive only: `short_equalities` re-examines rows on every bound change (4.16M examinations, 13.3M candidate checks for 243k attempts; CONT-200 alone 2.3M) | large under aggressive | needs a structural-revision argument | moderate |

Measured but not fixable without changing what is stored or computed:
Hessian tail shifts on removal (3.7M elements, about 1.1 ms of Q25FV47's
4.45 ms and 0.5–0.7 ms of STCQP1/2's 2.5 ms), the sparsification column
scan (15.9M units, about 25 ms of BOYD1), and the snapshot copies the tape
genuinely needs. See "Inherent costs" below.

If findings 1–4 and 6–9 all land, the structural saving is roughly 18–30 ms,
3.5–5.5% of the corpus, with no change to any reduction. Finding 5 adds
2–4 ms if the size comparison stays clean.

## Findings

### 1. `short_equalities` pays for candidates it rejects structurally

`src/rules/substitution.rs:45-58`.

```rust
let largest = row.iter().map(|(_, a)| a.abs()).fold(0.0, f64::max);   // line 45: full row walk
candidates.clear();
for (j, a) in row {
    let degree = self.a.column(j).len();                                // line 48: scattered 12-B load
    if (options.require_free_variable && self.bounds[j] != Bounds::FREE)
        || (options.require_linear_variable && !self.objective.p.column(j).is_empty())
        || degree < options.min_column_length
        || degree > options.max_column_length
```

What happens: every equality row of length 3..8 is examined at least once
(the constructor's `changed_row` seeds the queue with all of them:
332k rows examined, 305k pass the length filter). For each such row the
row is walked once to compute `largest`, then walked again for the
candidate loop, and for every entry the column list header is loaded
*before* the `require_free_variable` test, which rejects 89% of them
(1,261,333 of 1,420,815 candidate entries in the default preset;
CONT-300 446,108 of 446,108, CONT-200 198,005, CONT-201 197,408).

Why it costs: the column headers are `List` records of 12 bytes indexed by
column; for CONT-300 that is a 1.1 MB array touched in row order, so each
load is a scattered L2 access (about 4–5 ns). The second row walk is
sequential (nodes are laid out in row order) but still 2–3 ns per node. Per
CONT-300 pass: 446k × (5 + 2.5) ns ≈ 3 ms upper bound; realistically
1.5–2.5 ms of its 35 ms. Corpus: 1.26M wasted loads plus 1.4M extra node
visits, 3–6 ms.

Change: (a) test `bounds[j]` and the Hessian column before touching the
constraint column, and load `degree` only when both pass; (b) make
`largest` lazy (`Option<f64>` filled by `get_or_insert_with` on the first
candidate that reaches the pivot test at line 57). Both are pure
reorderings of side-effect-free reads; `structural_rejections` still counts
each rejected candidate once, so `EqualityStats` are unchanged too.

Effect: CONT-300 4–7%, CONT-200/201 5–8%, corpus 0.6–1.1%.
Risk: none. Effort: 10 lines.
Settles it: `run time --problem CONT-300` A/B (the effect is above the 5%
per-problem noise there).

### 2. Equation snapshots: one `Arc<Equation>` plus one `Vec` per snapshot

`src/model/mod.rs:257-269` (`equation`), `:562-583` (`substitute_on_side`),
`src/model/tape.rs:60-64`, `:139-144`, `:120-126`.

What happens: 188,361 snapshots are created in the default preset
(`EQ_MISS`), each an `Arc::new(Equation { row, entries: row.to_vec(), bounds })`,
i.e. two mallocs (the 64-byte Arc block and the entries buffer) and one row
traversal. The cache added yesterday pays off where it can (BOYD1: 18
snapshots serve 244k `TightenedBound` records; CONT-300: 8k hits vs 5.7k
misses), but most rows are snapshotted exactly once because the tape needs
one copy per row version: CONT-200 39,601 snapshots for 39,601 tightenings,
STOCFOR3 15,711, CONT-100 9,801, MAROS-R7 3,136 (one per substitution,
145k entries). Only about 11k snapshots are wasted (rejected substitutions
that had already snapshotted, and `singleton_rows` snapshots that end in a
row deletion without a tightening); the rest are required by postsolve.

Why it costs: 376k small mallocs at 20–30 ns each is 7–11 ms; the 1.9M
entry copies (sequential node reads) are another 2–4 ms but are needed in
any design that copies. CONT-200 spends an estimated 2 ms of its 16 ms
here; CONT-201 and CONT-100 similar fractions; STOCFOR3 about 0.8 ms of 9 ms.

Change: give `RecoveryTape` a snapshot arena:

```rust
pub struct Snapshot { row: u32, start: u32, len: u32, bounds: Bounds }   // 32 B
pub struct RecoveryTape { rules: Vec<Rule>, snapshots: Vec<Snapshot>, entries: Entries }
```

`equation(row)` appends the row's entries once and returns a `u32`
snapshot id (the per-row cache becomes `Vec<Option<u32>>`, still
invalidated in `changed_row`); `TightenedBound`, `Substituted` and
`Unlocked` store ids instead of `Arc<Equation>`/`Vec<Equation>`;
`substitute_on_side` appends its own snapshot with the chosen side (or, if
the cached snapshot is warm, a new `Snapshot` header pointing at the same
entry range). `recover`/`reduce_point` index `tape.entries[start..start+len]`
with the same arithmetic in the same order. `Rule::TightenedBound` shrinks
from 32 to 24 B of payload.

Effect: 5–9 ms (1–1.7%), concentrated on CONT-200/201/100 and STOCFOR3;
nothing on BOYD1/BOYD2 (18 and few snapshots).
Risk: none to reductions (the same numbers are stored and read); moderate
code churn in `tape.rs`, `model/mod.rs`, `rules/bounds.rs:82,95`,
`rules/rows.rs:35`, `rules/substitution.rs`. Effort: half a day with tests.
Settles it: A/B on CONT-200 alone; if the malloc estimate holds, it is a
≥1.5 ms (10%) drop on a 16 ms problem, far above noise.

A further step (not recommended yet): record `(row, version)` and copy the
row only when it is later modified (`replace_row` has the old row in
`row_scratch`; `fix`/`aggregate` know the removed entry) or at pack time.
That would also remove the 1.9M-entry copies (69% of snapshots belong to
rows never modified afterwards: `EQ_INVALIDATE_LIVE` 58.8k of 188k) but
makes every row-mutation path responsible for tape correctness.

### 3. `try_substitute` looks up every candidate pair by binary search

`src/model/objective.rs:165-223`, specifically line 174
`let old = self.p.get(a.column, b.column);`.

What happens: for each affected column `j` the inner loop visits `l`
in increasing column order (either `j..affected.len()` or a sorted index
list) and calls `p.get`, a binary search over row `a.column`. That is
763,708 pair visits and 995,960 `get` calls in the default preset. The DUAL
family (dense 75–111 column Hessians, one equality row) dominates: DUAL3
does 216,482 lookups in 111 attempts, DUAL4 132,959, DUAL1 109,852, DUAL2
80,694. Every one of those 367 attempts fails with `QuadraticFill`
(`TS_FAIL_QF`) after about 1950 lookups: the Hessian has a few hundred
structural zeros and `max_fill = 64` is exhausted after the 33rd missing
pair, but reaching it means probing most of the dense pairs first. The four
DUAL problems take 5.5 ms in total; the lookups are most of it (the
micro-measurement gives 1.4–2.0 ms for DUAL3's 216k searches alone,
against its 2.09 ms total). CVXQP2_L (83k lookups, 2251 `HessianGrowth`
rejections of 27 visits each), AUG3DC (55k), Q25FV47 (21k) are the other
users.

Change: since `b.column` increases monotonically within the inner loop,
take `let row = self.p.row(a.column)` once per `j` and advance a cursor:

```rust
while cursor < row.len() && row[cursor].0 < b.column { cursor += 1; }
let old = if cursor < row.len() && row[cursor].0 == b.column { row[cursor].1 } else { 0.0 };
```

The same `old` value is read, in the same visit order, so the fill
accounting, the `visits` deadline counter and the staged updates are
identical bit for bit. The closure needs the row slice and cursor passed in
(or the loops inlined); `self.p` is only read during the scan, so the
borrow is fine.

Effect: micro-measured 4× on the DUAL3 shape (2.0 ms → 0.5 ms for 216k
lookups); about 3 ms across DUAL1–4, 0.5–1 ms elsewhere; 0.6–0.8% of the
corpus. Under the aggressive preset the same loop does 42.7M lookups
(AUG2DC 13M), so the gain there is tens of ms.
Risk: none. Effort: 30 lines in `objective.rs`; the existing dense-reference
test covers it.
Settles it: A/B on DUAL3 (2 ms problem, expected −1 ms).

Not recommended: an early `HessianGrowth` reject from `|affected|²` versus
row lengths. The `added` count depends on exact cancellations
(`value == 0.0` for a missing position), so any bound computed without the
values can reject a substitution the current code accepts.

### 4. `remove_variable` returns a `Vec` that nobody reads

`src/matrix/sparse.rs:235-248`; callers `src/model/objective.rs:230` and
`:274` both write `self.p.remove_variable(k);` and drop the result.

What happens: `let entries = self.row(column).to_vec();` allocates and
copies the row (83,911 calls, 170,428 entries in the default preset; WOODW
4,399 calls, MAROS-R7 4,122, CVXQP2_L 4,501/20k entries, STCQP1 3,676/25k,
QSHELL 436/32.6k). `try_substitute` has already copied the same row into
`gradient` at line 111, so every curved fix or substitution copies its
Hessian row twice and mallocs twice.

Change: iterate the slot region by index. After `std::mem::take(&mut
self.slots[column])` the region `[start, start+len)` is dead space that
`set_entry(row, column, 0.0)` on *other* rows never touches (a removal only
shifts inside that row's own slot; it never relocates or compacts), so

```rust
let Slot { start, len, capacity } = std::mem::take(&mut self.slots[column]);
for i in start..start + len {
    let row = self.entries[i].0;
    if row != column { self.set_entry(row, column, 0.0); self.nnz -= 1; }
}
```

is safe and touches the same rows in the same order. Change the signature
to return nothing (the two callers discard it; the unit test at line 271
can read the row before removing). Bit-identical.

Effect: 84k malloc+copy pairs, 2–3 ms (0.4–0.6%); 324k calls under the
aggressive preset.
Risk: none. Effort: 10 lines.

### 5. `replace_row` re-queues retained columns whose state did not change

`src/model/mod.rs:339-356`.

What happens: after the lock merge (which already skips columns whose lock
contribution is unchanged, line 327), `replace_row` calls
`column_changed(j)` for **every** column of the old row and every new
column. `column_changed` loads the column list header (12 B, scattered) for
the length and pushes `unlocked_columns` (a flag byte plus possibly a push).
`RR_UNCHANGED_COLS` counts 587,341 of the 1,044,304 old-row columns as
retained with an unchanged lock contribution: DEGEN3 103,683 of 105,040,
UBH1 91,133 of 93,376, GREENBEA 77,919 of 83,074, STOCFOR3 12,448 of
30,305. These are the long rows that a substitution or a sparsification
batch rewrites while touching only a few of their entries. Under the
aggressive preset it is 59.8M of 64.0M.

Why it costs: two scattered accesses per column, 587k times, about 2–4 ms
(UBH1 and DEGEN3 each spend an estimated 0.5–0.8 ms of their 4.5 ms here).

Change: skip `column_changed(j)` for a retained column when `added ==
removed` in the lock merge. Argument: the consumers of the column queues
are `simple_dual_fix` (reads `alive`, `locks[j]`, `objective.diagonal(j)`,
`c[j]`, `bounds[j]`), `empty_columns` and `singleton_columns` (column
length). None of those change for a retained column with an unchanged lock
contribution: its length is unchanged (still present in this row and no
other row changed), and every mutation of `c[j]`/`P[j,:]` in a
substitution pushes `j` separately (`substitute_equation` line 713,
`fix` line 511, `eliminate_coupled` line 861, `aggregate` line 822).
So the *decision* the queue consumer would take for `j` is the same whether
or not it is re-queued now.

Why the risk is medium anyway: what can change is `j`'s *position* in the
round. If `j` is pushed by a later event in the same round, it now lands at
the later position; if the round contains a `remove_unlocked` or a
Hessian-coupled `fix` before that position, `j` is evaluated with slightly
different locks or `c[j]` than before, and a fix could move to the next
cleanup iteration, where `singleton_rows` may already have tightened its
bound. On the corpus this is probably invisible, but it is not provable
from the code alone. Run `run size` on both presets before keeping it.

Effect: 2–4 ms (0.4–0.8%); UBH1, DEGEN3, GREENBEA/B about 10% each.
Effort: 5 lines. Settles it: size comparison (0 changed) plus A/B on UBH1.

### 6. Singleton round reads `rows[i]` for every queued row

`src/rules/substitution.rs:185-194`.

```rust
for i in self.queues.changed_activities.take_singleton_round() {
    if self.rows[i] == RowDomain::Deleted || self.a.row_singletons(i) == 0 { continue; }
```

What happens: `set_bounds` pushes every incident row, so this loop sees
937,193 rows in the default preset (BOYD2 186,531, CONT-300 103,144,
CONT-201 49,340, CONT-200 39,601, STOCFOR3 34,523) and rejects 93% of them
(66,289 scanned). The `rows[i]` load is a 24-byte scattered read (4.5 MB
array on BOYD2) made before the 4-byte `row_singletons` read, and it is
redundant: a deleted row has no entries, and `LinkedMatrix::remove` keeps
`row_singletons` exact, so a deleted row's count is zero.

Change: `if self.a.row_singletons(i) == 0 { continue; }` (or at least test
the count first). Same rows scanned in the same order.

Effect: BOYD2 about 0.8 ms (1.7%), CONT-300 0.4 ms, corpus 1–2 ms.
Risk: none. Effort: one line.

### 7. Tape record size and growth

`src/model/tape.rs:85-174`, `:200-203`; `size_of::<Rule>() == 80`.

What happens: 780,524 records are pushed in the default preset, 667,122 of
them `TightenedBound` (payload `usize + Arc + Side + f64` = 25 B, 32 with
padding). The enum is 80 B because `Substituted` (8 + 8 + 32 + 24 + 1) and
`Fixed` (8 + 8 + 32 + 24) are 72–73 B. `rules` is a plain `Vec` grown by
doubling. BOYD1 pushes 244,376 records: 19.5 MB written plus about the same
copied through 18 doublings; the profile from yesterday attributes 4% of
BOYD1 to tape write and drop.

Change (no boxing; boxing `Substituted`/`Fixed` would add 82k mallocs and
eat the gain, as yesterday's report also concluded):

- `Substituted`: store `gradient.terms` and `other_rows` in one `Entries`
  with a `u32` split, and `gradient.constant` inline: 8 + 8 + 8 + 24 + 4 + 1
  = 53 → 56 B. `Fixed` likewise (52 → 56). `RowCombination` (64) can carry
  `activity` as `Option<u32>` column plus the terms in the same buffer as
  `targets` (56). With finding 2, `TightenedBound` is 24 B and `Rule`
  becomes 56–64 B.
- `RecoveryTape::default()` → `Vec::with_capacity(rows + columns)` in
  `from_parts`; untouched capacity costs nothing on macOS/Linux (lazy
  commit) and removes most doublings (BOYD1 still needs two).

Postsolve reads the same values; only where they sit changes.

Effect: BOYD1 1–2 ms (1.5–2.5%), corpus 1–3 ms.
Risk: none. Effort: mechanical, every `match` in `tape.rs`.

### 8. `sparsify_rows` scan table is allocated eagerly and is 24 B per row

`src/rules/sparsification.rs:21-25`, `:84`.

What happens: `vec![Scan::default(); m]` is allocated and filled (the
default is not all-zero, so it is a real fill, not calloc) before the
reference loop: 741,783 records × 24 B = 17.8 MB over the corpus. On BOYD2
(186,531 rows, 4.5 MB) only 101k candidate visits follow; on CONT-200
(39,601 rows) no reference row has 10 entries, so the 0.95 MB fill is the
whole cost of the rule.

Change: allocate on the first reference row that passes the length filter
(`SP_REF_LONG` is 27,694 rows across the corpus vs 695,847 references
examined), and shrink `Scan` to 16 B (`seen: u32, count: u32, ratio: f64`;
`m < u32::MAX` already holds for the linked matrix). The scan then touches
denser lines during its 1.5M candidate visits. Identical results.

Effect: 1–2 ms; BOYD2 0.3–0.5 ms. Risk: none. Effort: trivial.

### 9. `fix` allocates a temporary for the shifted domains

`src/model/mod.rs:496-504`.

`updates: Vec<(usize, RowDomain)>` is collected only to check that every
`shifted` result is `Some` before mutating, then zipped with the removed
entries. `shifted` is a pure function of `(rows[i], a * value)`, so a first
pass that only checks and a second pass that recomputes it while applying
gives the same domains and drops 50,996 mallocs (default preset).

Effect: about 1 ms. Risk: none. Effort: trivial.

### 10. `singleton_range` recomputes the residual activity from scratch

`src/rules/substitution.rs:152-163` → `residual_activity` → `Activity::compute`.

What happens: 70,535 calls walk 1,285,861 entries (FIT2P 24,000 calls /
383k entries for 3,000 substitutions; MAROS-R7 7,443 / 342k; 80BAU3B
75k). In 63% of the calls both residual extremes are infinite
(`SC_RANGE_BOTHINF` 27,565 of 43,885 on the sampled problems; FIT2P
23,974 of 24,000), so the arithmetic decides nothing.

Why nothing is proposed: the only cheaper formulations change
floating-point order. Using the cached row activity with
`Extreme::excluding` (as propagation does) computes `sum − term` instead
of the direct sum and can move an implied bound by an ulp, which changes
which sides count as implied. Prefix/suffix sums have the same problem. A
structural shortcut (skip when the exclusive infinite count is ≥ 1 on
both sides) would be exact, but the cache is cold when
`singleton_columns` runs (3,747 warm of 43,885 calls; FIT2P 0 warm) and
maintaining separate exact infinite counts per row costs a bounds load per
nonzero in `from_parts`, roughly what it saves. Listed so nobody re-derives
it; if the owner ever accepts an FP-order change here, this is worth
2–4 ms.

## Aggressive-preset findings

These do not move the default corpus number but dominate several
aggressive-preset problems.

### A1. Hessian arena: apply staged updates per row, not per entry

`src/model/objective.rs:235-237`, `src/matrix/sparse.rs:149-211`.

Counters (aggressive): 13.5M `set` calls, each doing `get` (binary
search) and then `set_entry` twice (each another binary search — the first
repeats the search `get` just did); 15.1M in-place inserts and 248M
elements shifted by `copy_within` (4 GB at 16 B; at the measured
0.6–1.5 ns per element that is 0.15–0.4 s); 744k slot relocations copying
22.4M entries; 66 compactions copying 14.8M; arena peaks 11M entries on
AUG2DC (176 MB) with 5.3M dead. AUG2DC, AUG2D, CVXQP1_L, CONT-300 carry it.

Change: `scratch.quadratic` is produced in `(j, l)` order with `l ≥ j`
over the sorted `affected` list, in two runs (`missing` true, then false).
Group the updates by row (for row `j`: the run of `(j, l)` entries; for
row `l`: the mirrored `(j, l)` entries with `j < l`, gathered by a counting
pass over `affected` indices) and merge each row's sorted update list into
its slot in one pass: O(len + updates) per row instead of O(updates ×
(log len + len/2)). The final arena content is a set of `(row, col, value)`
and each pair is staged at most once, so the order of application does not
change the result; rows stay sorted; `nnz` and `revision` bookkeeping is
the same. Also fold `get` into `set_entry` (return the old value from the
first search) for the remaining single-entry path.

Effect: aggressive only; expect AUG2DC/AUG2D to drop by a large fraction
of their Hessian time (they currently hit the 2 s time limit or close to
it). Default preset: 10k sets, nil.
Risk: none to values. Effort: a `merge_row_updates` on `SymmetricMatrix`
plus the grouping in `try_substitute`, with the randomized arena test
extended.

### A2. `equality_dependencies`: pivot membership by stamp instead of binary search

`src/rules/dependencies.rs:102-147`.

For each candidate row the loop tests every basis row (up to 64) with
`row.binary_search_by_key(&reference.pivot, ..)`: 1,514,849 searches for
12,548 hits (CONT-300 528k, CONT-201 228k). Keep a `Vec<u32>` stamp per
column, stamped with the candidate's id when the row is loaded and
re-stamped after each successful subtraction (12.5k times), so the test is
one array read. The `work` decrement at line 106 stays as is, so the
work-limit behaviour is unchanged. Effect: 10–20 ms on CONT-300 aggressive.
Risk: none. Effort: 15 lines.

### A3. `short_equalities` re-examination storm

With `require_free_variable == false`, `set_bounds` (`model/mod.rs:435`)
re-queues every longer equality on any bound change, so the aggressive
preset examines 4.16M rows and 13.3M candidate entries to make 243k
attempts; CONT-200 alone examines 2.3M rows (58 times per row) for 23k
substitutions. A bound change cannot create a candidate (the structural
filter reads degrees, Hessian columns and, in this preset, no bounds) and
does not change acceptance of a substitution (fill and Hessian checks do
not read bounds; only `effective`/`retained` do). A per-row memo of the
last structural revision at which the row produced no candidate would skip
almost all of these, but it needs a structural revision counter on the
model (bumped by every A or P edit) and a careful argument for the
`preserve_nonzeros` path, which does read bounds through `effective`.
Sketched, not designed. Effect: large on CONT-* aggressive.

### A4. Per-row `entries` buffers in `substitute_equation`

`src/model/mod.rs:662`: one `Vec` per other row per substitution.
Default preset 23k (yesterday's "not worth the plumbing" holds); aggressive
1.56M (CONT-300 666k, CONT-201 294k). A `Vec<Entries>` pool on the model,
handed to `replace_row` and returned, removes them. Aggressive only.

## Inherent costs (measured, no bit-identical fix found)

- **Hessian tail shifts on removal.** `remove_variable` deletes column `k`
  from each neighbour's sorted slot with `copy_within`; 3,739,247 elements
  shifted in the default preset (Q25FV47 1,698,770, STCQP1 723,450,
  STCQP2 440,716, QSHIP12S 175k). At the measured 0.6–1.5 ns per element
  that is about 1.1 ms of Q25FV47's 4.45 ms and 0.5–0.7 ms of STCQP1/2's
  2.5 ms; 3–5 ms over the corpus. Every alternative (tombstones, unsorted
  slots, deferred removal across fixes) either stops `row()` from being a
  slice or changes the order in which `gradient(k)` is recorded, which
  postsolve sums in that order. The per-`Vec` layout it replaced had the
  same shift cost, so this is not a regression, just the price of sorted
  contiguous rows.
- **Sparsification column scan.** 15.9M scan units (BOYD1 5.84M, its
  8×(nnz+bounds) work limit, about 25 ms of its 78 ms; FIT2D 1.19M;
  MAROS-R7 0.83M). The unit cost is the scattered linked column traversal
  (yesterday: about 4 ns per unit); the scan visits Σ column-length² over
  the reference row, which is the algorithm. BOYD1 does get 41k nonzeros
  out of it.
- **Snapshots the tape needs.** After finding 2 removes the mallocs, the
  1.9M-entry copies remain unless the tape stores row versions (see the
  note under finding 2).
- **Transactional rejections.** All 6,764 rejected substitutions in the
  default preset fail on the Hessian side (`TS_FAIL` = `SE_REJECT`;
  405 `QuadraticFill`, 3,451 `HessianGrowth`, 0 numerical/constraint-fill
  on the sampled problems), and the A-side merge work thrown away is only
  16,388 entries corpus-wide (Q25FV47 12,904). Reordering the Hessian scan
  before the A merge would save nothing measurable; the expensive part of a
  rejection is the scan itself (finding 3).
- **Doubleton candidates.** 7,793 of 16,775 queued rows are stale (no
  longer a doubleton equality) but the check is two loads; fine.

## Considered and rejected

- Warm-cache shortcut in `singleton_range`: cache is cold there (see 10).
- Early `HessianGrowth` rejection from sizes alone: not identical (see 3).
- Boxing large `Rule` variants: 82k extra mallocs (see 7).
- Splitting `try_substitute` into stage/commit so the A-side merge runs
  only after the Hessian check: A-side rejections are 0 by default and the
  wasted merge is 16k entries.
- Skipping `column_changed` for *value-changed but lock-unchanged* columns:
  the queue consumers do not read coefficient values, so the same argument
  as finding 5 applies, with the same order caveat; it is subsumed.
- Restricting the singleton round push in `set_bounds` to rows with
  singletons: rejected yesterday (order), still rejected.

## Interaction with the compiler profile

Findings 3 and 4 move work out of `SymmetricMatrix::get`/`set_entry`,
which are small cross-module helpers; if the profile agent adds
`#[inline]` there, the residual per-call overhead of the current code
shrinks but the algorithmic counts above do not. Finding 1's lazy loads are
inside one function and unaffected.

## Appendix A: default-preset counter totals

Corpus totals from the instrumented scratch copy (236 problems, one presolve
each). Top contributors in parentheses.

| counter | total | top problems |
|---|---|---|
| `try_substitute` calls / curved / pair visits / `p.get` | 88,929 / 30,881 / 763,708 / 995,960 | DUAL3 216k gets, DUAL4 133k, DUAL1 110k |
| `try_substitute` failures (QuadraticFill / HessianGrowth) | 6,764 (405 / 3,451 on sampled set) | CVXQP2_L 2,251, AUG3DC 1,200, DUAL1–4 367 |
| `SymmetricMatrix::set` / in-place inserts / relocations / compactions | 10,370 / 2,712 / 235 / 2 | Q25FV47 6,489 sets |
| arena elements shifted (`copy_within`) | 3,739,247 | Q25FV47 1.70M, STCQP1 723k, STCQP2 441k |
| `remove_variable` calls / entries copied | 83,911 / 170,428 | CVXQP2_L 4,501, WOODW 4,399, MAROS-R7 4,122 |
| `equation()` hits / misses / entries copied / live invalidations | 35,935 / 188,361 / 1,915,528 / 58,773 | BOYD1 526k entries in 17 misses; CONT-200 39,601 misses |
| `substitute_equation` calls / accepted / rejected / entries merged then discarded | 37,930 / 31,166 / 6,764 / 16,388 | MAROS-R7 3,136, FIT2P 3,000, STOCFOR3 2,801 |
| `replace_row` calls / old entries / retained with unchanged locks | 82,545 / 1,044,304 / 587,341 | DEGEN3 105k/104k, UBH1 93k/91k, GREENBEA 83k/78k |
| `replace_rows_batch` calls / rows / touched columns | 849 / 2,309 / 81,223 | BOYD1 41,551 touched |
| `fix` calls / curved / removed entries / gradient terms | 50,996 / 23,768 / 133,253 / 165,200 | QSHELL 32.6k gradient terms |
| `short_equalities` rows examined / past length filter / candidate entries / degree loaded then rejected on bounds | 332,392 / 305,173 / 1,420,815 / 1,261,333 | CONT-300 89k rows, 446k entries |
| singleton round rows / scanned / entries; singleton candidates / alive; `singleton_range` calls / entries; attempts | 937,193 / 66,289 / 786,956; 113,480 / 103,333; 70,535 / 1,285,861; 28,441 | BOYD2 186k round rows; FIT2P 24k ranges / 383k entries |
| doubleton candidates / stale / attempts | 16,775 / 7,793 / 8,981 | STOCFOR3 3,016 |
| sparsification references / ≥10 entries / scan units / candidates / subtractions / rows rewritten / `Scan` records allocated | 695,847 / 27,694 / 15,936,510 / 1,512,258 / 1,719 / 2,309 / 741,783 | BOYD1 5.84M units; QAP15 372k candidates |
| tape records / TightenedBound / Substituted / Fixed / RowCombination / DeletedRow | 780,524 / 667,122 / 31,166 / 50,996 / 849 / 28,165 | BOYD1 244k tightened |
| unique `Arc<Equation>` on tape / entries | 177,464 / 1,924,469 | |

`size_of`: `Rule` 80, `Equation` 48, `Gradient` 32, `Entries` 24,
`Bounds` 16, `Affected` 32, `Slot` 24.

## Appendix B: aggressive-preset highlights

| counter | total | top problems |
|---|---|---|
| `try_substitute` calls / pair visits / `p.get` / staged updates | 324,983 / 27.3M / 42.7M / 13.6M | AUG2DC 8.7M visits, 13M gets |
| `set` / in-place inserts / elements shifted / relocations (entries copied) / compactions (copied) | 13.5M / 15.1M / 247.8M / 744k (22.4M) / 66 (14.8M) | AUG2DC 57M shifted, CVXQP1_L 56M |
| arena peak / dead peak (entries) | 28.1M / 13.1M summed; AUG2DC 11.0M / 5.3M | |
| `equation()` hits / misses | 2.94M / 372k | CONT-200 2.11M hits |
| `substitute_equation` accepted / rejected / entries merged | 273,214 / 3,545 / 77.3M | CONT-300 63k, UBH1 24M merged |
| per-row `entries` buffers | 1,559,327 | CONT-300 666k |
| `replace_row` old entries / unchanged-lock retained | 64.0M / 59.8M | UBH1 24M, CONT-300 13M, QAP15 7.7M |
| `short_equalities` rows examined / candidates / attempts | 4,164,775 / 13,310,612 / 243,387 | CONT-200 2.30M / 8.17M / 23k |
| dependencies rows / basis tests / pivot hits / removed | 30,965 / 1,514,849 / 12,548 / 44 | CONT-300 528k tests, QAP15 27 removed |
| quadratic eliminations attempted / accepted | 3,768 / 951 | AUG3DC 3,154 / 738 |
| tape records / TightenedBound | 3,862,554 / 3,496,543 | CONT-200 2.15M |

## Appendix C: reproducing the counters

The instrumented copy lives in the session scratchpad
(`scratchpad/inst`, a `pub mod counters` of `AtomicUsize`s bumped at the
sites named above, plus `benchmark/src/bin/counters.rs`, which loads
problems with the benchmark's `data.rs` and prints one TSV column per
problem). It is not part of the tree; rebuilding it takes the patch list in
this session's transcript. Counts are deterministic, so a single run per
preset suffices.
