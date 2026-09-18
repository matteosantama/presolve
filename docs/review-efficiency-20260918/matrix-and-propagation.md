# Efficiency review: working matrix and propagation hot loops

Scope: `src/matrix/linked.rs`, `src/model/activity.rs`, `src/model/queues.rs`,
`src/rules/bounds.rs`, `src/rules/rows.rs`, `src/rules/variables.rs`, and the
parts of `src/model/mod.rs` they touch. Branch `new-presolve-rules`, commit
`2528782`. Default preset, `threads: 1`.

Method: read the code; counted events with relaxed atomic counters in a copy
of the tree (no timing was run, per the review rules); read the release
disassembly of the owner's `target/release/benchmark` (built 08:55 today from
this commit) to see what is and is not inlined; `size_of` for the structs.
Every count below is a corpus total over the 236 problems unless a problem is
named. Estimates of time are derived from the counts and the rule-level
profiles in `docs/benchmarks/new-rules-20260917/*` (propagation about 26-31%
of `run`, `set_bounds` 36% inclusive on BOYD2, linked traversal about 31%
inclusive on the Netlib set); they are order-of-magnitude and each finding
says what to measure.

## Summary

| rank | finding | est. corpus impact | risk to reductions | effort |
|---|---|---|---|---|
| 1 | `#[inline]` the hot cross-module helpers: `Model::activity` (12.0M out-of-line calls), `Activity::replace_bound` (10.7M), `ActivityRows::push` (9.1M), `Worklist::push` (4.8M), `Queues::column_changed`, `Model::changed_row` | 2-5% | none | trivial |
| 2 | One 64-byte row record (`activity`, `RowDomain`, kind byte, queue flags) instead of four scattered per-row arrays | 3-6% (BOYD2, BOYD1, CONT-*, PILOT*, QAP15) | none (layout only) | medium |
| 2b | Cheap first step of 2: fold `row_kinds` into `ActivityRows.flags` | 0.5-1% | none | low |
| 3 | `implied_bound`: test "not tighter than the current bound" before the gain arithmetic; 6.40M of 7.87M calls (81%) end that way | 0.5-2% | none (proved equivalent below) | trivial |
| 4 | `from_columns`: compute both column links from a one-entry lookahead, removing the random read-modify-write of the previous column node per nonzero; `offsets` as `u32` | 0.5-1% (build is 5-7% of the corpus) | none (identical arena) | low |
| 5 | `Activity::compute`: separate no-exclude and exclude instantiations; branch-free `Extreme::add` | 0.3-1% | none (argument below) | low |
| 6 | `fix`: drop the `Vec<(usize, RowDomain)>` collect (out-of-line `GenericShunt::next` per entry) with a check-then-apply pass | <0.5% | none | low |
| 7 | `singleton_rows`: check the row length before taking the `Arc<Equation>` snapshot | <0.2% | none | trivial |
| - | Cross-scope notes: `Arc` refcount traffic on 667k `TightenedBound` records; column fingerprints by row sweep (parallel.rs, 12.2M scattered steps); `Rule` at 80 bytes | 1-3% combined | none / needs care | medium |
| - | Measured and rejected: precomputing activities at build, merging double-side tightenings, local activity copy in the entry loop, arena compaction, queue changes | 0 | - | - |

## Counts that shaped the findings

Struct sizes (bytes): `Node` 32, `List` 12, `Activity` 32 (`Extreme` 16 =
`f64` + `usize`), `Locks` 16, `RowDomain` 24, `Bounds` 16, `Rule` 80,
`Option<Arc<Equation>>` 8, `Equation` 48.

Corpus totals (default preset; dual propagation, dominated columns and
equality dependencies are off by default, so they contribute nothing here):

| event | count |
|---|---|
| nonzeros built (`from_columns`) / packed | 4.47M / 3.56M |
| node visits, row axis / of which `next == id + 1` | 45.9M / 39.3M |
| node visits, column axis (always scattered: 0.15% within one node of the previous) | 40.5M |
| column visits by rule: parallel_columns / sparsification / propagation (`set_bounds`) / redundant bounds / singleton columns / cleanup | 12.2M / 12.2M / 8.05M / 7.55M / 0.13M / 0.23M |
| row visits by rule: propagation / parallel_rows / sparsification / build / pack / short equalities / singleton columns | 13.8M / 11.0M / 5.6M / 4.47M / 3.56M / 2.9M / 2.6M |
| `set_bounds` calls / incident rows visited | 669k / 8.18M |
| `relax_bound` calls / incident rows | 554k / 2.52M |
| `bound_implied` calls / rows visited / early `true` | 1.28M / 5.04M / 554k |
| `propagate_bounds` rows examined / entries scanned / skipped by the two-infinite shortcut | 883k / 7.89M / 237k |
| `implied_bound` calls: huge or infinite value / crossing but within tolerance / not tighter than the current finite bound / insufficient gain / accepted | 7.87M: 752k / 2 / **6.40M** / 49k / 667k |
| `Model::activity` cached hits / recomputations (entries) | 11.0M / 963k (6.19M) |
| recompute entries: first computation of an untouched row / of a row already shifted incrementally / later recomputations | 2.23M / 2.00M / 1.95M |
| residual recomputations in propagation (both sides for one entry) | 193k (0) |
| stale activity found mid-row in propagation | 9.6k |
| both sides tightened on one entry visit | 5.1k of 647k tightenings |
| `ActivityRows::push` / new in propagation list / new in singleton list | 9.11M / 1.14M / 1.16M |
| `Worklist::push` / inserted | 4.77M / 2.43M |
| `changed_row` calls; `fix` calls (entries); `replace_row` calls (old, new entries); `replace_row_bounds` calls (entries) | 969k; 51k (133k); 82.5k (1.04M, 0.84M); 2.1k (15.7k) |
| `equation` snapshots (entries copied) / cache hits | 188k (1.92M) / 36k |
| linked inserts / removes | 24k / 455k |
| round drains (`take_round`, `drain`) | 10.5k |

The four largest problems:

| | BOYD1 | BOYD2 | CONT-300 | CONT-200 |
|---|---|---|---|---|
| nnz | 559k | 424k | 449k | 198k |
| `set_bounds` incident rows | 1.48M | **4.45M** | 66k | 197k |
| `bound_implied` rows + `relax_bound` rows | 881k + 512k | 390k + 386k | 716k + 440k | 357k + 99k |
| propagation entries scanned | 1.58M | 424k | 503k | 198k |
| `implied_bound` calls / not tighter | 1.38M / 703k | 288k / 209k | 1.01M / 991k | 396k / 356k |
| `Model::activity` hits | 2.14M | 814k | 1.15M | 555k |
| `ActivityRows::push` | 1.48M | 4.64M | 158k | 237k |
| column node visits (propagation / parallel cols / sparsify / redundant) | 7.5M (1.48M / 1.21M / 3.43M / 1.39M) | 6.36M (4.45M / 0.94M / 0.19M / 0.78M) | 2.2M (0.06M / 0.98M / 0 / 1.16M) | 1.09M |
| row node visits | 8.2M | 2.9M | 4.6M | 1.8M |

## Findings

### 1. Cross-module helpers called out of line from the hot loops

Files: `src/model/mod.rs:294` (`activity`), `src/model/activity.rs:211`
(`replace_bound`), `src/model/queues.rs:100` (`ActivityRows::push`),
`src/model/queues.rs:20` (`Worklist::push`), `src/model/queues.rs:173`
(`Queues::column_changed`), `src/model/mod.rs:271` (`changed_row`).

What happens now. The release disassembly of `propagate_bounds` contains two
`bl` calls to `Model::activity` (the row header at `bounds.rs:27` and the
per-entry call at `bounds.rs:73`), and `bound_implied` one (`bounds.rs:159`).
`set_bounds` calls `Activity::replace_bound`, `ActivityRows::push`,
`Worklist::push` (twice) and `Queues::column_changed` out of line. None of
these carry `#[inline]`; they sit in `model/` while the callers sit in
`rules/`, so with the default 16 codegen units and no LTO they are not
inlined (the same hazard the refactor report traced a 5% regression to for
`Extreme::implied` and friends, fixed with `#[inline]`).

Why it costs. `Model::activity` is 82 instructions with four bounds checks and
returns a 32-byte struct through memory; it is called 12.0M times, 7.9M of
them once per scanned propagation entry. `replace_bound` (115 instructions)
is called 10.7M times, once per incident row of every bound change.
`ActivityRows::push` (55 instructions) 9.1M times, 87% of which find the row
already queued and only touch the flag. `Worklist::push` (71 instructions)
4.8M times. That is about 37M calls, each with prologue, argument spills,
the out-pointer store/reload for `Activity`, and a lost opportunity to hoist
loads (`self.bounds`, `self.settings`) across the call. The structures report
measured `Worklist::push` alone at 17% self time on BOYD2 before the
`ActivityRows` merge.

Change. Add `#[inline]` to the six functions above (`activity`,
`replace_bound`, `ActivityRows::push`, `Worklist::push`,
`Queues::column_changed`, `changed_row`); leave `Activity::compute` and
`tighten_bound` out of line (cold or large). Optionally split
`Model::activity` into an inlined fast path and a `#[cold]` `recompute`.

Expected effect. 2-5% of the corpus, concentrated on BOYD1/BOYD2/CONT-*/
PILOT*/QAP15 where the per-entry and per-incident-row calls dominate. This is
the cheapest item in the review.

Risk. None to reductions (no semantic change). Interaction: if the
compiler-profile review adopts `codegen-units = 1` or LTO this inlining happens
anyway; the attributes keep it robust under the default profile and in
`cargo test`.

Measurement. Fresh-process corpus timing, three runs each, min and median,
as in `docs/benchmarks/refactor-20260918/REPORT.md`; the in-process repeated
timing is expected to be flat, as it was there.

### 2. One row record instead of four scattered per-row arrays

Files: `src/model/mod.rs:38-58` (`rows`, `activities`, `row_kinds`),
`src/model/queues.rs:82-86` (`ActivityRows.flags`), consumers at
`src/model/mod.rs:418-444` (`set_bounds`), `src/rules/bounds.rs:154-176`
(`bound_implied`), `src/rules/bounds.rs:21-27` (propagation row header),
`src/rules/substitution.rs:186` (singleton scan), `src/model/mod.rs:271-289`
(`changed_row`).

What happens now. Per incident row of a bound change, `set_bounds` touches
the column node (32 B, scattered), `activities[i]` (32 B, scattered),
`changed_activities.flags[i]` (1 B, scattered) and `row_kinds[i]` (1 B,
scattered), plus `rows[i]` (24 B) for cone rows: four cache lines per row,
8.18M rows. `bound_implied` touches node, `rows[i]` and `activities[i]`:
three lines per row, 5.04M rows. The propagation row header touches
`rows[i]`, `lists[0][i]`, `activities[i]` and the flag byte cleared by
`drain`: four lines, 1.14M rows. The singleton scan touches `rows[i]` and
`row_singletons[i]` for 1.16M rows. `changed_row` touches `rows`, `activities`,
`equations`, `lists[0]`, `row_kinds` and the queue flags: six lines, 969k
calls. Column-axis traversal is the dominant memory pattern of the presolve
(40.5M scattered node visits), and every one of those visits in the model
code is followed by two or three more scattered lines for the row.

Why it costs. On problems whose row arrays exceed L1 (BOYD2: 186k rows, so
`activities` alone is 6 MB, `rows` 4.5 MB, the node arena 13.5 MB) each extra
line is an L2/SLC access of 5-20 ns. BOYD2 does 4.45M `set_bounds` row
visits; two avoidable lines per visit at 5-10 ns is 45-90 ms of the roughly
45 ms the problem takes, which says the accesses mostly overlap in the
out-of-order window; the honest estimate is a fraction of that.

Change. `struct RowState { activity: Activity, domain: RowDomain, kind: u8,
flags: u8 }` with `#[repr(align(64))]`, one `Vec<RowState>` on the model: 32 + 24 + 2 = 58
bytes of payload, padded to one cache line, versus the same 58 bytes spread
over four arrays today. (Flattening `Activity` to two `f64` sums and two
`u32` counts would make it 24 bytes and leave room for `row_singletons` or
the equation slot, but it is not needed for the 64-byte record.) `ActivityRows` keeps
its two entry lists and their first-push order; the flag logic moves to a
method taking `&mut u8`, so `set_bounds` reads and writes one line per row
for activity, kind and both queue flags, and the cone block comes from the
same line. `self.rows[i]` becomes `self.rows[i].domain` at 42 sites (grep
count); `RowDomain` stays the value type used by `shifted`,
`Locks::contribution` and the rules, so the rule code is a mechanical rename.

Expected effect. About 25M fewer scattered line touches corpus-wide
(16M in `set_bounds`, 5M in `bound_implied`, 2-4M in the row headers and
`changed_row`). 3-6% of the corpus, most of it on the seven problems above
20 ms; near zero on problems whose arrays fit L1.

Risk. None to arithmetic or queue order: the same values are read and
written, and each `ActivityRows` list keeps its own push order. The
`queues.rs` unit tests need adapting to the moved flags.

Effort. Medium: about 150 lines, one afternoon, plus `run size` on both
presets.

Measurement. This needs a prototype and a fresh-process corpus comparison;
BOYD2, CONT-300 and QAP15 are the sentinels. Suggested order: do 2b first
(below), measure, then the full record.

2b. Cheap first step: put the row kind in bits 2-3 of `ActivityRows.flags`
(`queues.rs:88-89`), set by `Queues::row_changed`, which already receives
the length and equality flag (`mod.rs:287`; pass the cone flag too), and have
`push` return the kind. `set_bounds` then makes one scattered byte access
instead of two per incident row (8.18M), `row_kinds` disappears. Same
guarantees, about 40 lines, 0.5-1%.

### 3. `implied_bound` computes the gain threshold before finding out the value is not tighter

File: `src/rules/bounds.rs:180-234`.

What happens now. Of 7.87M calls, 6.40M (81%) carry a finite, non-huge value
that is not tighter than the current finite bound (BOYD1 703k of 1.38M,
CONT-300 991k of 1.01M, FIT2D 270k of 271k). Each of these executes, in
order: the huge test, the load of `old`, the crossing test, then the gain
block at lines 212-226 (`side.value(old).is_finite()`, `value != opposite`,
the subtraction, two multiplications, an `abs`, a `max` and the compare) and
only then returns `Ok(false)` from the gain test.

Change. After the crossing test at line 200 and before the gain block,
return `Ok(false)` when `(side == Lower && value <= old.lower) || (side ==
Upper && value >= old.upper)`.

Equivalence. Keep the crossing test first, so nothing changes for inputs with
inverted bounds. Once `lower <= upper` holds, a not-tighter value has `gain <=
0`; the threshold is `max(factor * feasibility, relative * |old|)` with
`factor` and `relative` sanitized to finite and `>= 0`
(`settings.rs:274-287`) and `feasibility >= 0`, so `gain <= threshold` holds
and the gain block returns `Ok(false)` today; with `propagation == false`
the gain block is skipped and line 227 returns `Ok(false)`. Either way the
result and the model state are identical, and no `proof` closure runs in
either version.

Expected effect. About ten fewer floating-point and branch operations on 6.4M
calls; 0.5-2% of the corpus if the operations are on the critical path,
possibly less because they overlap with the loads of the next entry. Cheap
enough to do together with finding 1 and measure once.

Risk. None. Effort: three lines.

### 4. `from_columns` writes the previous column node on every nonzero

File: `src/matrix/linked.rs:185-211`.

What happens now. The second pass places node `id` (a 32-byte random write
into the row-major arena), reads and increments `offsets[row]` (8 bytes,
random), and then, at line 206, `out.nodes[col_list.tail as usize].next[1] =
id`: a random read-modify-write of the node placed for the previous entry of
the same column, which lives in a different row block. Three random accesses
per nonzero, 4.47M nonzeros, on a pass the schedule report timed at 2-3 ms
per large problem out of a build of 5-7% of the corpus.

Change. The next entry of column `j` (row `i'`) will be placed at
`offsets[i']` as it stands now, because no other placement happens between
the two. Iterate each column through `.filter(|(_, v)| *v != 0.0).peekable()`,
set `next[1] = peek().map_or(NONE, |(i', _)| offsets[i'])` and `prev[1] =
column_tail`, and write the node once with all six links. The `nodes[tail]`
write disappears; `offsets` can be `Vec<u32>` (ids are asserted below
`u32::MAX`), halving its footprint on the random access. The arena, ids and
links are identical, so `pack`, the cursors and every traversal order are
unchanged.

Expected effect. One of three random accesses per nonzero removed from the
node-writing pass: 0.5-1% of the corpus, largest on BOYD1/BOYD2/CONT-300
(6-7 ms builds). The `vec![Node; nnz]` fill (about 2% of the largest builds
per the structures report) stays; a zeroed default would need `NONE == 0`,
which is a larger change than it is worth.

Risk. None. Effort: low, with the existing dense-vs-linked test covering it.

### 5. `Activity::compute` per-entry branches

File: `src/model/activity.rs:225-240`, `78-84`.

What happens now. 8.77M entries pass through `compute`. 6.19M come from
`Model::activity` with `exclude == None`, yet each entry evaluates `Some(j) ==
exclude` (an `Option<usize>` compare, 16 bytes). Each of the two
`Extreme::add` calls branches on `is_finite`; rows mixing finite and
infinite terms (free variables among bounded ones) mispredict.

Change. Two instantiations: `compute(row, bounds)` and
`compute_excluding(row, bounds, j)`; `Model::activity` uses the first,
`residual_activity` the second. In `add`, `let finite = term.is_finite();
self.sum += if finite { term } else { 0.0 }; self.infinite += usize::from(!finite);`.

Equivalence of the branch-free `add`. `sum` starts at `+0.0` and only ever
has finite terms added; a sum of finite terms is `-0.0` only if every operand
is `-0.0`, which `+0.0` as the start excludes, so `sum` is never `-0.0` and
`sum + 0.0 == sum` bit for bit. The `replace` path (`activity.rs:56`) is left
alone: its early returns and cancellation tests are the semantics, not
overhead.

Expected effect. 0.3-1%; `Activity::compute` was 10% of CONT-300 and 19% of
MAROS-R7 in the structures profile, but most of that is the scattered
`bounds[j]` loads, which stay. Needs measurement on CONT-300.

Risk. None. Effort: low.

### 6. `fix` collects a temporary `Vec<(usize, RowDomain)>` through an out-of-line iterator

File: `src/model/mod.rs:496-504`, `517-529`.

What happens now. `collect::<Option<Vec<_>>>()` allocates a 32-byte-per-entry
vector and, per the disassembly, drives the column iterator through
`GenericShunt::next` calls that are not inlined; then `remove_column`
allocates `entries` (needed by the tape) and the loop at 517 zips the two.
51k calls, 133k entries, two allocations per call.

Change. First pass: `self.a.column(column).iter().all(|(i, a)|
shifted(self.rows[i], a * value).is_some())`; on success take `entries =
remove_column(column)` and, in the loop, recompute `shifted(self.rows[i], a
* value).unwrap()` (pure function of the same inputs, so the same domain).
One allocation per fix instead of two, and no shunt.

Expected effect. Under 0.5% (about 50k allocations and 133k out-of-line
iterator steps); relevant to STOCFOR3, WOODW, MAROS-R7. Risk: none. Effort:
low.

### 7. `singleton_rows` snapshots rows before checking their length

File: `src/rules/rows.rs:35-40`.

`equation(i)` allocates an `Arc<Equation>` with a copied row for every popped
row, including rows that are no longer singletons and are then skipped at
line 38. Check `self.a.row(i).len() == 1` before calling `equation`. The
snapshot for a true singleton row is still needed by the tape. Effect under
0.2%; risk none; one line.

## Cross-scope notes (owned elsewhere, listed because the counts came out of this review)

- Tape refcounts (`src/model/tape.rs`, `src/rules/bounds.rs:80-84`). Every
  tightening clones an `Arc<Equation>` (667k atomic increments) and every
  tape drop decrements it; 188k snapshots allocate an `Arc` plus a `Vec` and
  copy 1.92M entries. An arena of equations in the tape with a `u32` index in
  `TightenedBound` would remove both atomics per record and shrink `Rule`
  from 80 bytes (the structures report estimated 2-3 ms on BOYD1 for the
  size alone). About 1-2% combined; medium effort in postsolve code.
- Column fingerprints (`src/rules/parallel.rs:23`, 291-294). `parallel_columns`
  walks every column list (12.2M scattered node visits, the largest single
  column-axis consumer, 1.2M on BOYD1, 0.98M on CONT-300). The fingerprint is
  a two-pass in-order fold (max and first sign, then djb2), so it can be
  computed for all columns by two row-major sweeps of the contiguous arena
  with per-column accumulators, visiting each column's entries in the same
  ascending-row order and giving the same `(u32, u32)`. Whether sequential
  arena reads plus ascending scattered accumulator writes beat the column
  walk depends on the problem shape (BOYD1's 18-row columns are already
  prefetcher-friendly); the sequential executor path would need a
  measurement on CONT-300 and QAP15 before adopting.
- Sparsification's 12.2M column visits (`sparsification.rs:131`) are the
  algorithm's candidate search and were already measured near the memory
  floor; nothing in the matrix layer changes that.

## Measured and rejected

- Precompute activities in the build's lock loop (`mod.rs:181-186`) instead
  of the first propagation round. Only 2.23M of the 4.47M nonzeros belong to
  rows whose activity is first computed without an earlier incremental shift
  (cleanup fixes and singleton-row bounds touch the rest before the first
  propagation); the fused loop would add a `bounds[j]` load and the term
  arithmetic for all 4.47M entries and save the traversal for 2.23M. Net
  about zero, and keeping bit-identical arithmetic needs a third
  "precomputed, untouched" state. Not recommended.
- Merging a lower and upper tightening of the same entry into one
  `set_bounds` walk: the arithmetic and queue order would be identical, but
  it happens 5.1k times in 647k tightenings.
- A second residual recomputation for the same entry never happens (0 of
  193k); a local copy of the row activity in the entry loop saves only the
  `STALE` compare and copy on 7.9M entries, with 9.6k mid-row stales; finding 1
  removes the call overhead, which is the part that costs.
- Arena locality: 85.5% of row steps go to `id + 1`, and the remainder are
  mostly row ends (BOYD2 averages 2.3 entries per row and has zero removals);
  455k removes and 24k inserts against 4.47M nodes leave the row-major
  layout intact. No compaction or reordering is warranted; `pack` reads rows
  sequentially and scatters into CSC, which is the right direction.
- Queues: 87% of `ActivityRows::push` and 49% of `Worklist::push` calls are
  duplicates, but the dedup byte is the only cost and finding 2 puts it on a
  line that is loaded anyway; 254k rows popped by propagation are deleted by
  then and cost one line each; 10.5k round drains allocate 10.5k vectors.
  Nothing to change.
- Row-side changes (`replace_row_bounds`, 2.1k calls) mark the activity stale
  although coefficients and variable bounds did not change; preserving the
  cached value would replace a fresh sum by an incrementally shifted one and
  is deliberately avoided (comment at `mod.rs:744`). Left alone.
- A 24-byte singly linked node: saves 25% of row-traversal bytes but every
  `remove` would need a column walk to find the predecessor (455k removes,
  average column length 10-30, so 5-10M extra scattered steps). Not worth it
  without a measurement, and column steps already cost one line each at 32
  bytes.
- `Locks` as two `u32` (8 bytes): touched about 2M times in `replace_row`
  and `fix`; the line count per access does not change. Skipped.

## Suggested order and verification

1. Findings 1 and 3 together (ten minutes), `cargo test --release`, `run size`
   on both presets (expect 0 changed), fresh-process corpus timing three runs
   each.
2. Finding 2b, then finding 2 as a prototype; same verification, with BOYD2,
   CONT-300, QAP15 and PILOT87 as sentinels.
3. Findings 4-7 in one batch; they are individually below the noise floor and
   should be judged on the corpus sum.
