# Scheduler and per-phase overhead: efficiency pass

Scope: the scheduler (`src/core/schedule.rs`) and the per-phase overhead of the
existing rules. `src/core/rules/dual_propagation.rs` was not touched. All
changes preserve the exact sequence of reductions; the size comparison over the
236 Netlib + Maros-Meszaros problems is identical to the baseline.

Machine: Apple M1 Pro (8 cores), shared with two other engineers benchmarking at
the same time. The load average moved between 5 and 29 during this work, so
every absolute number below is noisy. The per-problem A/B tables use
interleaved runs (A B A B ...) and report the minimum and median of 9 trials;
structural savings (fewer visits, fewer scattered loads) were preferred over
anything that could only be justified by a few percent under this noise.

## 1. Where the time goes (baseline)

Method: temporary `thread_local` timers/counters around every rule call in
`schedule.rs` plus counters inside the hot rules (rows/entries visited, queue
pushes, residual recomputations, fingerprint vs sort split, activity cache
misses). The instrumentation was stripped before producing `optimization.patch`.
A `samply` profile of the whole corpus (symbolicated with `atos`) confirmed the
attribution: MPS parsing dominates the harness process and presolve itself is
about 25% of samples; within presolve the self-time leaders were
`dual_propagation`, `sparsify_rows`, `propagate_bounds`, `Activity::compute`,
`fingerprint`, `remove_redundant_bound`, `singleton_columns`, `set_bounds` and
`Worklist::push`.

Corpus total with rule-level timers only (in-process size run, one call per
problem; `run` = `Model::run` = 1088 ms in that run; nested timers overlap, e.g.
`propagate_rounds` and `sparsify_cleanup` contain `propagate_bounds`):

| rule | ms | % of run | calls |
|---|---|---|---|
| propagate_bounds | 281 | 26% | 903 |
| sparsify_rows | 198 | 18% | 236 |
| dual_propagation (other engineer) | 184 | 17% | 406 |
| remove_redundant_bounds | 104 | 10% | 236 |
| singleton_columns | 98 | 9% | 1028 |
| parallel_columns | 87 | 8% | 406 |
| build_model (matrix + queue seeding, outside `run`) | 81 | 7% | 236 |
| parallel_rows | 52 | 5% | 406 |
| doubleton_equalities | 33 | 3% | 1028 |
| cleanup loop (all cheap rules; 5118 iterations) | 30 | 3% | 4211 |
| short_equalities | 16 | 1% | 1434 |
| pack (outside `run`) | 12 | 1% | 211 |
| coupled_dual_fix | 3 | 0% | 406 |

Counters that shaped the work (corpus totals):

- The cleanup loop is not a problem: 4211 calls, 5118 iterations (1.2 per
  call), 30 ms in total; the queues already make it incremental, so "repeated
  full passes" do not occur. `work_size()` is O(1). Extra propagation rounds are
  rare (104 in the corpus). `sparsify_cleanup` is a fixpoint loop (PILOT87: 16
  propagation calls, GREENBEA: 57) but each iteration only processes queued rows.
- `propagate_bounds`: 885k rows examined, 6.1M entries scanned, 272k residual
  recomputations (2.65M entries), and 1.22M `set_bounds` calls visiting 10.7M
  rows. On BOYD2 the tightened columns average 32 rows, so `set_bounds` (cached
  activity update plus queue pushes per incident row) is most of its 47 ms;
  on BOYD1 it is the 336k tightenings themselves. This is the algorithm's
  inherent work; bit-identical arithmetic rules out shortcuts.
- `singleton_columns`: two thirds of its time was the first loop, which walks
  every entry of every row whose bounds changed (914k rows, 7.67M entries) to
  find columns of length one. 845k of those rows (6.83M entries, 89%) contain
  no singleton column at all.
- `remove_redundant_bounds`: 724k candidate columns, 1.28M one-sided checks,
  553k bounds relaxed. Every relaxation went through `set_bounds`, pushing rows
  into three work queues that no rule drains afterwards (it is the final pass),
  and the candidate sort recomputed the column length in every comparison
  (14 ms of sort).
- Parallel detection: fingerprint + sort 106 ms (64 fingerprint, 42 sort) for
  1.95M items; the proportionality comparisons cost 3.5 ms. Only 16% of the
  fingerprinted items (estimated 14 ms) belong to a repeated medium phase.
- `sparsify_rows` runs at about 4 ns per unit of its work budget (15.8M
  units); it is close to the memory-bound floor of linked column traversal.
- `build_model` on the largest models: ~1 ms counting pass, 2-3 ms node
  writing, ~2 ms lock/queue seeding; EXDATA's 12.8 ms is entirely the dense
  Hessian conversion (`SymmetricMatrix::from_upper_columns`).

Slowest problems in the same run (ms; `run` excludes build and pack):

| problem | run | propagate | sparsify | redundant bounds | parallel rows+cols | singleton cols | dual prop | build |
|---|---|---|---|---|---|---|---|---|
| BOYD1 | 84.9 | 32.4 | 25.4 | 12.0 | 11.6 | 2.5 | 0.5 | 6.3 |
| BOYD2 | 63.2 | 28.4 | 1.4 | 10.6 | 12.4 | 0.8 | 9.1 | 6.9 |
| CONT-300 | 44.8 | 9.2 | 0.3 | 10.9 | 10.0 | 0.8 | 10.4 | 6.6 |
| CONT-201 | 19.3 | 4.5 | 0.1 | 5.4 | 3.8 | 0.4 | 3.5 | 3.2 |
| MAROS-R7 | 52.5 | 0.7 | 4.1 | 1.0 | 4.9 | 20.4 | 21.2 | 2.7 |
| PILOT87 | 86.1 | 37.8 | 5.3 | 12.6 | 1.2 | 3.8 | 25.0 | 0.7 |
| FIT2P | 40.2 | 0.2 | 25.9 | 0.8 | 3.0 | 4.7 | 5.5 | 0.6 |
| FIT2D | 35.9 | 2.4 | 21.6 | 4.6 | 3.5 | 0.6 | 3.1 | 1.2 |
| PILOT | 45.2 | 11.0 | 25.6 | 1.7 | 1.5 | 1.4 | 3.7 | 0.3 |
| GREENBEA | 57.7 | 18.1 | 12.1 | 0.7 | 2.9 | 1.5 | 4.5 | 0.3 |
| DFL001 | 20.3 | 1.9 | 1.0 | 2.0 | 3.9 | 2.2 | 9.0 | 0.5 |
| STOCFOR3 | 15.5 | 2.9 | 0.2 | 2.1 | 2.3 | 1.0 | 5.4 | 0.6 |
| QAP15 | 14.4 | 2.8 | 4.3 | 2.5 | 1.6 | 0.1 | 3.0 | 0.8 |
| EXDATA | 3.0 | 0.1 | 0.1 | 0.1 | 2.6 | 0.0 | 0.1 | 12.8 |

## 2. Changes

### 2.1 Per-row count of singleton columns (`src/matrix/linked.rs`, `src/core/rules/substitution.rs`)

`LinkedMatrix` keeps `row_singletons[i]`: the number of entries of row `i`
whose column has exactly one entry. A column's length crosses one only while
it has at most two entries, so `insert` and `remove` maintain the count in
O(1): on insert, new length 1 increments this row and new length 2 decrements
the row of the previously lone entry; on remove, new length 0 decrements the
removed node's row and new length 1 increments the row of the remaining entry.
`from_columns` initializes it with one pass over the column lists.

`singleton_columns` skips a queued row whose count is zero. That loop's only
effect is `singleton_columns.push(j)` for entries whose column has length one,
so a row without such entries contributes nothing, and the push sequence
(hence the candidate order and every later substitution) is unchanged. The
existing randomized dense-vs-linked test now also checks the counts after
every edit.

Effect (instrumented runs): scan entries 7.67M -> 0.85M, scan time 70 -> 8 ms
on the corpus, `singleton_columns` 108 -> 46 ms. BOYD2 11.3 -> 0.3 ms,
CONT-300 6.9 -> 0.2 ms, PILOT87 2.6 -> 0.1 ms.

### 2.2 `remove_redundant_bounds` (`src/core/rules/bounds.rs`, `src/core/model.rs`)

- The candidate order `(column length, index)` is materialized once and sorted
  as stored keys instead of recomputing `column(j).len()` (a scattered 16-byte
  load) in every comparison. Keys are unique, so the order is identical.
- `relax_bound` (single caller: `remove_redundant_bound`) no longer goes
  through `set_bounds`. It updates the bound, the cached activities of the
  incident rows (the same incremental `replace_bound` with the same
  cancellation fallback, so later columns of the pass see bit-identical
  activities) and the revision. The queue pushes (`changed_activities`,
  `singleton_activity_rows`, `unlocked_columns`, cone blocks, doubleton and
  short-equality re-queues) are dropped: `remove_redundant_bounds` runs after
  `run_phases` returns and nothing drains a queue afterwards, so they were dead
  work. The doc comment states this precondition. `revision` is kept because
  `presolve_owned` uses `revision == 0` to return `Unchanged`.

### 2.3 Packed fingerprint keys (`src/core/rules/parallel.rs`)

The candidate sort compared `((u32, u32), Option<(u32, u32)>)` keys. Packing
each `(support, coefficients)` pair into one `u64` (`support << 32 |
coefficients`) keeps the lexicographic order exactly (and `None < Some`), so
groups and their member order are unchanged, while column elements shrink from
32 to 24 bytes and each comparison is a single word. `fingerprint` and its
tests are untouched; packing happens at the two call sites.

## 3. Safety: same reductions

`cargo run --release -p benchmark -- run size --name <n> --rule all
--results-dir benchmark/results/agent` was taken before any edit (`base`) and
after each change (`c1` = 2.1, `c2` = 2.1 + 2.2, `c3` = final):

```
compare size base c1: 236 matched cases; 0 changed. 0 unmatched or invalid cases.
compare size base c2: 236 matched cases; 0 changed. 0 unmatched or invalid cases.
compare size base c3: 236 matched cases; 0 changed. 0 unmatched or invalid cases.
```

Variables, rows, nonzeros, bound sides and outcome are identical for every
problem. None of the changes alters an arithmetic expression, a candidate
order or a tie-break; the equivalence argument for each is given above.

## 4. Measured times

Interleaved A/B runs of the baseline binary against the candidate,
`run time --rule all --trials 3 --pool-mode reused --problem <p>`, 3 rounds in
alternating order (9 trials per binary); min and median in ms.

Baseline vs changes 2.1 + 2.2:

| problem | base min | new min | ratio | base med | new med | ratio |
|---|---|---|---|---|---|---|
| BOYD1 | 91.01 | 86.96 | 0.955 | 92.66 | 87.87 | 0.948 |
| BOYD2 | 65.97 | 63.13 | 0.957 | 66.91 | 63.99 | 0.956 |
| CONT-300 | 49.81 | 44.78 | 0.899 | 59.71 | 45.56 | 0.763 |
| CONT-201 | 20.88 | 19.47 | 0.932 | 21.63 | 19.84 | 0.918 |
| MAROS-R7 | 24.23 | 24.35 | 1.005 | 25.63 | 24.67 | 0.963 |
| PILOT87 | 19.79 | 17.86 | 0.902 | 20.51 | 18.52 | 0.903 |
| FIT2P | 18.50 | 18.49 | 0.999 | 19.39 | 18.93 | 0.976 |
| PILOT | 14.06 | 12.70 | 0.903 | 14.36 | 12.85 | 0.895 |
| GREENBEA | 12.32 | 11.56 | 0.938 | 13.02 | 12.09 | 0.929 |
| DFL001 | 13.71 | 12.85 | 0.937 | 14.07 | 13.38 | 0.951 |
| STOCFOR3 | 13.48 | 12.99 | 0.964 | 14.27 | 13.33 | 0.935 |
| QAP15 | 13.35 | 12.52 | 0.938 | 13.85 | 12.90 | 0.931 |

Changes 2.1 + 2.2 vs adding 2.3 (packed keys):

| problem | c2 min | c3 min | ratio | c2 med | c3 med | ratio |
|---|---|---|---|---|---|---|
| BOYD1 | 87.50 | 86.37 | 0.987 | 88.06 | 87.68 | 0.996 |
| BOYD2 | 61.96 | 61.28 | 0.989 | 63.95 | 61.99 | 0.969 |
| CONT-300 | 46.37 | 45.32 | 0.977 | 48.23 | 48.22 | 1.000 |
| MAROS-R7 | 24.00 | 24.69 | 1.029 | 25.42 | 25.43 | 1.000 |
| FIT2P | 18.69 | 18.46 | 0.988 | 19.46 | 19.09 | 0.981 |
| DFL001 | 13.08 | 13.16 | 1.006 | 13.47 | 13.79 | 1.024 |

MAROS-R7 and FIT2P are expected to be neutral: their time is in
`dual_propagation`, `singleton_columns` substitutions and `sparsify_rows`,
none of which these changes touch.

Whole-corpus harness comparison (`run time --trials 5 --pool-mode reused`,
sum of per-problem medians):

Two baseline runs (`base-time3`, `base-time4`) and two candidate runs
(`final3`, `final4`) were taken back to back in the order base, final, base,
final (`compare time <a> <b>` sums the per-problem medians of 5 trials):

| comparison | first | second | change |
|---|---|---|---|
| base-time3 -> base-time4 (noise reference, same binary) | 826.5 ms | 724.4 ms | -12.4% |
| final3 -> final4 (noise reference, same binary) | 691.0 ms | 683.4 ms | -1.1% |
| base-time4 -> final4 | 724.4 ms | 683.4 ms | -5.7% |
| base-time4 -> final3 | 724.4 ms | 691.0 ms | -4.6% |
| base-time3 -> final3 | 826.5 ms | 691.0 ms | -16.4% (inflated by the load spike in base-time3) |

Pooling both runs per binary (10 trials per problem): sum of per-problem
minima 706.1 ms -> 663.5 ms (ratio 0.940), sum of per-problem medians
764.6 ms -> 684.8 ms (0.896); 220 of 236 problems have a lower minimum with
the candidate. The first `compare time base-time final` pair, taken while
another engineer's job pushed the load average to 18, is not usable: its
same-binary control (`base-time2`) came out +62% over `base-time`.

Per-problem minima over those 10 trials (ms):

| problem | base | final | ratio |
|---|---|---|---|
| BOYD1 | 92.17 | 86.88 | 0.943 |
| BOYD2 | 65.42 | 61.77 | 0.944 |
| CONT-300 | 49.73 | 44.25 | 0.890 |
| MAROS-R7 | 25.16 | 24.44 | 0.971 |
| CONT-201 | 20.48 | 19.01 | 0.928 |
| CONT-200 | 20.02 | 18.88 | 0.943 |
| PILOT87 | 19.84 | 18.00 | 0.907 |
| FIT2P | 18.69 | 18.48 | 0.989 |
| FIT2D | 16.80 | 15.94 | 0.949 |
| STOCFOR3 | 13.95 | 13.03 | 0.934 |
| PILOT | 13.87 | 12.86 | 0.927 |
| DFL001 | 13.77 | 12.81 | 0.931 |
| QAP15 | 13.28 | 12.73 | 0.959 |
| GREENBEA | 12.39 | 11.39 | 0.920 |
| GREENBEB | 11.70 | 10.90 | 0.932 |

Honest reading: a 5-6% corpus-level saving, 6-11% on the large bound-heavy
models (CONT-*, PILOT*, BOYD*), neutral where the time is in
`dual_propagation`, substitutions or sparsification (MAROS-R7, FIT2P). The
instrumented per-rule runs attribute the saving to `singleton_columns`
(-60 ms), `remove_redundant_bounds` and the parallel candidate sort; the
structural counts (7.67M -> 0.85M scan entries, 553k relaxations without queue
pushes, 1.95M sort keys of one word) do not depend on the load.

Instrumented per-rule runs before (`profile3`, baseline code) and after (final
code with the same instrumentation re-applied; both in-process size runs under
load, so unrelated rules move by +-10% and two stalls are visible: BOYD2's
redundant-bound pass and FIT2P's `dual_propagation` in the after run):

| rule | before ms (% of run) | after ms (% of run) |
|---|---|---|
| run | 971 | 1001 |
| propagate_bounds | 297 (31%) | 343 (34%) |
| remove_redundant_bounds | 176 (18%) | 172 (17%) (CONT-300 27.6 -> 16.9, BOYD1 26.8 -> 18.3; BOYD2 21.9 -> 60.0 is a stall) |
| singleton_columns | 131 (13%) | 50 (5%) |
| dual_propagation (untouched) | 111 (11%) | 156 (16%) (FIT2P 5.0 -> 20.7 stall) |
| sparsify_rows | 90 (9%) | 102 (10%) |
| parallel_columns | 72 (7%) | 71 (7%) |
| parallel_rows | 42 (4%) | 44 (4%) |
| build_model | 76 (8%) | 72 (7%) |
| cleanup | 24 (2%) | 30 (3%) |

The instrumentation itself (counters in `set_bounds`, `activity`, `equation`)
inflates `propagate_bounds` and `remove_redundant_bounds` in both columns
relative to the rule-level-only table in section 1.

What is left, per the profile: `propagate_bounds` (the tightenings and their
`set_bounds` row visits), `dual_propagation` (other engineer), `sparsify_rows`
(already ~4 ns per work unit), the parallel fingerprint/sort (64 + 42 ms, of
which only ~14 ms is repeat-phase work), and `build_model` (~7%, spread over
three passes).

## 5. Tests, lints, formatting

- `cargo test --release`: all passing (31 + 11 + 4 + 5 + 6 + 1 tests).
- `cargo clippy --all-targets --release`: no warnings.
- `cargo fmt --all --check`: clean.
- No README, public API or default `Settings` changes.

## 6. Ideas tried or evaluated and rejected

- **Caching fingerprints across medium phases** (invalidate per row/column on
  edit so repeated `parallel_rows/columns` scans only rehash changed items):
  measured that only 16% of fingerprinted items (estimated 14 ms of 64 ms) are
  in repeated phases, and the slowest problems (BOYD1/2, CONT-300) have a
  single medium phase. Not worth per-edit invalidation bookkeeping.
- **Restricting parallel scans to changed rows/columns**: a changed row can
  join a group of unchanged rows and the group base and comparison order must
  match a full scan, so every fingerprint is needed anyway; same 14 ms bound.
- **Pre-filtering parallel candidates by unique row length**: not exactly
  equivalent under a support-hash collision across different lengths (the
  collision bin's canonical ordering would change which pair is compared
  first), so rejected on the zero-difference requirement.
- **Skipping the `singleton_activity_rows` push in `set_bounds` when the row
  has no singleton column**: replaces one flag access by another (no saving),
  and the equivalence argument is subtle (a row can gain a singleton before the
  queue is drained). The consumption-side skip (2.1) is trivially equivalent.
- **Local copy of the row activity in the `propagate_bounds` entry loop**,
  refreshed only after a tightening: equivalent, but it saves one `Option`
  check and a 32-byte copy per entry, estimated below 5 ms for the 6.1M scanned
  entries; unmeasurable under the noise, so left out.
- **Skipping the residual recomputation when the cached sum equals the term**:
  not equivalent (the direct recomputation can yield a small nonzero residual
  that the cached subtraction rounded away), so it could change bounds.
- **Boxing large `Rule` variants** so the 336k `TightenedBound` records on
  BOYD1 push 40 instead of ~88 bytes: about 2-3 ms on BOYD1 and it touches
  postsolve matching code outside this area.
- **Radix sort on the packed candidate keys**: the sort is 42 ms corpus-wide,
  but the parallel path uses `par_sort_unstable` and the packed-key change
  already captures the cheap part.
- **Packing `seen/ratios/counts` in `sparsify_rows` into one struct**: the
  scan is already ~4 ns per work unit and dominated by linked column traversal
  on the problems where it matters (FIT2P, BOYD1, PILOT); the arrays fit in L1.
- **Reusing `Arc<Equation>` snapshots across propagation calls**: 2.35M copied
  entries corpus-wide (about 5 ms); would need a per-row cache invalidated in
  `changed_row`.
- **`build_model`**: measured its parts; no single fat target, and EXDATA's
  cost is the Hessian conversion in `matrix/sparse.rs`, outside this area.
- **Skipping the `rows[i]` load in `set_bounds` when there are no cones**:
  `rows[i]` is needed anyway for the doubleton re-queue check.

## 7. Files

- `optimization.patch`: `git diff HEAD` with the three changes (5 files,
  no instrumentation).
- Benchmark results: `benchmark/results/agent/size/{base,c1,c2,c3}.jsonl` and
  `benchmark/results/agent/time/*.jsonl`.
