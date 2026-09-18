# Efficiency review: discovery rules, scheduler, executor

Branch `new-presolve-rules` at 2528782, September 18, 2026. Scope: `src/rules/mod.rs`,
`src/rules/parallel.rs`, `src/rules/dominated_columns.rs`, `src/rules/dual_propagation.rs`,
`src/rules/dual_fixing.rs`, `src/rules/cones.rs`, `src/executor.rs`, `src/model/queues.rs`.

No corpus timing runs were made (other agents were benchmarking). Every number below comes
from one of three load-independent or load-robust sources: (a) the existing
`benchmark/results/plain-csc-20260918/time/new-*.jsonl` runs (corpus sum of per-problem
medians, min over the four runs: **527 ms**); (b) an instrumented copy of the tree in my
scratchpad (`git archive 2528782`, counters in every in-scope function, run over all 236
problems under both presets: counts are deterministic); (c) microbenchmarks taken as the
minimum of 9 interleaved repetitions on the real matrices or the real fingerprint keys, so
ratios are trustworthy even though the machine was loaded (load average 5). Absolute
microbenchmark times are probably 10-20% high.

## Summary

| Rank | Finding | Est. corpus impact (default preset) | Risk to reductions | Effort |
| --- | --- | --- | --- | --- |
| 1 | Candidate sort: replace `RADIX_THRESHOLD` + 4-pass LSD radix with an MSD radix (adaptive digit, constant-digit skip, `is_sorted` pre-check, LSD fallback on skewed keys) | −1.5% to −2% (sort work −31% measured on all 472 real key sets; −3.4% BOYD2, −5% CONT-300, −6% CONT-200/201) | none: result is full `Ord` order by construction, existing unit test covers it | 2-3 h |
| 2 | Skip a parallel scan when the previous scan of the same kind was fruitless and its inputs are unchanged (rows: structural revision bumped in `changed_row`; columns: `revision`) | −0.7% default (64 of 404 row scans, 44 of 404 column scans are exact repeats); about −1.2% of the 6.4 s aggressive corpus (152 + 195 repeats) | none for reductions; `Stats::parallel_comparisons` no longer counts the repeated comparisons | 1 h + audit (done below) |
| 3 | Hash Hessian columns lazily: sort by the A key first, hash P only inside runs of equal A key | −0.6% corpus, −19% on EXDATA (1.8 ms of 9.6); column sort elements shrink 32 → 16 B | none: order is provably identical | 1-2 h |
| 4 | `fingerprint`: keep the scalar short-row path out of the SIMD function body (or drop `wide`) | −0.3% to −0.5%, uncertain; short-row problems lose 10-15% of fingerprint time to the SIMD build, long-row problems gain 4-8% from it | none (same arithmetic; tests assert equality) | 30 min + quiet A/B |
| 5 | `coupled_dual_fix`: iterate a maintained list of columns with 2..=33 Hessian entries instead of scanning all `n` every medium phase | −0.2% to −0.3% (878k column visits, QPs only) | none if the list stays in index order | 1 h |
| 6 | `groups_from`: iterate runs of `entries` in place instead of allocating `Vec<Vec<usize>>` | −0.1% to −0.2% (31.5k allocations) | none | 30 min |
| 7 | Dual propagation internals (aggressive only): flat CSR proofs instead of `Vec<Vec<u32>>`, dense weights + heap instead of `BTreeMap` | ≤ 0.1% of the aggressive corpus | none | 2 h, low value |
| — | Scheduler fixed costs (cleanup loop, `Instant::now`, round allocations, `work_size`) | measured negligible (< 0.1%) | — | do nothing |

Findings 1-3 together are worth roughly 3% of the default corpus and are all
bit-identical by construction. Findings 4-6 are small and cheap. Nothing in the scheduler's
per-phase bookkeeping is worth touching. Several plausible ideas were measured and rejected
(section "Tried and rejected"), including the column-sweep row fingerprint (2x slower) and
pre-sizing the candidate vector (slower on macOS).

## Where the time is (in scope)

Instrumented counts over the corpus, default preset, and per-unit costs from the
microbenchmarks (fingerprint on the linked matrix 3.3 ns/entry, sort 13.9 ns/item,
`proportional` about 50 ns/comparison):

| Component | Count (default preset) | Estimated cost | Share of 527 ms |
| --- | --- | ---: | ---: |
| Row + column A fingerprints | 10.2 M entries (808k rows, 1.15 M columns) | 34 ms | 6.5% |
| Hessian column fingerprints | 3.73 M entries (531k columns; EXDATA 2.25 M) | 3.4 ms | 0.6% |
| Candidate sorts | 808 sorts, 1.95 M items (9 radix sorts, 712k items) | 27 ms | 5% |
| `proportional` comparisons | 98.5k (FIT2P 27k, 80BAU3B 22k, FIT2D 14k) | 5 ms | 1% |
| `coupled_dual_fix` scans | 205 scans, 878k column visits | 1-2 ms | 0.3% |
| Cleanup loop | 4196 calls, 5087 iterations, 3675 first-iteration no-ops | 0.2 ms | 0.04% |
| `Worklist::take_round` / `ActivityRows` drains with allocation | 5189 | 0.25 ms | 0.05% |
| `Instant::now` in `run_phases` | 930 | 0.02 ms | 0 |

Parallel detection is about 13% of the corpus, consistent with the 12.8% attributed to it in
`docs/benchmarks/new-rules-20260917/schedule-REPORT.md`. The five problems with n or m
≥ 32768 (BOYD1/2, CONT-200/201/300) are 36% of the corpus; problems under 1 ms are 10.7%
and under 0.25 ms only 1.2%, so per-phase fixed costs cannot move the sum (they move the
geometric mean, which the owner also tracks).

Phase counts: 85 problems run one medium phase, 136 two, 13 three, 2 four. Every medium
phase re-fingerprints and re-sorts every row and column.

## Finding 1: the candidate sort

**Where.** `src/executor.rs:62-110` (`RADIX_THRESHOLD = 32_768`, `sort_by_radix_key`),
called from `src/rules/parallel.rs:119-131` (`sorted_candidates`).

**What happens now.** Below 32768 candidates the vector is sorted with `sort_unstable`
(pdqsort/ipnsort, 24 B elements for rows, 32 B for columns); above it, four 16-bit LSD
counting passes over `(key, position)` pairs, each allocating a 512 KB `starts` array and
scanning it for a constant digit, then a gather and a tie-break sort of equal-key runs.

**Why it costs.** The packed key is `support_hash << 32 | coefficient_hash`, both djb2. On
the structured problems the coefficient hash is nearly constant (all coefficients ±1 give
the same quantized values) and the support hash of short rows has nearly constant top bits
(`5381 * 33^k + Σ j * 33^…` stays small for k ≤ 3). Distinct values per 16-bit window
(bits 0-15, 16-31, 32-47, 48-63) on the initial matrices:

| Keys | n | distinct per window |
| --- | ---: | --- |
| BOYD2 rows | 186531 | 10, 10, 65536, 53 |
| BOYD2 columns | 93263 | 698, 109, 52642, 57 |
| CONT-300 rows | 90298 | 7, 6, 57563, 57454 |
| LISWET1 rows | 10000 | 1, 1, 10000, 173 (and the keys arrive already sorted) |
| BOYD1 columns | 93261 | 48444, 40457, 249, 29 (discrimination only in the low word) |
| FIT2P columns | 13525 | 25, 25, 3025, 26 |

So two of the four LSD passes order almost nothing on BOYD2/CONT-* (the counter
`sort_radix_passes` confirms all 8 passes ran for BOYD2's two sorts: the constant-digit
skip never fires because 10 distinct values is not 1), and the 65536-bucket scatter is
cache-hostile at 40k-186k items (each bucket receives 1-3 elements).

**Measured.** Standalone microbenchmark on the packed keys of all 236 problems (472 key
sets, one initial scan each; min of 9; production semantics for "current"):

| Variant | Total over 472 key sets | n ≥ 32768 subset | Files > 5% slower |
| --- | ---: | ---: | ---: |
| current (`sort_unstable` < 32768, else 4×16 LSD) | 24.1 ms | 12.3 ms | — |
| MSD radix, density 4, constant-digit skip, `is_sorted` pre-check, n ≥ 4096 | 17.3 ms (−28%) | 7.6 ms | 3 (BOYD1 cols +0.8 ms, FIT2P cols +57 µs, DTOC3 cols +19 µs) |
| same, falling back to the current LSD when `distinct(bits 32-47) < n/64` | 16.6 ms (−31%) | 7.2 ms | 2 (FIT2P cols +57 µs, DTOC3 cols +19 µs) |

Per-problem (one scan): BOYD2 rows 3.05 → 1.82 ms, BOYD2 cols 1.57 → 1.21, CONT-300 rows
1.54 → 0.51, CONT-300 cols 1.55 → 0.66, CONT-200/201 rows and cols 0.73-0.75 → 0.19-0.28,
STOCFOR3 cols 0.30 → 0.19, CONT-100 cols 0.17 → 0.07, LISWET rows 7 → 4.5 µs (sorted
pre-check). The synthetic-key benchmark (uniform random keys) showed the same ordering:
MSD 2x faster than `sort_unstable` at 2k-40k and 1.6x faster than the LSD at 93k-186k.

**Design.** `msd_radix(values, scratch, consumed_bits)`: if `len ≤ 64` or all 64 bits are
consumed, `sort_unstable`; else pick `bits = clamp(log2(len/4), 4, 16)`, histogram the next
`bits` most significant unconsumed bits (u32 counters), if one bucket holds everything
consume the digit and retry without scattering, otherwise scatter into `scratch`, copy back,
recurse into buckets longer than one. Before that, `if values.is_sorted() { return }` (44 of
the 472 key sets, 175k items, arrive sorted: LISWET, UBH1, AUG2D rows and similar) and use
`sort_unstable` below 4096. Fallback: the level-1 histogram and a second 16-bit histogram
of bits 32-47 are computed in the same pass; if `distinct(bits 32-47) < n/64` the keys are
BOYD1-like (few supports, discrimination in the coefficient word) and the existing LSD is
used. One scratch `Vec<T>` replaces the current three (`keyed`, `scratch`, `sorted`) and
the four 512 KB `starts` allocations.

**Expected effect.** Sort work −31% on the initial scans; scaled to the 1.95 M items
actually sorted, about 8-9 ms, i.e. 1.5-2% of the corpus, concentrated on BOYD2 (−1.6 ms,
3.4%), CONT-300 (−1.9 ms, 5.4%), CONT-200/201 (−1 ms each, 6%), STOCFOR3, DFL001, QAP15,
AUG2D*, CONT-100/101. BOYD1 unchanged (fallback). Residual regressions to verify on a
quiet machine: FIT2P columns +57 µs per medium phase (2 phases, 12 ms problem: +0.9%) and
DTOC3 +19 µs. A variant with density 16 and a "poorly spread digit → `sort_unstable`"
guard makes FIT2P cols 113 µs (faster than today) but costs elsewhere; if FIT2P must not
regress, apply that guard only when the level-2 histogram has fewer than `len/64`
non-empty buckets and re-measure.

**Risk to bit-identical reductions.** None. The output is the unique full-`Ord` order of
distinct values (index is part of every value), which is what `sort_unstable` produces;
`radix_ordering_matches_comparison_ordering_with_duplicate_prefixes` already asserts this
across the threshold and should be extended with a pre-sorted input and a BOYD1-like key
distribution. The parallel executor path (`par_sort_unstable`) is untouched.

**Effort.** 2-3 hours including tests. **Measurement that settles it:** quiet A/B on
BOYD2, CONT-300, CONT-201, STOCFOR3, FIT2P, DTOC3, LISWET1.

## Finding 2: skip exact repeat scans

**Where.** `src/rules/mod.rs:130-135` (the medium phase calls `parallel_rows` and
`parallel_columns` unconditionally).

**What happens now.** Under the default `Progress::Nonzeros` policy a cycle that reduced
nonzeros by 5% is followed by another cycle whose fast phases and propagation may change
nothing, after which both scans run again on an unchanged model. Instrumented counts
(default preset): of 404 `parallel_rows` calls, **37** start at the very same `revision` at
which the previous `parallel_rows` call started and that call changed nothing; with a
structural revision that ignores variable-bound changes (see below) it is **64** calls
(282k entries fingerprinted, 58.8k items sorted). Of 404 `parallel_columns` calls, **44**
are exact repeats (346k A+P entries, 70k items). Under the aggressive preset (`AnyChange`)
it is 152 of 600 row scans (3.9 M entries) and 195 of 600 column scans (13.3 M entries).

Problems affected (default): STOCFOR3, FIT2P, WOOD1P, CVXQP1/2/3_L, FIT1P, SHIP04L/08L,
QSHIP04L/08L, STOCFOR2, and about 40 smaller ones; each pays one useless scan of 0.1-0.5 ms
(3-5% of those problems).

**The change.** Record, per scan kind, the input revision at entry and whether the scan
changed the model (`revision` unchanged at exit). At the next call, return `Ok(0)` when the
previous scan was fruitless and the input revision is unchanged. For columns use
`self.revision` (the column scan reads bounds, `c`, `P`, `alive`, so every bump matters).
For rows add `rows_revision`, incremented in `Model::changed_row` (`src/model/mod.rs:271`),
because `parallel_rows`/`merge_parallel_rows` read only `rows[i]`, `a.row(i)` and settings,
and variable-bound tightenings (the most common edit between phases) never touch those.

Audit that every row-domain or matrix mutation ends in `changed_row`: writes to
`self.rows[..]` at `model.rs:238` (`replace_rows_batch`, `changed_row` at 252), `:520`
(`fix`, 527), `:743` (`replace_row_bounds`, 745) and `cones.rs:257` (cone→cone block
rewrite, irrelevant to linear rows); matrix mutations `replace_rows` (242 → 252),
`replace_row_into` (309 → 358), `remove_column` in `fix` (516 → 527) and in `aggregate`
(818 → 820); `add_column` (217) adds an empty column and changes no row. So `rows_revision`
is exact for the row scan's inputs.

**Why it is exact.** A scan is a deterministic function of the model state it reads
(serial executor; the parallel executor sorts the same values). A fruitless scan at state
S followed by a scan at the same S returns the same empty result. Only
`Stats::parallel_comparisons` changes (it no longer counts the repeated comparisons, e.g.
FIT2P's 13.5k per phase); the harness does not record that statistic.

**Expected effect.** Default: about 1.75 ms (rows) + 2.1 ms (columns) ≈ 3.9 ms ≈ 0.7% of
the corpus, spread over ~50 mid-size problems at 3-5% each. Aggressive: about 60-75 ms,
1% of the 6.4 s aggressive corpus. **Effort:** one hour. **Measurement:** quiet A/B on
STOCFOR3, FIT2P, WOOD1P, CVXQP3_L.

Not recommended beyond this: a per-row/column fingerprint cache with invalidation in
`changed_row`/`column_changed`. The repeat-phase fingerprint work is bounded by 1.3 M
entries (fp entries minus initial nnz on the 151 multi-phase problems) ≈ 4.4 ms, of which
this finding already removes about half; the remainder needs a byte store per edit and two
cache arrays. That agrees with the rejection in the schedule report.

## Finding 3: lazy Hessian fingerprints

**Where.** `src/rules/parallel.rs:285-300` (`parallel_columns` key closure hashes
`self.objective.p.column(j)` for every column with Hessian entries), consumed at
`dominated_columns.rs:245` (`key.1.is_some()`).

**What happens now.** 531k Hessian columns, 3.73 M entries, are hashed per medium phase in
the default preset; EXDATA alone hashes 1500 columns of 1500 entries (2.25 M) in its single
medium phase, Q25FV47 211k, QSHIP12L 176k, QSHIP08L 136k, CVXQP3_L 110k. The contiguous
slice hash costs 0.8 ns/entry on EXDATA (measured, SIMD path) and 3.2-3.4 ns/entry on
diagonal Hessians (one entry per column, per-call overhead), so EXDATA spends about 1.8 ms
of its 9.6 ms here and the corpus about 3.4 ms.

**The change.** Sort by `(a_key, index)` first. Then, inside each run of equal `a_key`
with at least two members, compute `p_key = (!p.is_empty()).then(|| packed(fingerprint(p)))`
and `sort_unstable` the run by `(p_key, index)`. `groups_from` and
`dominated_support_groups` then operate on `(a_key, p_key, index)` as today;
`dominated_support_groups` only needs "P non-empty", which is an O(1) slot read.

**Why it is order-identical.** The current order is lexicographic on
`(a_key, p_key, index)` with `None < Some`. Sorting by `(a_key, index)` and then each
equal-`a_key` run by `(p_key, index)` yields exactly that order; runs of length one need no
`p_key` because their group membership cannot depend on it. On the corpus the runs are
almost empty: EXDATA, QSHIP12L, QSHIP08L, CVXQP1/3_L have **zero** columns inside runs of
equal A key, Q25FV47 has 18, BOYD1 6276 (of 93k, all single-entry P). Essentially every
Hessian hash disappears.

**Bonus.** The first sort's element becomes `(u64, usize)` (16 B) instead of
`((u64, Option<u64>), usize)` (32 B); the column sorts then run at the row-key speed
(BOYD2 cols 1.21 → about 1.0 ms with finding 1, CONT-300 cols 0.66 → about 0.5).

**Expected effect.** −0.6% corpus, −19% EXDATA, −2-4% on the QSHIP/CVXQP/Q25FV47 family.
**Risk:** none. **Effort:** 1-2 hours (restructure `parallel_columns`, keep
`dominated_support_groups` signature).

## Finding 4: the SIMD path and the short-row path

**Where.** `src/rules/parallel.rs:20-70` (`fingerprint`).

**What happens now.** Rows and columns shorter than 8 take the scalar fold; longer ones
take the `f64x4`/`u32x4` blocks. Counts: 9% of row items but 50% of row entries, 7% of
column items and 29% of column entries go through the SIMD blocks.

**Measured** (all rows and columns of every initial matrix, alternating order, min of 9;
identical results asserted): SIMD build 31.0 ms, scalar-only build 29.5 ms. Short-row
structured problems are 10-15% slower with the SIMD build (BOYD2 3.26 vs 2.97 ms, CONT-300
2.77 vs 2.41, CONT-200/201 1.09/1.13 vs 0.92/0.95); long-row problems are 3-10% faster
with it (MAROS-R7 0.92 vs 1.00, PILOT87 0.60 vs 0.66, WOOD1P 0.56 vs 0.59, FIT2D 0.85 vs
0.88). Since both builds run the identical scalar code on short rows, the slowdown is
codegen: the scalar fold lives in a function that also carries the SIMD constants and
blocks, and inlines worse. On the linked matrix the block gather is still a pointer chase,
so the vector arithmetic saves little; on contiguous Hessian slices it helps (0.8 ns/entry
on EXDATA), but finding 3 removes those hashes anyway.

**The change.** Either (a) split `fingerprint` into an `#[inline]` scalar function for
`len < 8` and a separate `#[inline(never)]` SIMD function for long items, so the hot
short path compiles as in the scalar-only build; or (b) drop the SIMD path and the `wide`
dependency. (a) keeps the long-row gain; (b) costs MAROS-R7 about +75 µs per phase.

**Expected effect.** −0.3% to −0.5% for (a) (about 1.5 ms per full scan set on the
structured problems: BOYD2 −0.3 ms, CONT-300 −0.35 ms), uncertain until A/B'd; no change on
long-row problems. **Risk:** none (`fingerprints_preserve_scalar_rounding_tails_and_overflow`
asserts equality). **Effort:** 30 minutes plus a quiet A/B on BOYD2, CONT-300, MAROS-R7.

## Finding 5: `coupled_dual_fix` scans every column

**Where.** `src/rules/dual_fixing.rs:73-120`, called every medium phase from `mod.rs:123`.

**What happens now.** For any problem with Hessian entries, all `n` columns are visited
and `self.objective.p.row(j).len()` read (a 24 B slot) to find columns with 2..=33 Hessian
entries: 205 scans, 878k visits in the default preset (BOYD1 93k, CONT-300 91k, CONT-200/201
40k each, AUG2D* 20k…). On diagonal Hessians (BOYD1, CONT-*, LISWET, AUG2D*) no column
qualifies, so the scan does nothing.

**The change.** Keep a `Vec<usize>` of columns whose Hessian row length is in 2..=33,
rebuilt (in index order) when `objective.p.revision` differs from the value at the last
build, and iterate that list. Same columns in the same order, so the same `fix` calls.

**Expected effect.** 1-2 ms, 0.2-0.3% of the corpus (BOYD1 −0.2 ms, CONT-300 −0.2 ms).
**Risk:** none. **Effort:** one hour.

## Finding 6: per-group allocations

**Where.** `src/rules/parallel.rs:134-148` (`groups_from`), used at `:162` and `:301`.

**What happens now.** Every group of ≥ 2 equal keys becomes a `Vec<usize>` inside a
`Vec<Vec<usize>>`: 171 row groups and 31 339 column groups per corpus run (FIT2P 6000,
TRUSS 4399, FIT2D 2935, 80BAU3B 1405 per phase). About 1 ms of allocation.

**The change.** Yield `(start, end)` runs over the sorted `entries` slice, as
`dominated_support_groups` already does, and index `entries[at].1`. Group order and member
order are unchanged. **Effect:** 0.1-0.2%. **Risk:** none. **Effort:** 30 minutes.

## Finding 7: dual propagation internals (aggressive preset only)

**Where.** `src/rules/dual_propagation.rs:33-57` (`DualScratch`), `:141` (`proofs`
resized to `2m` empty `Vec`s), `:227-233` (a `Vec<u32>` per tightened `(row, side)`),
`:364-396` (`BTreeMap` walk in `apply_direction`).

**Counts** (aggressive, all 235 calls): 72 early returns (all-curvature QPs); 382k column
visits over 2.55 M entries in 563 rounds; 352k queue pushes; 39k records into 21k proof
slots (21k inner-`Vec` allocations); 1915 directions tried, 22.7k records walked, 46k
`BTreeMap` operations, 51k verification entries, 1856 reductions, 59 rejections (all in the
row check). MAROS-R7: 855 directions, 20k records walked, 40k map operations.

The extraction/verification machinery is small (≈ 100k entry-level operations corpus-wide);
the rule's cost is the propagation sweep itself plus the two conclusion sweeps (282k row
and 233k column candidates, 75k `Activity::compute` recomputations for boxed columns, which
are never queued). Changing `proofs` to a flat CSR built after propagation (one allocation,
same `partition_point` lookups) and `weights` to a dense `Vec<f64>` indexed by sequence
number plus a max-heap of touched sequences (identical accumulation order, identical pop
order) would save perhaps 3-4 ms of the 6.4 s aggressive corpus. `DualScratch` is created
per presolve call, so its reuse only helps if the rule ran more than once per call, which it
no longer does. Low priority; listed for completeness.

`dominated_support_groups` (aggressive only): 77k runs, 418k pairs, 642k equality
fingerprints, 1.2 M merge steps, 1208 fixes; `begin_pass` resizes two `n`-length arrays per
medium phase (negligible). The general search (`dominated_columns`, 5.7 M visits) is the
known 7% and is off by default.

## Scheduler and executor: measured and negligible

- `cleanup()` (`mod.rs:39-68`): 4196 calls, 5087 iterations, 3675 first iterations that
  change nothing (six empty-queue pops each). About 0.2 ms.
- `run_phases` `Instant::now()` (`mod.rs:89,174,181,198`): 930 calls per corpus run.
- `Worklist::take_round` and `ActivityRows` drains that allocate: 5189 per corpus run
  (about 0.25 ms). `swap_round` already covers the round-based rule.
- `work_size()` is O(1); `significant_progress` is a few flops.
- `propagate_rounds` (`mod.rs:206-244`): 106 extra rounds in the default preset; the cost
  estimate loop reads one list header per queued row.
- The aggressive `AnyChange` re-queue of equalities (`mod.rs:153-165`) walks every row per
  cycle: 1.56 M row visits, about 3 ms of the 6.4 s aggressive corpus.
- `Executor::filter_map_sorted` serial path: no indirection cost (generic, monomorphized).
  `enough_work` is computed and ignored serially (trivial).
- No phase runs with nothing queued and still pays setup, except the two parallel scans
  (finding 2) and `coupled_dual_fix` (finding 5).

## Tried and rejected

- **Row fingerprints by a column sweep** (accumulate each row's djb2 state while walking
  columns in order, so row traversal of a column-major node arena is avoided). Implemented
  and verified bit-identical on all 236 matrices; it is **2x slower** (12.4 ms → 23.2 ms
  corpus-wide; BOYD2 1.22 → 2.47 ms). Row traversal is only 6.5 ns per row on BOYD2: the
  arena is built column-major but the problems are banded, so a row's nodes are close in
  memory, and the second pass for `max` is L1-resident. The scattered per-row state writes
  of the sweep cost more.
- **`Vec::with_capacity(count)` for the candidate vector** (`executor.rs:51`). Measured
  slower (186k × 16 B: 131 µs grown vs 174 µs pre-sized; 40k: 31 vs 61 µs): a fresh multi-MB
  allocation page-faults on first touch, whereas macOS grows large blocks cheaply. Leave
  `collect()` alone.
- **Better hash mixing** (avalanche the djb2 words so the top bits spread and any radix
  works). Changes the sorted order of groups, hence tape order and queue push order: not
  bit-identical. The sort must adapt to the keys, not the other way round.
- **Memoizing failed `proportional` pairs across phases** (FIT2P re-compares 13.5k pairs
  per phase for nothing). The pairs are 98.5k comparisons ≈ 5 ms corpus-wide; finding 2
  removes the fully repeated phases, and per-pair state would cost more than the rest.
- **Per-item fingerprint cache** (see finding 2).
- **Skipping the dual-propagation sweep on problems that cannot tighten anything**: a free
  singleton column already fixes a multiplier, so the cheap pre-checks do not exist without
  per-entry work; the sweep is inherent.

## Method and files

- Instrumented copy: scratchpad `repo/` (from `git archive 2528782`), counters in
  `src/probe.rs`, probe binary `benchmark/src/bin/probe.rs` with modes `counts`, `fp`,
  `simd`, `pfp`, `keys`. Nothing under the project tree was modified.
- Corpus timing distribution: `benchmark/results/plain-csc-20260918/time/new-*.jsonl`.
- Sort benchmarks: scratchpad `sortbench{,2,3,4,5}.rs` over `keys/*.txt` (472 files).
- Aggressive corpus total for context: 6.4 s
  (`benchmark/results/refactor-20260918/time/t1-final-aggr.jsonl`).
- All microbenchmarks: min of 9 interleaved repetitions, load average about 5.
