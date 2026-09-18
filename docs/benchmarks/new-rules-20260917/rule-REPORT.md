# Dual propagation efficiency pass

Scope: `src/core/rules/dual_propagation.rs` and its scheduling. `src/core/schedule.rs`
is unchanged; the "skip when nothing changed" idea lives inside the rule (change 1),
which covers every caller and keeps the schedule readable.

All code changes are in `optimization.patch` (`git diff HEAD -- src`, 224 lines, two
files: `src/core/rules/dual_propagation.rs`, `src/core/queues.rs`).

## Where the time went (baseline)

An env-guarded instrumentation build (not part of the patch) timed every segment of
each `dual_propagation` call during the timed presolve of all 236 problems
(`--rule all`, reused pool):

| segment | corpus total | share of rule |
| --- | --- | --- |
| `effective_bounds` (activity copy, sort, two sweeps) | 77.4 ms | 58% |
| propagation loop | 39.2 ms | 29% |
| conclusion sweeps (rows + columns) | 17.5 ms | 13% |
| rule total | 134.6 ms | 16.9% of the 794.5 ms presolve total |

Counts: 406 calls; 1.51 M queue pushes for 82 k multiplier tightenings (MAROS-R7:
392 k pushes for 8.6 k tightenings, dense rows re-pushed at every tightening);
706 k `Activity::compute` calls in the loop plus one per dual column in the
conclusion sweep; 22 calls (3.6 ms) were repeats at an unchanged revision after a
fruitless pass; the work limit was hit 7 times. Inside `effective_bounds` the
`sort_unstable_by_key` on `(len, j)` cost 0.1-1.9 ms per call on the large problems
(BOYD2 1.9 ms, CONT-300 1.6 ms, CONT-201 0.66 ms, QAP15 0.45 ms).

## Changes and why the reductions are identical

1. Skip a pass at a settled revision (`DualScratch::settled`). After a pass that
   changed nothing (`self.revision` unchanged at exit) the revision is remembered; a
   later call at that revision returns 0 immediately. The pass is a pure function of
   the model state (rows, bounds, `alive`, coefficients, `c`, `P` columns, cones) and
   every mutation of that state bumps `revision` (audited every `self.rows[..] =`,
   `alive[..] = false`, `objective.c[..] =` site and the cone-block rewrite in
   `cones.rs`, which bumps revision at the end of its transaction). Saves the 22
   repeated calls (STOCFOR3 and FIT2P ran the identical pass twice).

2. Counting sort instead of `sort_unstable_by_key` for the visiting order. Columns
   are bucketed by length while scanning `j` upward, which yields exactly the
   `(len, j)` order the comparison sort produced (unique keys, so stability is
   irrelevant). Same order, same dropped sides. Removes about 5 ms of comparison
   sorting over the corpus.

3. Queue the columns of a tightened row once per round, not once per tightening.
   A row that tightened several times in one round walked its whole row each time
   (`Worklist::push` deduplicates, but the linked-list traversal was still paid).
   Tightened rows now go into a small row `Worklist`; their columns are pushed when
   the round ends. Next-round order is unchanged: in the baseline a column's queue
   position is set by its first push, i.e. by the first tightening of the first row
   (in tightening order) containing it; the row worklist keeps first-tightening
   order and pushes those rows' columns in the same row order, so the deduplicated
   column sequence is identical. Removes most of the 1.5 M row walks.

4. Reuse the round vector (`Worklist::swap_round`). `take_round` allocated a fresh
   `Vec` per round (up to 10 rounds per pass); the scratch now keeps one vector and
   swaps it with the queue. Added with a unit test in `queues.rs`.

5. Reuse the loop's column activity in the conclusion sweep. The last visit of
   column `j` computed `Activity::compute(column, y)`; when the queue drained
   normally every later tightening re-queued and re-visited `j`, so the stored value
   is over the final `y`. The conclusion sweep changes neither `y` nor the column's
   entries (`fix` on other columns and `implied_equality` touch other columns and
   row sides only), so the stored value equals what the baseline recomputed. When
   the work limit ended the loop early (`complete == false`) or the column was never
   visited, the sweep recomputes as before. Saves one `Activity::compute` per dual
   column (BOYD2: 93 k, CONT-300: 33 k).

Nothing changes the arithmetic: the same activities, residuals, implication and
tightening tests run in the same order on the same values, so the proved signs and
the resulting `implied_equality` / `fix` calls are the same.

## Verification

- `cargo run --release -p benchmark -- compare size base final --results-dir benchmark/results/agent`:
  `236 matched cases; 0 changed.` `0 unmatched or invalid cases.`
  (also 0 changed for the intermediate candidates `cand1` and `cand2`).
- `cargo test --release`: all 6 test binaries pass (31 + 11 + 4 + 5 + 6 + 1 tests,
  0 failed), including the new `swapped_rounds_reuse_the_previous_round_and_unmark_entries`.
- `cargo clippy --all-targets --release`: no warnings.
- `cargo fmt --all --check`: clean.
- No README, public API or default-settings changes; SPDX headers intact.

## Measured timings

Caveat: two other engineers were benchmarking on this machine (load average 15-32
during the last runs); two identical baseline runs of the whole corpus differ by
about 5%. Baseline and candidate were always interleaved (A/B/A), with the baseline
built from `HEAD` in a separate worktree. Values are sums of per-problem medians in
ms (5 trials, reused pool).

Whole corpus, `--rule all`:

| run | base | final | base again |
| --- | --- | --- | --- |
| first A/B/A (cand1 = final library code) | 747.5 | 711.5 | 727.1 |
| second A/B/A | 726.0 | 731.2 | 766.7 |
| initial isolated baseline (`base-time`, taken first) | 917.2 | | |

The whole-corpus effect (about -1% to -3%) is inside this machine's noise. The rule
is about 17% of presolve time and the change removes roughly a third of it, so a
3-6% total improvement is the structural expectation.

Rule only, `--rule dual_propagation` (rule plus cleanup), median ms, 7 trials,
interleaved base / final / base on the problems where the rule is heaviest:

| problem | base | final | base again | change |
| --- | --- | --- | --- | --- |
| BOYD2 | 15.48 | 14.36 | 15.21 | -6% |
| CONT-300 | 15.81 | 15.53 | 16.00 | -2% |
| D6CUBE | 1.16 | 0.93 | 1.19 | -20% |
| DFL001 | 2.57 | 2.30 | 2.67 | -12% |
| FIT2D | 4.19 | 4.10 | 4.23 | -3% |
| FIT2P | 3.23 | 2.69 | 3.31 | -18% |
| MAROS-R7 | 11.01 | 9.78 | 11.11 | -12% |
| PILOT87 | 2.25 | 1.94 | 2.23 | -14% |
| QAP15 | 2.75 | 2.56 | 2.94 | -10% |
| STOCFOR3 | 3.32 | 3.15 | 3.49 | -7% |
| WOOD1P | 2.33 | 1.98 | 2.25 | -14% |

Whole corpus, `--rule dual_propagation`, A/B/A/B at load average 16-32:
base 171.1, final 161.4, base 171.6, final 175.5 (the last run coincided with the
load peak; an earlier attempt at load 30+ was discarded because the baseline itself
doubled between its two runs).

Focus problems, `--rule all`, first A/B/A (base / final / base): MAROS-R7 33.9 /
23.2 / 25.1, DFL001 17.2 / 13.5 / 14.0, FIT2P 19.0 / 15.9 / 19.8, D2Q06C 11.7 / 10.5
/ 11.1, STOCFOR3 14.3 / 13.7 / 15.4, 80BAU3B 6.17 / 6.13 / 6.14, MODSZK1 1.35 / 1.27
/ 1.33, PILOT87 20.6 / 20.1 / 20.8, BOYD1 104.5 / 94.2 / 94.3 (BOYD1 has 18 rows;
its swing is scheduler noise).

## Tried and rejected

- Per-row "usable side" flags in `effective_bounds` (skip rows whose needed side is
  infinite or whose extreme already has two unbounded contributions before loading
  the row). Correct and size-identical (`cand2`, 0 changed), but the whole-corpus
  A/B/A was flat (707.2 / 707.3 / 715.7 ms for cand1 / cand2 / cand1) and rule-only
  timings split evenly (BOYD2, CONT-300, STOCFOR3 slightly worse; D6CUBE, FIT2D,
  MAROS-R7 slightly better). The extra per-row byte load costs about what the
  skipped branch saves, so it was dropped to keep the code simpler.
- Skipping `effective_bounds` when bounds and columns did not change since the last
  pass. In the corpus every call at a changed revision followed bound propagation or
  coupled dual fixing, which edit exactly that state, so the check would never fire;
  the unchanged case is covered by change 1.
- Replacing the activities copy with the model's cached activities. Not possible
  without changing results: the sweeps update the copy as sides are dropped so later
  implications never rest on a dropped side, while the model cache must keep the
  real bounds.
- Caching per-column activities incrementally in the loop. A column is recomputed
  only when queued, and it is queued only when a multiplier in it changed, so an
  incremental update touches the same entries; direct recomputation is also what
  avoids cancellation drift, and changing it could change which tightenings pass the
  relative-gain filter.
- Tracking strictly signed multipliers to shrink the row conclusion sweep. That
  sweep is a linear scan over `rows` / `y` with no matrix access and measured only
  0.3-1.4 ms on the largest problems; the column sweep was the expensive part and
  is covered by change 5 without extra bookkeeping.
