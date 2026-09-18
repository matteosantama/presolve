# Efficiency review — September 18, 2026

Four reviewers each covered one part of the pipeline at commit 2528782, using
counters and phase timers in scratch copies of the tree, `size_of` checks, and
the release disassembly; no corpus timing runs were made while they overlapped.
Their reports are beside this file:

- `matrix-and-propagation.md`: the linked arena, activity caches, queues, and bound propagation
- `substitution-and-hessian.md`: substitution, the symmetric arena, sparsification, dependencies
- `discovery-and-scheduling.md`: parallel rows/columns, dominated columns, dual propagation, scheduler
- `boundaries-and-build.md`: CSC in/out, model construction, the tape, the release profile

Where the pipeline's 527 ms (default preset) goes, from their measurements:
entry and exit boundary about 13%, parallel detection about 13% (fingerprints
34 ms, sorts 27 ms), scheduler fixed costs under 0.1%, problems under 1 ms
only 10.7% of the sum. Everything below is bit-identical unless marked.

## Tier 1: low effort, no reduction risk

| Finding | Est. impact | Source |
| --- | ---: | --- |
| `#[inline]` on the hot cross-module helpers the 16-CGU build leaves out of line (`Model::activity` 12.0M calls, `Activity::replace_bound` 10.7M, `ActivityRows::push` 9.1M, `Worklist::push` 4.8M, `Queues::column_changed`, `Model::changed_row`, `Iter<0>::next`, `Queues::row_changed`) | 2–5% | matrix #1, boundaries #3 |
| `Model::from_parts`: walk the arena by index for locks and bulk-initialise activities, equations, and queues in the same push order | 1.5–2% | boundaries #1 |
| Candidate sort: MSD radix with constant-digit skip and an `is_sorted` pre-check (djb2 keys have a near-constant low word on ±1 matrices; two of four LSD passes order nothing) | 1.5–2% | discovery #1 |
| `implied_bound`: test "not tighter than the current bound" before the gain arithmetic (81% of 7.87M calls end there) | 0.5–2% | matrix #3 |
| `short_equalities`: load the column degree and row maximum after the structural rejections that drop 89% of candidates | 0.6–1.1% (CONT-300 4–7%) | substitution #1 |
| `from_columns`: one-entry lookahead for both column links, `u32` offsets | 0.5–1% | matrix #4, boundaries #4 |
| Skip a parallel row/column scan when the previous one was fruitless and the model's structural revision is unchanged (64 + 44 exact repeats of 404 scans) | 0.7% | discovery #2 |
| `try_substitute`: cursor walk over the affected Hessian row instead of a binary search per candidate pair (996k lookups) | 0.6–0.8% | substitution #3 |
| Lazy Hessian fingerprints: sort by the A key, hash P only inside equal-key runs | 0.6% (EXDATA −19%) | discovery #3 |
| `remove_variable` returns a `Vec` both callers discard (84k mallocs) | 0.4–0.6% | substitution #4 |
| `Activity::compute`: separate exclude/no-exclude instantiations, branch-free `Extreme::add` (never produces −0.0) | 0.3–1% | matrix #5 |
| Singleton round: drop the scattered `rows[i]` read (`row_singletons == 0` implies not deleted) | 0.2–0.4% | substitution #6 |
| `fix`: no temporary `Vec` of shifted domains (51k mallocs) | 0.2–0.5% | substitution #9, matrix #6 |
| `sparsify_rows`: allocate the `Scan` table lazily (17.8 MB memset) | 0.2–0.3% | substitution #8 |
| `size()` twice per call; `before` in O(1), `after` on the survivors walk | 0.2% | boundaries #5 |
| `coupled_dual_fix`: maintained list of Hessian-coupled columns instead of an O(n) scan per phase | 0.2–0.3% | discovery #5 |
| `groups_from` without `Vec<Vec<usize>>` (31.5k allocations); `singleton_rows` length check before the snapshot | 0.1–0.2% each | discovery #6, matrix #7 |

Estimates overlap (several target the same loops), so the realistic total for
this tier is 6–10% of the corpus rather than the column sum.

## Tier 2: medium effort, no reduction risk

| Finding | Est. impact | Source |
| --- | ---: | --- |
| Tape storage: `Rule` is 80 B while 85% of records need 32 B; each `TightenedBound` clones an `Arc<Equation>` (667k clone/drop pairs, 188k snapshots = 376k mallocs, mostly one per row). Snapshot arena with `(start, len)` handles, `u32` indices, packed rare variants. | 1.5–3% (BOYD1 −2–3%, CONT-200/201/100 5–12%) | boundaries #2, substitution #2 and #7, matrix cross-scope |
| One 64-byte row record (activity, `RowDomain`, kind, queue flags) instead of four scattered per-row arrays; `set_bounds` touches four cache lines per incident row today. First step: fold `row_kinds` into `ActivityRows.flags`. | 3–6% on large problems; first step 0.5–1% | matrix #2 |
| Exit path: seven separate O(n)/O(m) passes, identity `input_linear` map, `usize` maps in `pack` | 0.3–0.6% | boundaries #6, #7 |
| Nine separate queue flag arrays and a non-zero `stale()` fill at construction | 0.1–0.2% | boundaries #9 |

## Needs a reduction check before adoption

- `replace_row` re-queues 587k retained columns whose coefficient and locks
  did not change (UBH1, DEGEN3, GREENBEA about 10% each). Skipping them
  changes first-push order within a round; `run size` on both presets decides
  (substitution #5).
- `fingerprint` SIMD path: short-row problems lose 10–15% of fingerprint time
  to the SIMD build, long-row ones gain 4–8%; splitting the scalar short path
  out is arithmetic-identical but needs a quiet A/B (discovery #4).

## Aggressive preset only

- Hessian arena applies 13.5M staged updates one entry at a time (248M
  elements shifted, 744k relocations on AUG2DC/AUG2D/CVXQP1_L); a per-row
  sorted merge fixes it (substitution A1).
- `equality_dependencies` does 1.5M binary searches for 12.5k pivot hits; a
  stamp array (substitution A2).
- `short_equalities` re-examines rows on every bound change (4.16M
  examinations for 243k attempts; CONT-200 alone 2.3M) (substitution A3).
- Dual propagation internals (flat proofs, dense weights) are ≤ 0.1% (discovery #7).

## Build settings

No `[profile.release]` means the benchmark measures what a default user
build gets, and a library's profile is ignored downstream anyway, so the
portable fix is `#[inline]` on the hot leaves (tier 1). `codegen-units = 1`
removes the same out-of-line calls (text −10%) and is worth adopting only as
a measurement-stability setting, documented as such. `lto = "thin"` is a
no-op here, `panic = "abort"` is safe (no `catch_unwind`) but benchmark-only,
and `target-cpu=native` gains nothing on `aarch64-apple-darwin`, which already
targets `apple-m1`.

## Measured and rejected

Scheduler fixed costs (cleanup loop, `Instant::now`, round allocations:
< 0.1%); precomputing activities at build; merging double-side tightenings
(5.1k of 647k); arena compaction (85% of row steps are already contiguous);
column-sweep row fingerprints (2× slower); pre-sizing candidate vectors
(slower on macOS); better hash mixing (changes group order); the Hessian
copy (bandwidth-bound; only a structure-of-arrays arena would help);
`singleton_range` shortcuts (every cheaper form changes floating-point order).

## Suggested order

1. The trivial tier-1 items in one commit each where they touch different
   code: inlining, `implied_bound` early exit, `short_equalities` ordering,
   the four allocation removals, the singleton-round load.
2. `Model::from_parts` bulk initialisation and the `from_columns` lookahead.
3. The MSD radix sort, repeated-scan skip, and lazy Hessian fingerprints.
4. Tape storage, then the row record.

Each step: identical reductions on both presets, then a matched fresh-process
timing pair on a quiet machine.
