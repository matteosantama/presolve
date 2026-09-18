# New presolve rules — September 17, 2026

Branch `new-presolve-rules`. Three rule families were added on top of `main`
(509e43f) and evaluated one at a time on the 98 Netlib and 138 Maros–Mészáros
problems: dual propagation, dominated columns, and quadratic elimination of
coupled free columns. Two further candidates from the initial survey were
dropped before implementation: forcing rows (2107 rows exist in the raw
problems but the existing pipeline already removes every one) and row
domination (13 rows across 9 problems). The implied-free pivot extension
(25 equality rows) was left out as well; its gain did not justify a change to
the substitution candidate filter.

Before the timing runs, three agents did an efficiency pass over the whole
pipeline with the constraint that reductions stay identical; their patches
are integrated here and their reports are under `docs/benchmarks/new-rules-20260917/`.

## Reductions (default settings, `run size`)

Totals over the 233 problems that stay `Reduced` in both runs (`analyze.py
size main-size final-size`):

| Measure | main | new rules | change |
| --- | ---: | ---: | ---: |
| variables | 770667 | 767377 | -3290 (-0.43%) |
| linear rows | 545118 | 544569 | -549 (-0.10%) |
| A nonzeros | 3556472 | 3528655 | -27817 (-0.78%) |
| P nonzeros | 833833 | 833833 | 0 |
| finite bound sides | 719888 | 716746 | -3142 (-0.44%) |

GENHS28, HS51 and HS52 are now solved outright by presolve. 56 problems
change size; the largest: MAROS-R7 (-860 variables), 80BAU3B (-422 variables,
-308 rows), STANDATA/STANDGUB/QSTANDAT (-337, -337, -325 variables), WOODW
(-242), MODSZK1 (-162, -2 rows), D2Q06C (-91, -67 rows), WOOD1P (-84).

Per rule, from the in-process scanner (kept locally under `benchmark/results/new-rules-20260917/scan/`, which git ignores), cumulative:

| Rule | Variables | Rows | Where it fires |
| --- | ---: | ---: | --- |
| Dual propagation (verified) | -1481 | -474 | 35 problems: MAROS-R7 860 fixed columns; 80BAU3B 377 columns and 300 rows; D2Q06C 71/65, FINNIS 29/25, FIT1P, SCRS8, STANDATA, SCFXM1-3 |
| Dominated columns (identical support) | -1638 | -67 | 30 problems: STANDATA/STANDGUB/QSTANDAT 324 each, WOODW 242, 80BAU3B 45, WOOD1P 39, SHELL 35, AGG2/AGG3 32, BNL2 25, VTP.BASE 23 |
| Quadratic elimination | -6 | 0 | GENHS28, HS51, HS52 (solved); STCQP1/2, MOSARQP1/2, QSEBA have bounded coupled columns and stay |

The general dominated-column search (`dominated_columns.general_search`,
enabled by the aggressive preset) finds a further ~330 variables in DFL001,
D2Q06C and PILOT.JA at a cost of roughly 7% of corpus time, so it is not on by
default.

## Timing (quiet machine, fresh process per trial, 5 trials, `--rule all`)

`analyze.py time base-c ...`. `base-c/d/e` are the `main` executable built in
a separate worktree; the other runs are this tree with rule families disabled
through the new `--without` option. Sums are per-problem medians over all 236
problems; the geometric mean weights small problems equally.

| Run | Rules on top of the optimized pipeline | Sum of medians | vs `main` | Geomean ratio vs `main` |
| --- | --- | ---: | ---: | ---: |
| base-c | `main` | 634.6 ms | — | 1.000 |
| base-d | `main` (repeat) | 632.6 ms | -0.3% | 1.000 |
| base-e | `main` (repeat) | 634.1 ms | -0.1% | 0.988 |
| opt-none | none | 532.8 ms | -16.0% | 0.885 |
| opt-default | quadratic elimination (the new default) | 534.0 ms | -15.9% | 0.892 |
| opt-dom-quad | + dominated columns | 544.6 ms | -14.2% | 0.911 |
| opt-dual | dual propagation only | 564.2 ms | -11.1% | 0.943 |
| opt-dual-dom | dual propagation + dominated columns | 570.4 ms | -10.1% | 0.958 |
| opt-all | all three | 570.7 ms | -10.1% | 0.962 |

Marginal cost of each rule on the optimized pipeline:

| Rule | Sum of medians | Geomean | Decision |
| --- | ---: | ---: | --- |
| Quadratic elimination | +0.2% | +0.8% (within the ±1% baseline noise) | on by default |
| Dominated columns (identical support) | +2.0% | +2.1% | aggressive preset only |
| Dual propagation | +5.9% | +6.5% | aggressive preset only |

Dual propagation was first scheduled once per medium phase, which cost about
17% of presolve time (58% of it in the implied-bound sweep); moving it to a
single pass after redundant-bound removal, where implied-free columns are
already exposed, kept the reductions and cut the cost to the 6% above, most of
it on problems where it proves nothing (the pass is one sweep over the dual
rows). Dominated columns first searched from each column's shortest row (+8 to
+13%); the identical-support test inside the existing parallel-column sort
keeps 82% of the reductions at +2%.

## Dual propagation rebuilt on verified directions

The first version applied its conclusions on the strength of the dual
argument alone, which assumes an optimal solution exists and, for conic
problems, a zero duality gap. It was replaced: every multiplier tightening now
records the dual row that produced it, each conclusion's proof is accumulated
into a primal direction, and the direction is verified exactly on the current
model (variables move away from finite bounds, row activities stay on their
allowed sides with exact zeros on equalities, `cᵀd ≤ 0`). Only then is the
column fixed or the row restricted, by the same shift argument dual fixing and
dominated columns rely on. Exact verification drops directions whose
coefficient ratios do not reproduce zeros in floating point.

| Version | Variables | Rows | Note |
| --- | ---: | ---: | --- |
| assumption-based | -1627 | -423 | MODSZK1 -262/-102, BNL2 -34 |
| verified directions | -1481 | -474 | MODSZK1 -1 (decimal ratios fail the exact test); 80BAU3B, D2Q06C, FINNIS, STANDATA gain rows |

In-process cost of the verified pass on the corpus: +7% sum, +8.5% geometric
mean (MAROS-R7 14.8 → 20.0 ms for 860 extracted directions), so the rule
remains in the aggressive preset only. `Stats::dual_reductions` and the
certificate sign check no longer exist.

## Efficiency pass (reductions identical to the pre-pass tree)

The three agent reports are in `agents/`. Integrated changes: skip fruitless
scans and reuse round buffers in dual propagation; a per-row singleton count so
the singleton-column scan skips 89% of the entries it walked; redundant-bound
removal without dead queue pushes; packed 64-bit fingerprint keys with a radix
sort; a const-generic axis for linked-matrix traversal; smaller list headers
and activity cache entries; merged activity queues; a per-row equation
snapshot cache; an arena-backed symmetric Hessian; scratch reuse in
`replace_row`. Together: -16.0% corpus presolve time with no rule enabled.

Two later changes alter reductions slightly and are noted for the record:
`fix()` now updates cached activities incrementally instead of invalidating
them (GREENBEB +2 and PILOT87 -1 variables from roundoff-level propagation
differences), and quadratic elimination changes the aggressive preset's
substitution order on AUG2DC (-204 variables, +931k Hessian nonzeros) while
improving AUG3DC (-739 variables, -11.7k Hessian nonzeros). The aggressive
preset otherwise gains -3930 variables and -686 rows over `main`
(`analyze.py size main-aggressive-size aggressive3-size`).

## Method

- Release build on macOS/arm64; `main` at 509e43f built in a git worktree.
- Two back-to-back `main` runs agree within 0.3% in the sum of medians and
  about 1% in the geometric mean; per-problem changes below 5% are noise.
- Sizes: `run size` for both presets, compared with `analyze.py` since the
  built-in `compare` refuses runs whose settings differ.
- Candidate-rule detection before implementation: `scan/scan.rs` (an
  in-process detector run on the reduced model of each problem), with its
  detection tables in `scan/`.
