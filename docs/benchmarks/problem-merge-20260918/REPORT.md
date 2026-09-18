# Merging `ProblemData` into `Problem` — September 18, 2026

Branch `new-presolve-rules`, on top of 6cd59c1. `ProblemData` and the
accessor layer on `Problem` were replaced by a single `Problem` struct with
public fields. Only the two matrix fields are opaque: `a: ConstraintMatrix`
and `p: Option<QuadraticMatrix>`, each built from a `CscMatrix` with `into()`
and exported with `as_csc()` (a free borrow when the storage is CSC) or
`into_csc()` (one pack). `QuadraticRef` is gone; `Problem::into_csc()` now
packs both matrices in place. The lazy export and the zero-copy unchanged path
are unchanged, so no work moved.

## Reductions (`run size`, both presets)

`compare size base-default new-default`: 236 matched cases, 0 changed.
`compare size base-aggr new-aggr`: 235 matched cases, 0 changed (CONT-200
reaches the time limit on both binaries, as before).

## Timing (fresh process per trial, 5 trials, `--rule all`, default preset)

The baseline binary is 6cd59c1 built in a separate worktree. Four pairs ran
back to back, alternating which binary went first:

| Run | Sum of medians | Geomean vs base-1 |
| --- | ---: | ---: |
| base-1 | 525.4 ms | 1.000 |
| new-1 | 521.3 ms | 0.991 |
| new-2 | 522.9 ms | 0.995 |
| base-2 | 519.5 ms | 0.984 |
| base-3 | 521.8 ms | 0.985 |
| new-3 | 520.3 ms | 0.981 |
| new-4 | 504.4 ms | 0.956 |
| base-4 | 520.1 ms | 0.982 |

Mean of the baseline sums 521.7 ms; of the new sums 517.2 ms, or 521.5 ms
without the unusually quiet `new-4`. The ten largest problems move within ±5%
in both directions. Raw runs stay in the untracked
`benchmark/results/problem-merge-20260918/`.
