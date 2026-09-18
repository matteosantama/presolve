# Plain CSC problem fields — September 18, 2026

Branch `new-presolve-rules`, on top of 7940006. `Problem` now holds plain
`CscMatrix` fields: `a`, and `p` as the upper triangle of the symmetric
Hessian on input and output. The `ConstraintMatrix` and `QuadraticMatrix`
wrappers, their private storage enums, the lazy export (`as_csc`, `into_csc`)
and `from_columns` are gone. Presolve copies both matrices into its working
structures, and a reduced outcome packs the survivors straight back into CSC
inside `presolve()`. An unchanged outcome still returns the caller's buffers
untouched, and the input Hessian is reused when it was neither edited nor lost
a column. The constraint pack derives its column pointers from the arena's
column lengths instead of a counting pass over the rows.

## Reductions (`run size`, both presets)

`compare size base-default new-default`: 236 matched cases, 0 changed.
`compare size base-aggr new-aggr`: 235 matched cases, 0 changed (CONT-200
reaches the time limit on both binaries, as before).

## Timing (fresh process per trial, 5 trials, `--rule all`, default preset)

The pack used to happen outside the timed region, in the caller's
`into_csc()`. For a fair comparison the baseline binary (7940006) was built
with its benchmark harness patched to call `into_csc()` on reduced outcomes
inside the timed region, which is what a CSC caller pays. Four alternating
pairs:

| Run | Sum of medians | Geomean vs base-1 |
| --- | ---: | ---: |
| base-1 | 561.4 ms | 1.000 |
| new-1 | 538.9 ms | 0.953 |
| new-2 | 543.3 ms | 0.978 |
| base-2 | 551.0 ms | 1.011 |
| base-3 | 550.1 ms | 1.005 |
| new-3 | 541.9 ms | 0.976 |
| new-4 | 543.8 ms | 0.976 |
| base-4 | 553.5 ms | 1.024 |

Mean of the baseline sums 554.0 ms; of the new sums 542.0 ms (-2.2%), with
the geometric mean about 0.975. The largest problems are within noise. Note
that `run time` now includes the pack for every reduced problem, so its sums
are not directly comparable with earlier reports (about 521 ms for the same
tree without the export). Raw runs stay in the untracked
`benchmark/results/plain-csc-20260918/`.
