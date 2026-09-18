# Source reorganization — September 18, 2026

Branch `new-presolve-rules`, on top of d2619ed. The private `core` module was
split into `model` (the working model, its caches and queues, and the recovery
tape, formerly `postsolve/tape.rs`) and `rules` (the rule families, with the
scheduler as the module root). The rayon wrapper moved to `src/executor.rs`.
Every file is a `git mv`; the public modules `matrix`, `problem`, `postsolve`,
`result`, and `settings` are unchanged, so the public API is identical.

## Reductions (`run size`, both presets)

`compare size base-default new-default`: 236 matched cases, 0 changed.
`compare size base-aggr new-aggr`: 235 matched cases, 0 changed (CONT-200
reaches the time limit on both binaries, as before).

## Timing (fresh process per trial, 5 trials, `--rule all`, default preset)

The baseline binary is d2619ed built in a separate worktree. Runs alternated
between the two binaries. `new-1` and `base-2` ran back to back and were both
disturbed (sums 40 ms high, single problems spiking 170–225%), so they are
excluded; the six remaining runs, in the order they ran:

| Run | Sum of medians | Geomean vs base-1 |
| --- | ---: | ---: |
| base-1 | 530.4 ms | 1.000 |
| new-2 | 534.4 ms | 1.007 |
| base-3 | 537.3 ms | 1.000 |
| new-3 | 536.3 ms | 1.006 |
| new-4 | 536.1 ms | 1.002 |
| base-4 | 531.6 ms | 0.999 |

Mean of the baseline sums 533.1 ms, of the new sums 535.6 ms (+0.5%), with
the ranges overlapping. The ten largest problems move within ±5% in both
directions. The `__TEXT` section of the two benchmark binaries has the same
size (1,949,696 bytes), as expected for a change that only renames module
paths.

Raw runs stay in the untracked `benchmark/results/reorg-20260918/`.
