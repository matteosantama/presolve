# Tier 1 refactors — September 18, 2026

Four refactors from `docs/review-20260918/` were applied one commit each on
top of `d65c4e2`, with `run size` on both presets after every commit
(236 default and 235 aggressive cases matched, 0 changed, every time):

1. `7d76acb` one implied-bound kernel (`Extreme::implied`,
   `Activity::implied`, `Activity::implied_side`); `bound_implied` with a
   recompute flag replaces `open_sides`; `linear_column` on the model.
2. `926e87f` model bookkeeping helpers: `column_changed`, `shift_cached`,
   `clear_row`, `point`, `primal_certificate`, `dual_certificate`.
3. `998207d` one `settings` field on the model, `schedule::Limits` removed,
   fill parameters dropped from the rule entry points.
4. `a9e3d01` entry point split into `working_model`, `finish`, `unchanged`;
   `pack` into `survivors`, `build_postsolve`, `compact_problem`;
   `OriginalMap` and `Coordinates` own index translation; the problem's
   `Matrix` enum keeps only `Csc` and `Linked`.

## Timing

Fresh-process runs, 5 trials each, all on a quiet machine. The first
measurement after the four commits read +5% (min over runs) and the
bisection put all of it in commit 1: the kernel's small generic helpers live
in another module and were not inlined into the propagation hot loop across
codegen units. `#[inline]` on the six activity helpers (`d832184`, which also
moves the survivor vectors into the postsolve map instead of cloning them)
removed it. In-process repeated timing of the corpus was flat throughout
(-0.8% sum, geomean 0.999), which is what isolated the effect as cold-start.

Final matched comparison, pre-refactor binary versus `d832184`:

| Pairing | Sum of medians | Geomean |
| --- | ---: | ---: |
| three runs each, min over runs | 516.4 → 519.6 ms (+0.6%) | 1.023 |
| three runs each, median over runs | 527.9 → 528.3 ms (+0.1%) | 1.017 |
| same session, base9 versus e-b/e-c, min | 528.7 → 520.3 ms (-1.6%) | 0.995 |

The ten largest problems are within run-to-run noise in every pairing
(BOYD1 75.9 → 75.3 ms, BOYD2 45.0 → 45.3, CONT-300 34.6 → 34.7). The
geometric mean is driven by problems under 100 µs, where fresh-process
jitter is several percent and the sign flips between sessions.

`t1-base8.jsonl` was a partial run (killed and not overwritable) and was
deleted; `t1-base6/7/9` are the complete baseline runs of the same binary.
