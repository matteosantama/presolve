# Efficiency review, tier 1 — September 18, 2026

Branch `new-presolve-rules`. The tier-1 items of the efficiency review
(`docs/review-efficiency-20260918/REPORT.md`) were landed in eight commits on
top of d8471f1, one group at a time. Each group was gated on identical
reductions on both presets (`run size`, 236 and 235 matched cases, 0 changed
every time) and on two alternating fresh-process timing pairs against the
previous commit (5 trials, `--rule all`, default preset, quiet machine).

| Commit | Group | Pair means (base → new) |
| --- | --- | ---: |
| 2739e2c | `#[inline]` on the hot leaf helpers called across codegen units | 547.5 → 539.9 ms (-1.4%) |
| c023232 | `implied_bound` tests for a tighter bound before the gain threshold | 537.2 → 534.9 ms (-0.4%; CONT-300/201 -3%) |
| 9d9274d | short-equality candidates rejected before their column loads | 522.0 / 543.4 → 515.8 ms |
| cb46cc5 | temporaries removed from `fix`, `remove_variable`, sparsification, singleton scans, parallel groups | 512.6 → 507.9, 528.6 → 524.2 ms (-0.9%) |
| 7e41e43 | model seeded from the arena in storage order; arena builder lookahead | 529.9 → 524.3 ms (-1.0%) |
| f218496 | MSD radix candidate sort with a skewed-key fallback | 525.5 → 520.3 ms; CONT-200/201 -5 to -7% |
| ad91cea | lazy Hessian fingerprints; exact repeat scans skipped | flat sum, geomean 0.975; FIT2P -9%, EXDATA -19% |
| e2a570b | Hessian row cursor in the pair scan; split, branch-free activity kernel; O(1) sizes; coupled-fix repeat skip | 505.6 → 495.9 ms (-1.9%) |

## End to end (d8471f1 → e2a570b)

| Run | Sum of medians | Geomean vs base-1 | Faster > 5% | Slower > 5% |
| --- | ---: | ---: | ---: | ---: |
| base-1 | 546.0 ms | 1.000 | | |
| new-1 | 493.9 ms | 0.893 | 185 | 9 |
| new-2 | 494.9 ms | 0.895 | 175 | 3 |
| base-2 | 547.7 ms | 0.993 | 35 | 32 |

Sum of medians -9.5%, geometric mean 0.89; the "slower" entries are
sub-100 µs problems inside run-to-run noise. Largest problems (median ms,
base-1 / new-1): BOYD1 79.8 / 76.1, BOYD2 46.5 / 42.7, CONT-300 37.2 / 31.7,
CONT-200 17.2 / 14.5, CONT-201 15.9 / 13.6, PILOT87 15.4 / 14.2, FIT2P
12.5 / 11.2, EXDATA 9.6 / 7.7.

## Notes

- One verification pair was invalid and repeated: group H's edits had been
  drafted in the working tree before group G was committed with `git add
  src`, so the first "H" pair compared identical binaries. The commit was
  split (G alone as ad91cea) and H re-verified against it.
- Two single-run spikes (BOYD1 +20% once under group F, CONT-300 +12% once
  under group G) were re-timed on the same binaries and found flat.
- Not landed from tier 1: the `from_columns` node fill (needs `unsafe` or
  `bytemuck`), and the SIMD short-row split (net zero per the review).
  Tier 2 (tape storage, the 64-byte row record, exit-path passes) is open.

Raw runs stay in the untracked `benchmark/results/t1eff-*/`.
