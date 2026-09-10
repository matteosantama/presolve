# Aggressive presolve configuration

The recommended starting point is `Settings::aggressive(Duration::from_secs(2))`.

It prioritizes remaining variables and constraints while allowing arbitrary fill and net Hessian growth. The existing default configuration remains unchanged. All candidate and resource limits are independently configurable.

```rust
use presolve::Settings;
use std::time::Duration;

let settings = Settings::aggressive(Duration::from_secs(2));
```

## Configuration

- Unlimited substitution fill; Hessian growth allowed.
- Generalized equality substitution: row length at most 16, pivot column degree 1–16, bounded and quadratic pivots allowed, no net-nonzero restriction, unlimited estimated work allowance.
- Continue fast phases and complete exploration cycles after any model edit, even when nonzeros increase. Revisit equality candidates whose pivot degree or curvature may have changed.
- No extra propagation round or work cap. Existing numerical thresholds remain in effect.
- Equality-only sparsification, using its default work allowance; no auxiliary activity variables.
- All rule families enabled; one execution thread; two-second soft budget per problem.

The 16-entry caps restrict which substitutions are attempted, not how much fill they may create. Larger candidates can consume the time budget before cheaper reductions are exhausted. To admit every row and column length, set both `settings.equalities.max_row_length` and `settings.equalities.max_column_length` to `usize::MAX`.

## Confirmation on 236 problems

98 Netlib LPs and 138 Maros–Mészáros QPs, on macOS/aarch64. Three trials for each problem and configuration (1,416 measurements), alternating configuration order deterministically. Each trial used a fresh process and one thread, with a two-second soft budget. Timings include `Presolver` construction and presolve, excluding file loading, process startup, and destruction. There is no warm-up in this comparison.

Counts below are totals of remaining model sizes. Timings use the median of the three trials for each problem.

| Metric | Default | Aggressive | Change |
| --- | ---: | ---: | ---: |
| Variables | 938,510 | 702,001 | -25.2% |
| Linear constraints | 706,897 | 565,181 | -20.0% |
| Finite variable-bound sides | 727,592 | 563,936 | -22.5% |
| Constraint nonzeros | 4,037,402 | 5,407,762 | +33.9% |
| Hessian nonzeros (full symmetric) | 3,262,060 | 11,545,175 | +253.9% |
| Sum of per-problem median times | 0.748 s | 10.438 s | 13.96× |
| Median problem time | 0.918 ms | 1.504 ms | |
| 95th-percentile problem time | 12.361 ms | 176.416 ms | |

235 problems finished within the budget. CONT-200 reached the soft budget in all three aggressive trials; its partial reductions are included in these counts. The longest measured aggressive call was 2.018 seconds. The limit is soft, and it is not a memory limit.

Model sizes and outcomes were identical across the three trials for every configuration/problem pair. The aggressive preset improved at least one of variable or linear-constraint count without worsening the other on 157 problems; 62 kept the same dimensions; 17 had a larger count in at least one dimension. Five previously unchanged QPs became reduced. Neither configuration reported infeasibility or unboundedness. These are presolve measurements, not downstream solver-time measurements.

## Candidate-limit exploration

Single trials per problem, with the same two-second budget. Both row and column length caps use the listed value; other aggressive settings are the same. These runs guided selection; the table above is the repeated confirmation.

| Candidate limit | Remaining variables | Remaining linear constraints | Total presolve time | Problems reaching budget |
| --- | ---: | ---: | ---: | ---: |
| 8 | 749,302 | 589,179 | 6.947 s | 0 |
| 12 | 719,508 | 572,517 | 8.858 s | 1 |
| 16 | 702,001 | 565,181 | 10.693 s | 1 |
| 32 | 713,157 | 593,063 | 15.940 s | 2 |
| 64 | 762,323 | 628,135 | 24.557 s | 7 |
| 128 | 808,893 | 648,339 | 28.121 s | 8 |
| unlimited | 854,412 | 671,358 | 49.894 s | 21 |

Lifting only the substitution fill cap and Hessian restriction, while retaining default searches, reduced aggregate variables by 0.4% and linear constraints by 0.2%. Broader candidate eligibility and continued exploration provide most of the additional reduction.

At limit 16, disabling sparsification left 795 more variables and 228 more linear constraints. Allowing inequality references removed 200 additional variables overall, but left 31 more linear constraints and 932 more finite bound sides; it also increased runtime. Equality-only sparsification is the selected compromise.

## Correctness and compatibility checks

- All 236 default outcomes and model sizes matched the pre-change release binary.
- Tests cover Hessian growth versus fill limits, row/column eligibility, retained bounds with either pivot sign, objective equivalence, primal feasibility and KKT stationarity, postsolve recovery, propagation budgets, progress on bound-only changes, candidate revisiting, auxiliary-variable control, and transactional deadline rejection.
- The initial unrestricted exploration exposed invalid infeasibility certificates caused by cancellation roundoff becoming a later pivot. A substitution now rejects nonzero updated constraint coefficients below `1e-10` of the larger contributing term; it preserves exact cancellations. A regression test covers near-dependent equalities. After this fix, the complete exploration and confirmation runs reported no infeasible or unbounded outcomes. The initial `screen.jsonl` is diagnostic data from before the fix and must not be used as a performance recommendation.
- Outcome consistency is not an independent proof of every reduced benchmark model; targeted tests check the transformations and recovery algebra.

## Reproducing measurements

```sh
cargo build --release -p benchmark
target/release/benchmark run time --name default-2s --rule all --trials 3 --preset default --time-limit-ms 2000
target/release/benchmark run time --name aggressive-2s --rule all --trials 3 --preset aggressive --time-limit-ms 2000
```

These commands save all results even when a case reaches its budget; the benchmark exits nonzero to flag those cases. The comparison command rejects different algorithm settings, so inspect the saved measurements when comparing presets. To reproduce the interleaved ordering used here, the local experiment runner is `results/aggressive-20260909/confirm.py` (choose a fresh output filename before rerunning).

Additional benchmark overrides: `--preset fill`, `--preset unrestricted`, `--equality-row-limit N`, `--equality-column-limit N`, and `--sparsification all|equalities|off`. Full settings are recorded in the benchmark run metadata.

Local raw measurements and analysis are in `benchmark/results/aggressive-20260909/`: `refine.jsonl`, `final-screen.jsonl`, `confirm.jsonl`, and their summary JSON files. The results directory is ignored by Git.
