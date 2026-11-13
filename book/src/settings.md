# Settings

Configure a `Settings` value before constructing a `Presolver`. `Presolver::default()` uses the default configuration; `Presolver::new(settings)` returns a result because worker-pool creation can fail. Settings are reused for each presolve call, while the time budget starts afresh for each call.

This example is for a Cargo project that already depends on this library; it cannot run in the standalone Rust Playground.

```rust,ignore
use presolve::{Presolver, Settings};
use std::time::Duration;

let mut settings = Settings::default();
settings.time_limit = Duration::from_secs(5);
settings.rules.dominated_columns = true;
settings.dominated_columns.general_search = true;
let presolver = Presolver::new(settings).expect("initialize presolver");
// let result = presolver.presolve(problem);
```

## Rule selection

`settings.rules` provides independent switches. `Rules::none()` disables every family and can be used as a starting point for selecting a subset. Disabling a consequence rule can leave a constraint or variable for the solver even when another reduction exposes it.

| Enabled by default | Disabled by default |
| --- | --- |
| Fixed variables, empty rows and columns, quadratic elimination | Dual propagation |
| Singleton rows and columns, doubleton and short equalities, implied-free equalities | Bound shifting |
| Bound propagation, redundant-bound removal | Equality dependency detection |
| Parallel rows and columns | Dominated columns |
| LP folding, convex dominance, network equalities, sparsification, cone rules | |

A switch permits a rule to run; eligibility, numerical, work, and time checks can still prevent a reduction. Some families apply only to linear, quadratic, or conic models. See [Rule scheduling](./scheduling.md) for dependencies and ordering.

The documented rules use these fields within `settings.rules`:

| Rule | Switch |
| --- | --- |
| Fixed variables | `fixed_variables` |
| Singleton rows / columns | `singleton_rows` / `singleton_columns` |
| Empty rows / columns | `empty_rows` / `empty_columns` |
| Bound propagation | `bound_propagation` |
| Doubleton / short equalities | `doubleton_equalities` / `short_equalities` |
| Parallel rows / columns | `parallel_rows` / `parallel_columns` |
| Dominated columns | `dominated_columns` |
| Dual fixing / propagation | `dual_fixing` / `dual_propagation` |
| Implied-free equalities | `implied_free_equalities` |
| Redundant-bound removal | `redundant_bounds` |

## Time, progress, and threads

| Setting | Default | Meaning |
| --- | --- | --- |
| `time_limit` | 60 seconds | Soft per-call budget including model construction and export preparation, excluding presolver initialization. |
| `progress` | `Progress::Nonzeros { minimum_reduction: 0.05 }` | Repeat phases and cycles only after a strictly greater than 5% reduction in constraint plus Hessian nonzeros. |
| `threads` | `1` | Serial execution. `0` lets Rayon choose; a positive value caps the pool size. |
| `final_parallel_scan` | `true` | Revisit parallel rows after final reductions when the parallel-row rule is enabled. |

`Progress::AnyChange` continues after any model edit, including bound changes and substitutions that increase nonzeros. It does not guarantee a mathematical fixed point. Nonzero thresholds are clamped to `[0, 1]`; NaN uses `0.05`.

The time limit is soft: an in-progress transaction can finish before the limit is observed. A time limit ends exploration, not validity of reductions already accepted. Inspect `result.stats.time_limit_reached`, elapsed time, phase counts, and before/after sizes to assess a run.

Threads accelerate candidate fingerprinting and sorting when scans are large enough. Reductions and recovery records are applied in deterministic serial order; more threads do not change the mathematical rules.

## Numerical tolerances

| Field in `numerics` | Default | Purpose |
| --- | --- | --- |
| `feasibility` | `1e-9` | Relative margin for applicable contradiction checks. |
| `parallel` | `1e-12` | Relative coefficient comparison for potential proportional rows or columns. |
| `huge_bound` | `1e7` | Reject general propagation candidates with absolute magnitude at least this value. |

The common feasibility test separates a proposed lower value \\(L\\) from an upper value \\(U\\) when

\\[
L>U+\varepsilon\bigl(1+\max(|L|,|U|)\bigr).
\\]

This is not a general rounding instruction. Fixed-variable detection still requires exactly equal finite bounds, geometric rules require exact structural zeros, and parallel-column aggregation requires exact scaled coefficient agreement. Parallel-row merging has additional guards against tolerance-only contradictions and newly created equalities.

## Equality substitution

`substitution_fill`, default `64`, caps newly introduced constraint and quadratic coefficients per shared substitution. `usize::MAX` removes that cap. `allow_hessian_growth`, default `false`, separately controls growth in total Hessian nonzeros.

`allow_unit_doubleton_fill`, default `true`, permits linear unit-slope doubleton replacements beyond the ordinary gross fill cap. Disable it to enforce that cap for every substitution.

The `equalities` settings control short-equality candidates and are also used by implied-free equality elimination:

| Field | Default | Effect |
| --- | --- | --- |
| `max_row_length` | `8` | Maximum short-equality length; this rule starts at three entries. |
| `min_column_length`, `max_column_length` | `2`, `8` | Eligible candidate column lengths. |
| `require_free_variable` | `true` | Require no finite variable bounds for short-equality candidates; ignored by implied-free elimination. |
| `require_linear_variable` | `true` | Exclude quadratic participation, including cross terms. |
| `relative_pivot` | `1.0` | Minimum pivot magnitude divided by the row's maximum coefficient magnitude. |
| `max_pivot_attempts` | `1` | Maximum retained candidates per visit to an equality; zero disables attempts. |
| `cost_aware` | `false` | Rank by column length, or by estimated constraint and quadratic work when enabled. Short-equality elimination first prefers free variables; implied-free elimination has no such preference. |
| `preserve_nonzeros` | `true` | Further cap fill by the number of coefficients removed. |
| `work_limit` | `WorkLimit::Default` | Twice current constraint nonzeros for a short-equality pass; sixteen times for the implied-free aggregation pass. |

`relative_pivot` outside `(0, 1]` falls back to `1.0`. The default avoids slopes larger than one in magnitude; smaller ratios permit amplification. Short-equality substitution retains necessary variable bounds as constraints. Implied-free elimination instead requires all explicit bounds to follow from the defining equality and the other variables' bounds, independently of `require_free_variable`. Its row and column length caps apply only to candidates with more than two column entries. Affected equalities can be revisited during the same implied-free pass.

More candidate attempts or a larger work allowance do not bypass fill or numerical checks. Short-equality statistics distinguish structural, pivot, work, numerical, constraint-fill, quadratic-fill, Hessian-growth, and deadline rejections; these detailed counters do not cover implied-free elimination.

## Bound propagation

| Field in `propagation` | Default | Effect |
| --- | --- | --- |
| `minimum_relative_gain` | `0.01` | Relative gain threshold against the old finite bound magnitude. |
| `minimum_gain_factor` | `1e4` | Absolute gain floor multiplied by `numerics.feasibility`. |
| `additional_rounds` | `3` | Extra rounds after the initial propagation pass. |
| `work_limit` | `Default` | Shared allowance of `max(A nonzeros / 4, 256)` for extra rounds. |

A finite, non-fixing update must exceed the larger of the relative threshold and absolute floor. An update exactly meeting the opposite bound bypasses this gain filter, but not the earlier finite-value, huge-bound, and feasibility checks. Invalid negative or nonfinite gain settings use their defaults. Singleton-row transfers bypass the general propagation gain and huge-bound filters.

## Dominated columns

`rules.dominated_columns` defaults to `false`. Enabling it permits the support-group search within an enabled parallel-column scan. `dominated_columns.general_search`, also default `false`, additionally enables the shortest-row search; both flags must be true for that general pass to run. Its `work_limit` defaults to twice the current constraint nonzeros per pass.

The search considers only variables without quadratic participation or incident cone rows. The general search uses rows of at most 256 entries; support groups compare at most the next 32 neighbors per member. These search widths are not adjustable settings. A bounded search can miss a valid dominance relation.

## Other rule budgets

`WorkLimit::Default` selects a family-specific allowance, `Entries(n)` supplies an explicit estimated-entry budget, and `Unlimited` removes the work cap while leaving time limits in effect. These count estimated sparse work, not milliseconds.

| Family | Default controls |
| --- | --- |
| `dependencies` | Scratch row length 128, basis size 64, work allowance four times constraint nonzeros. |
| `dual_propagation` | Work allowance four times constraint nonzeros. |
| `bound_shift` | Fill 64, column length 64, ratio limit `1e4`, work allowance 64 times constraint nonzeros. |
| `folding` | At most 16 refinement rounds, work allowance 64 times constraint nonzeros plus linear scans. |
| `sparsification` | Inequality references may introduce auxiliary variables; work allowance `max(8 * (A nonzeros + bound entries), 1024)`. |
| `convex_dominance` | Require a strictly greater than 5% structural saving before applying a batch. |

## Aggressive preset

`Settings::aggressive(time_limit)` favors dimensional reduction within the supplied soft budget. It allows unlimited substitution fill and Hessian growth, permits bounded and quadratic equality candidates, raises equality length caps to 16, and uses `Progress::AnyChange`.

It also enables equality dependency detection, dual propagation, and dominated columns including general search. Propagation uses a relative gain of `0.005`, unlimited extra rounds and work, and convex dominance accepts any nonempty saving. Sparsification uses equality references only, introducing no auxiliary activity variables. Execution still defaults to one thread. Every field remains configurable; this preset can spend more time or produce denser reduced objectives.
