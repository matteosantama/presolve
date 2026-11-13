# Rule scheduling

Presolve alternates inexpensive cleanup with broader searches because one reduction can expose another. For example, substitution can create a singleton row, bound transfer can fix its variable, and fixing can expose an empty row. Enabling a family allows it to participate at the stages below; applicability and budgets can still skip it.

## Initial LP folding

When enabled and time remains, LP folding runs before ordinary substitutions so that those substitutions do not first break the original symmetry. For eligible linear models without cone blocks, enabled singleton-row elimination first transfers row-encoded bounds to variables. Folding itself excludes quadratic and conic models.

## Cleanup

Cleanup repeats the following order until a pass makes no model edits:

1. Fixed-variable elimination.
2. Cone simplification, when applicable.
3. Empty-column elimination, including eligible quadratic elimination.
4. Simple dual fixing.
5. Singleton-row elimination.
6. Empty-row elimination.
7. Cone simplification again, when applicable.

Each entry respects its switch. A variable fixed by singleton-row elimination can therefore be removed on the next cleanup pass. This local repetition is independent of the outer progress threshold, and the overall time budget is soft rather than a hard interruption of every cleanup operation.

## Fast phase

Each outer iteration begins with cleanup. A fast phase then runs:

1. Singleton-column elimination.
2. Cleanup.
3. Implied-free equality elimination.
4. Doubleton equality elimination.
5. Short equality elimination.
6. Cleanup.

Fast phases repeat while the selected progress policy reports sufficient progress for that phase. When progress is insufficient, the next iteration moves to the medium phase rather than stopping immediately.

## Medium phase

After its initial cleanup, a medium phase runs:

1. Bound propagation, including permitted extra rounds.
2. Coupled dual fixing.
3. Cleanup.
4. Implied-free and short equality elimination.
5. Cleanup.
6. Parallel-row elimination.
7. Parallel-column elimination, including the enabled dominated support-group search.
8. General dominated-column search, only when both the rule and `general_search` are enabled.
9. Cleanup.

The end of the medium phase completes an exploration cycle. If the cycle made sufficient progress, exploration restarts with fast phases. Otherwise the main loop ends and final passes begin. Under `Progress::AnyChange`, remaining doubleton and longer equality candidates are also reconsidered on the next cycle because other edits may have changed their eligibility.

### Extra propagation rounds

A propagation stage always starts with one round. It can perform up to `propagation.additional_rounds` extra rounds while tightenings continue, affected-row work fits the shared allowance, and time remains.

Before each extra round, cleanup and enabled singleton-column, implied-free, doubleton, and short equality rules process the consequences of tighter bounds, followed by cleanup again. These extra rounds do not repeat global parallel scans.

## What counts as progress?

The default policy measures constraint plus Hessian nonzeros, rather than only the number of variables or rows. Fast phases compare their own before/after sizes. Full cycles compare against the baseline established before the first cycle or after the previous medium phase.

For size \\(N\\) before and \\(N'\\) after, the default requires

\\[
N>0,\qquad N'<(1-0.05)N.
\\]

The reduction must be strictly greater than 5%. This avoids spending repeated searches on small structural gains and accounts for substitutions that shrink the constraint matrix but densify the Hessian.

`Progress::AnyChange` instead continues after any model edit, including bound-only changes. It is more persistent but still subject to candidate restrictions, numerical checks, work allowances, and time limits. Neither policy promises a mathematical fixed point. See [Settings](./settings.md#time-progress-and-threads) for configuration.

## Final passes

If the time limit has not been reached, the following enabled families run in order:

1. Equality dependency detection.
2. Sparsification.
3. Redundant variable-bound removal.
4. Implied-free equality elimination, revisiting opportunities exposed by removed bounds.
5. Bound shifting.
6. A single dual-propagation pass on the final model.
7. Network equality reductions.
8. Convex dominance.
9. A final parallel-row scan when `final_parallel_scan` and the parallel-row rule are enabled.

These passes do not restart the entire exploration cycle. Successful dependency or sparsification reductions trigger local repetition of cleanup, propagation, and equality substitutions. Successful final implied-free, bound-shift, dual-propagation, network, or parallel-row changes trigger local cleanup and substitutions without another propagation or redundant-bound sweep. Convex dominance triggers cleanup when it changes the model. Local repetition stops when there are no edits or the deadline is observed.

## Stopping and returned results

The main loop stops when cycle progress is insufficient or the time budget is observed to have expired. Rule-specific work caps can stop an individual search earlier. Valid reductions already applied remain usable when a budget is exhausted; that is not a claim that no further reductions exist.

A certificate can terminate exploration before later rules run. The public result distinguishes an unchanged problem, a reduced problem with postsolve mapping, a solved problem, infeasibility, and unboundedness with both a feasible point and improving direction. An intermediate improving direction alone does not establish feasibility.

Use `result.stats.time_limit_reached` to distinguish observed time exhaustion from ordinary completion, and inspect phase counts and reduction statistics when deciding whether broader settings are useful. Initialization of the reusable presolver is excluded from the per-call budget; construction of the working model and export preparation are included.
