# presolve

Sparse presolve and postsolve for convex quadratic and conic optimization in Rust.
The library simplifies a problem before it reaches a solver, then maps solutions
and certificates back to the original variables and constraints.

This README catalogs the implemented rules, their applicability, and how they
interact. Fourteen of the fifteen rule families are enabled by default; the bounded
equality-dependency pass is enabled by the aggressive preset.

## Problem model and terminology

The native problem is

```text
minimize    ½ xᵀ P x + cᵀ x + c₀
subject to  lᵣ ≤ A x ≤ uᵣ
            lₓ ≤ x   ≤ uₓ
            G x + s = h,  s ∈ K
```

`P` is positive semidefinite; omitting it gives a linear objective. Linear and
conic rows share the `ProblemData.a` sparse matrix and are distinguished by
`Constraint::Linear` and `Constraint::Cone`. The notation `A` and `G` above separates
these two row domains for readability.

- A **singleton** row or column has one nonzero; a **doubleton** has two.
- A **free** variable has neither a finite lower bound nor a finite upper bound.
- A row's **activity interval** is the minimum and maximum of its linear
  expression over the current variable bounds.
- A direction is **locked** when moving a variable in that direction can violate
  a constraint. An incident conic row conservatively locks both directions.
- **Fill** means newly introduced nonzero coefficients during a transformation.

The caller validates dimensions, sparse structure, finite coefficients, bound
consistency, convexity, and cone geometry. Supply only the upper triangle of `P`.
Empty `variable_bounds` means all variables are free. See
[problem data](src/problem/mod.rs) and [cone conventions](src/problem/cone.rs).

## Rule catalog

The switch names below are fields of [`settings::Rules`](src/settings.rs).
An enabled rule still checks its structural, numerical, and work limits before
applying. Several rules can perform the same underlying operation: for example,
turning off `fixed_variables` does not prevent `empty_columns` or `dual_fixing`
from fixing a variable themselves.

| Category | Rule switch | Main effect |
| --- | --- | --- |
| Local cleanup | `fixed_variables` | Substitute a variable whose bounds agree. |
| Local cleanup | `empty_columns` | Optimize and remove an independent variable. |
| Local cleanup | `empty_rows` | Delete a satisfied constant linear row or detect infeasibility. |
| Local cleanup | `singleton_rows` | Convert a one-variable linear row to variable bounds. |
| Bounds and optimality | `bound_propagation` | Derive bounds, remove redundant row sides, and detect contradictions. |
| Bounds and optimality | `redundant_bounds` | Remove variable-bound sides implied by retained rows. |
| Bounds and optimality | `dual_fixing` | Use objective derivatives and unlocked directions to eliminate variables. |
| Substitution | `singleton_columns` | Eliminate a variable occurring in one linear row. |
| Substitution | `doubleton_equalities` | Eliminate one variable from a two-variable equality. |
| Substitution | `short_equalities` | Substitute from equalities under configurable candidate and fill limits. |
| Redundancy | `equality_dependencies` | Remove exact linear combinations of equalities using a bounded scratch basis. |
| Parallel structure | `parallel_rows` | Merge proportional linear constraints. |
| Parallel structure | `parallel_columns` | Aggregate interchangeable variables or exploit objective dominance. |
| Matrix sparsity | `sparsification` | Cancel shared coefficients between linear rows. |
| Cone geometry | `cones` | Simplify constant blocks and recognized structural cone patterns. |

### Local cleanup

**Fixed variables — `fixed_variables`.** When a variable has exactly equal,
finite bounds, substitute that value into every incident linear or conic row and
into the objective, then remove the column. Quadratic cross terms update the
remaining linear objective, and constant terms update `c₀`. For example, fixing
`x = 2` changes `3x + y ≤ 10` to `y ≤ 4`. Fixing introduces no new coefficients;
nonfinite transformed values cause the operation to be skipped.
Source: [variables.rs](src/core/rules/variables.rs),
[shared transformations](src/core/model.rs).

**Empty columns — `empty_columns`.** A variable absent from all linear and conic
rows can be optimized separately if it has no off-diagonal Hessian couplings.
For its scalar objective `½ p x² + c x`, the rule chooses `-c/p` clipped to the
bounds when `p > 0`. With `p = 0`, it chooses the lower bound for `c > 0`, the
upper bound for `c < 0`, or zero clipped to the bounds for `c = 0`. A finite
choice is substituted and removed. A linear objective with an infinite improving
side produces a recession ray. Variables coupled through `P` are retained.
Source: [variables.rs](src/core/rules/variables.rs).

**Empty rows — `empty_rows`.** A linear row with no coefficients reduces to
`l ≤ 0 ≤ u`. The rule deletes it if satisfied within the feasibility tolerance,
or produces an infeasibility certificate if a side excludes zero by more than
that tolerance. Conic constant rows are handled as blocks by `cones`.
Source: [rows.rs](src/core/rules/rows.rs).

**Singleton rows — `singleton_rows`.** Convert `l ≤ a x ≤ u` into bounds on `x`,
reversing the sides when `a < 0`, and intersect them with the existing bounds.
For example, `2 ≤ -2x ≤ 6` gives `-3 ≤ x ≤ -1`. Delete the row only after the
variable bounds fully represent its restriction. Contradictions beyond tolerance
produce a certificate; overflow or a small inconsistent interval can leave the
row in place. An equality can expose a fixed variable for the next cleanup pass.
Source: [rows.rs](src/core/rules/rows.rs).

### Bounds and optimality

**Bound propagation — `bound_propagation`.** For each changed linear row with
at least two nonzeros, compute its activity interval from the variable bounds.
The rule performs three related reductions:

1. Detect infeasibility if the attainable activity lies beyond a row side by more
   than the feasibility margin.
2. Remove a row side already guaranteed by the activity interval. Delete the
   entire row when both sides are redundant. Redundancy requires actual
   containment, without a tolerance-based relaxation.
3. Isolate each variable using the other terms' activity interval. For example,
   `x + y ≤ 5` and `y ≥ 2` imply `x ≤ 3`.

For a positive coefficient `aⱼ` and residual activity `[r_min, r_max]`, the implied
bounds are `(l - r_max)/aⱼ ≤ xⱼ ≤ (u - r_min)/aⱼ`; negative coefficients reverse
the bound directions. Unbounded residuals may prevent an implication. The rule
skips very large or insignificant changes as described under numerical controls.
Source: [bounds.rs](src/core/rules/bounds.rs).

**Redundant variable bounds — `redundant_bounds`.** After the main phases, remove
a finite bound side if a retained linear row and the other variables' current
bounds imply a restriction at least as strong. For example, `x + y ≤ 5` and
`y ≥ 2` make an explicit `x ≤ 4` unnecessary. Implications are recomputed as
bounds disappear to avoid circular proofs. This final pass lets explicit bounds
help earlier rules, then removes unnecessary bound constraints before solving.
It is skipped when the main run reports a time limit.
Source: [bounds.rs](src/core/rules/bounds.rs).

**Dual fixing — `dual_fixing`.** Use objective derivatives to prove that an
optimum can be chosen at a variable bound when constraints do not block movement
toward it. This switch controls two passes:

- **Simple dual fixing:** for an uncoupled objective term, fix to a finite lower
  bound if decreasing the variable is unlocked and the derivative there is
  nonnegative. Apply the symmetric test at an upper bound. A nonzero linear cost
  with an infinite improving side gives a recession ray. For a flat objective,
  choose a finite unlocked bound; if that side is infinite, remove the variable
  and its incident rows, saving them so postsolve can choose a finite feasible
  value later.
- **Coupled quadratic dual fixing:** bound the derivative using the other
  variables' bounds. A resolved positive derivative at an unlocked lower bound,
  or negative derivative at an unlocked upper bound, permits fixing. This pass
  requires a positive diagonal and a Hessian row with 2–33 nonzeros. It rejects
  nonfinite derivative bounds and signs unresolved by its roundoff margin.

For example, minimizing `x` with `x ≥ 0` and only upper-sided constraints having
positive coefficients of `x` permits fixing `x = 0`.
Source: [dual_fixing.rs](src/core/rules/dual_fixing.rs).

### Substitution

All substitution rules solve an equality for a pivot variable, replace that
variable in other rows and the objective, and record how to recover it. For
`x = Tz + d`, the quadratic update is

```text
P′  = Tᵀ P T
c′  = Tᵀ(c + P d)
c₀′ = c₀ + cᵀd + ½ dᵀP d
```

Bounds on an eliminated variable that are not proved implied become a linear
row on the remaining variables. Thus a substitution may remove a column while
retaining a transformed row. All these rules respect `substitution_fill`, reject
nonfinite arithmetic, and forbid a net increase in Hessian nonzeros unless
`allow_hessian_growth` is enabled. `substitution_fill = usize::MAX` removes the
allocation cap; it does not by itself allow net Hessian growth.
Sources: [substitution.rs](src/core/rules/substitution.rs),
[model.rs](src/core/model.rs), [objective.rs](src/core/objective.rs).

**Singleton columns — `singleton_columns`.** Consider a variable appearing in
exactly one row of the shared constraint matrix, where that row is linear and
has at least two nonzeros.

- In an equality, substitute the variable after discarding bounds implied by
  the equality and the other variables' bounds. If both bound sides remain,
  defer to doubleton or generalized equality substitution. A variable appearing in `P` is also deferred if
  its coefficient is smaller in magnitude than another coefficient in the row.
- In an inequality or ranged row, the variable must be absent from `P`. Its
  linear cost selects a preferred finite row side: lower when `cⱼ/aⱼ > 0`,
  upper when `cⱼ/aⱼ < 0`. A zero cost selects an available finite side.
  Substitute on that side only when the implied-bound checks establish that it
  is reachable without a blocking variable bound.

For example, a free singleton `x` in `2x + y = 6` can be replaced by `3 - y/2`.
A curved objective generally cannot justify making an inequality tight, because
its optimum may be in the interior.

**Doubleton equalities — `doubleton_equalities`.** An equality
`a x + b y = r` allows replacing one variable by an affine function of the other.
Any retained bounds become a singleton row and may then become bounds on the
survivor. With quadratic involvement, prefer the larger-magnitude pivot to avoid
amplifying curvature. Otherwise prefer a singleton column, an integral
substitution ratio, then the shorter column. The chosen ratio magnitude must
lie in `[1e-7, 1e7]`. For a purely linear pair, the other pivot is tried if the
preferred substitution fails its fill or arithmetic checks.

**Short equalities — `short_equalities`.** By default, consider rows with 3–8
nonzeros and a free pivot absent from `P`, appearing in 2–8 constraint rows.
The pivot must have maximum coefficient magnitude in the equality. Prefer free
variables, then shorter columns. Default substitutions remove both the variable
and equality, leave `P` unchanged, and cap fill by the removed entries.

`settings.equalities` can widen the row and column limits, admit bounded and
quadratic pivots, remove the no-growth restriction on total constraint and
Hessian nonzeros, and change the work allowance. Necessary bounds on a bounded
pivot become a ranged row; bounds proved implied are omitted. `relative_pivot`
defaults to `1.0` (maximum magnitude); lower values permit amplification and
should be evaluated for the intended models. `max_pivot_attempts` defaults to
`1`; raising it tries other eligible pivots after transactional rejection.
`cost_aware` ranks candidates using estimated constraint and Hessian work.
It defaults to `false`: measured speed and Hessian-size improvements came with
mixed dimensional results. Work budgets include quadratic update visits under
either ordering. Finite-arithmetic and cancellation checks remain in effect.
`result.stats.equalities` reports candidate decisions and rejection reasons.
The historical rule name is retained even for long equalities.

**Equality dependencies — `equality_dependencies`.** After the ordinary phases,
construct a bounded scratch basis of equality rows. Delete a row only after
exact arithmetic establishes that both its coefficients and right-hand side
are a linear combination of surviving equalities. No variables or Hessian
entries change in this pass. An inconsistent exact relation can produce a
Farkas certificate. Postsolve transfers deleted-row multipliers for warm starts.

Products and differences must pass exactness checks; very small products whose
rounding residual could underflow are rejected. This is a conservative subset
of numerical rank detection and can miss dependent decimal-valued rows.
`dependencies.max_row_length` (128), `max_basis_rows` (64), and `work_limit`
(default four times A nonzeros) bound scratch work and memory. The pass runs
once and drains consequences through the existing cleanup rules. It is enabled
by `Settings::aggressive`; enable it explicitly elsewhere only when the extra
reductions justify its measured cost.

### Parallel structure

**Parallel rows — `parallel_rows`.** Identify linear rows with matching support
and proportional coefficients, scale their bounds into the same coordinates,
intersect the intervals, and remove the duplicate row. For example,
`x + y ≤ 4` and `2x + 2y ≤ 6` reduce to `x + y ≤ 3`. Negative ratios reverse the
sides. Hashes propose candidates; coefficient comparisons verify proportionality
using `numerics.parallel`. A contradiction certificate requires exact
proportionality and a gap beyond the feasibility margin. Creating a new equality
from inequalities also requires exact proportionality.
Source: [parallel.rs](src/core/rules/parallel.rs).

**Parallel columns — `parallel_columns`.** Suppose the shared constraint columns
satisfy `M[:, k] = r M[:, j]`, where `M` includes both linear and conic rows.
The rule requires exact proportionality and the exact curvature relation
`P[:, k] = r P[:, j]`. It then considers:

- **Aggregation:** if `cₖ = r cⱼ`, constraints and objective depend only on
  `v = xⱼ + r xₖ`. Replace the pair with `v`, whose bounds are the interval sum
  of the original bounds with sign-aware scaling. For instance, two identical
  columns with bounds `[0, 2]` and `[1, 3]` aggregate to bounds `[1, 5]`.
- **Objective dominance:** if the costs differ, movement along
  `(Δxⱼ, Δxₖ) = (-r d, d)` leaves constraints and curvature unchanged while
  changing the linear objective. If one variable can move indefinitely in the
  improving direction, fix the other at its limiting bound. If both can move
  indefinitely, produce an improving recession ray.

Approximate column relations are insufficient for either reduction, even when
they pass the initial candidate comparison. Postsolve splits an aggregate value
back into variables satisfying their original bounds.
Source: [parallel.rs](src/core/rules/parallel.rs).

### Matrix sparsity

**Row sparsification — `sparsification`.** Search for substantial coefficient
cancellation between linear rows, using reference rows with at least 10 nonzeros.
For an equality `aᵀx = b`, replace a target expression by
`a_targetᵀx - α aᵀx` and shift its bounds by `-αb`.

For a reference inequality or ranged row `l ≤ aᵀx ≤ u`, introduce an activity
variable `t` with bounds `[l, u]` and the equality `aᵀx - t = 0`. Target rows can
then use `(a_target - αa)ᵀx + αt`, preserving their original bounds. This may add
a variable while reducing the total coefficient count in a conic export.

The rule requires each rewritten target to be shorter and a positive total
saving after accounting for the activity variable and finite bound sides. It
limits search work, coefficient growth, and scaling; rejects overflow and tiny
nonzero cancellation residuals; and does not change `P`. Cleanup runs afterward
to exploit any new singletons or bound implications.
Source: [sparsification.rs](src/core/rules/sparsification.rs).

### Cone geometry

All reductions below share the **`cones`** switch. A coordinate is
**structurally zero** only when its constraint row is empty and its right-hand
side is exactly zero, so its slack is identically zero. Near-zero values do not
qualify. Sources: [cones.rs](src/core/rules/cones.rs),
[membership and separation](src/problem/cone.rs).

| Cone or pattern | Implemented reduction |
| --- | --- |
| Constant recognized block | When every row is empty, test `h ∈ K`. Delete a block classified inside; return an infeasibility certificate when a sufficiently separated dual witness is found. An inconclusive classification does not itself remove the block. |
| Zero cone | Convert every coordinate to a linear equality `Gᵢx = hᵢ`. |
| Nonnegative cone | Convert every coordinate to an upper-sided linear row `Gᵢx ≤ hᵢ`. |
| Second-order cone: negative constant head | Detect infeasibility when the constant head is below `-numerics.feasibility`. |
| Second-order cone: zero head | Since `0 ≥ ‖z‖₂` forces every tail coordinate to zero, replace the tail by linear equalities and remove the head. |
| Second-order cone: zero tail coordinates | Remove structurally zero tail coordinates and shrink the cone. |
| Second-order cone: dimension 1 or 2 | Convert a lone head to nonnegativity. Convert `(t, z)` with `t ≥ |z|` into the two linear inequalities `t + z ≥ 0` and `t - z ≥ 0`. |
| Rotated second-order cone | Remove structurally zero tail coordinates. If only the two heads remain, replace the cone by their two nonnegativity constraints. |
| PSD cone: all diagonals structurally zero | PSD then forces every off-diagonal to zero. Replace off-diagonals by linear equalities and remove the diagonal coordinates. |
| PSD cone: independent principal blocks | Split along connected components of structurally nonzero off-diagonal expressions. Keep larger components as smaller PSD cones, turn scalar components into nonnegativity constraints, and remove zero coordinates between components. |
| Exponential and power cones | Apply constant-block membership/separation only; there are no additional structural reductions. |
| Opaque cones | Preserve the caller's geometry, identifier, and coordinates; do not attempt membership tests or cone-specific simplification. |

PSD coordinates use upper-triangular, column-major `svec` storage with
off-diagonals multiplied by `√2`. Partial PSD zero faces are deliberately retained:
removing them can prevent finite dual recovery. The all-zero-diagonal face has a
dedicated finite dual reconstruction. Rotated cones do not have a corresponding
zero-head face reduction in this implementation.

## Scheduling and numerical controls

The [scheduler](src/core/schedule.rs) runs rules in the following order. The
thresholds below describe the default configuration:

1. **Cleanup to stability:** fixed variables, cones, empty columns, simple dual
   fixing, singleton rows, empty rows, and cones again.
2. **Fast exploration:** singleton columns, doubleton equalities, and short
   equalities, with cleanup between groups. Repeat fast phases while each reduces
   `nnz(A) + nnz(G) + nnz(P)` by more than 5%.
3. **Medium exploration:** bound propagation, coupled dual fixing, short
   equalities, parallel rows, and parallel columns, interleaved with cleanup.
   Propagation permits up to three extra rounds subject to work and time limits.
   Start another fast/medium cycle only if the completed cycle reduced the same
   nonzero measure by more than 5%.
4. **Final passes:** row sparsification and its follow-up cleanup, then redundant
   variable-bound removal, provided the run has not reported a time limit.

The progress threshold and bounded searches mean presolve need not exhaust every
possible reduction. For this measure, Hessian off-diagonal entries count twice. `Progress::AnyChange`
continues fast phases and full cycles after any model edit, including bound-only
changes or substitutions that increase nonzeros. Rule-specific limits and the
time budget can still prevent further reductions.

| Setting | Default | Meaning |
| --- | --- | --- |
| `time_limit` | 60 seconds | Per-call soft budget including model construction and export preparation, excluding `Presolver` initialization; an in-progress transaction can finish before the limit is observed. |
| `substitution_fill` | `64` | Maximum new constraint and Hessian coefficients per substitution; `usize::MAX` removes the cap. |
| `allow_hessian_growth` | `false` | Permit a net increase in Hessian nonzeros; the fill cap still applies. |
| `equalities.relative_pivot` | `1.0` | Minimum coefficient magnitude divided by the row maximum; invalid values use 1.0. |
| `equalities.max_pivot_attempts` | `1` | Maximum alternative transactions per row and pass; 0 disables them. |
| `equalities.cost_aware` | `false` | Prefer estimated A and P work after free-variable status, instead of column length. |
| `equalities.max_row_length` | `8` | Maximum equality row length for the generalized substitution rule. |
| `equalities.min_column_length` / `max_column_length` | `2` / `8` | Candidate pivot column degree range; use `1` / `usize::MAX` to admit all nonempty columns. |
| `equalities.require_free_variable` | `true` | Require a variable with no explicit bounds. |
| `equalities.require_linear_variable` | `true` | Require a pivot absent from the Hessian. |
| `equalities.preserve_nonzeros` | `true` | Cap new coefficients by removed constraint entries, preventing net growth in total constraint and Hessian nonzeros. |
| `equalities.work_limit` | `Default` | Estimated entry visits per pass; default is twice the current constraint nonzeros. |
| `propagation.minimum_relative_gain` | `0.01` | Relative threshold for finite non-fixing bound changes; aggressive uses `0.005`. |
| `propagation.minimum_gain_factor` | `1e4` | Absolute threshold multiplier applied to `numerics.feasibility`. |
| `propagation.additional_rounds` | `3` | Extra propagation rounds after the initial pass; `usize::MAX` removes the cap. |
| `propagation.work_limit` | `Default` | Work allowance across extra rounds; default is `max(constraint nonzeros / 4, 256)`. |
| `progress` | `Nonzeros { minimum_reduction: 0.05 }` | Fractional nonzero decrease required to continue, or `AnyChange` to continue after any edit. |
| `sparsification.allow_auxiliary_variables` | `true` | Allow inequality references that introduce activity variables; `false` uses equality references only. |
| `sparsification.work_limit` | `Default` | Work allowance per pass; default is `max(8 * (constraint nonzeros + bound entries), 1024)`. |
| `threads` | `1` | Thread count for fingerprinting and sorting: `1` is serial, `0` lets Rayon select automatically. |
| `numerics.feasibility` | `1e-9` | Relative feasibility margin for separation checks; it does not round bounds or cone coordinates to zero. |
| `numerics.parallel` | `1e-12` | Relative coefficient tolerance for parallel candidates; operations requiring exact relations still enforce them. |
| `numerics.huge_bound` | `1e7` | Propagation skips candidate bounds with absolute value at least this large. |

For propagation to improve an already finite bound, the gain must exceed
`max(minimum_gain_factor * feasibility, minimum_relative_gain * abs(old_bound))`,
using the propagation settings (defaults `1e4` and `0.01`), except when the new bound exactly
meets the opposite bound. Singleton-row bound extraction does not use this gain
filter or `huge_bound`. Numerical scaling and cancellation safeguards remain
internal constants. Each `WorkLimit` accepts `Default`, `Entries(n)`, or
`Unlimited`. These work budgets are independent of `time_limit`. The time limit
is soft: matrix preparation, an in-progress operation, and result packing may
extend elapsed time beyond it; it is not a memory limit.

Select rule families through `Settings`, for example:

```rust
use presolve::{Settings, settings::Rules};

let settings = Settings {
    rules: Rules {
        fixed_variables: true,
        empty_columns: true,
        empty_rows: true,
        singleton_rows: true,
        ..Rules::none()
    },
    ..Settings::default()
};
```

### Aggressive dimensional reduction

Use the measured aggressive preset when fewer variables and constraints matter
more than matrix sparsity:

```rust
use presolve::Settings;
use std::time::Duration;

let settings = Settings::aggressive(Duration::from_secs(2));
```

This enables unrestricted substitution fill and Hessian growth, admits bounded
and quadratic equality pivots, continues after any model edit, and removes the
extra propagation round and work caps. It enables bounded exact equality
dependency checks and lowers the relative propagation gain threshold to 0.005.
It uses equality-only sparsification and keeps the default numerical tolerances
and one-thread execution.

Candidate equality rows and pivot columns are limited to 16 entries. These are
search limits, not fill limits: a selected substitution can introduce any number
of coefficients. In the Netlib and Maros–Mészáros comparison, this restriction
produced more dimensional reduction within the budget than admitting every
candidate. Set `equalities.max_row_length` and `equalities.max_column_length` to
`usize::MAX` to remove these limits too. All preset fields remain editable.

The preset is a starting point, not a guarantee of the smallest model on every
problem. Use the benchmark commands below to compare dimensions, finite bound
sides, matrix nonzeros, runtime, and time-budget outcomes on your models.

Equality substitution also rejects near-cancellation when a nonzero updated
coefficient is smaller than `1e-10` times the larger contributing term. This
prevents roundoff from becoming a later pivot. Exact zeros remain valid; rejected
substitutions leave the model unchanged.

### Optional parallel execution

```rust
use presolve::{InitError, Presolver, Problem, Settings};

fn process(problems: impl IntoIterator<Item = Problem>) -> Result<(), InitError> {
    let presolver = Presolver::new(Settings {
        threads: 4,
        ..Settings::default()
    })?;
    for problem in problems {
        let result = presolver.presolve(problem);
        // Use result.outcome and result.stats here.
    }
    Ok(())
}
```

The first parallel implementation computes row/column fingerprints and sorts
candidate keys concurrently. Proportionality checks, reductions, cleanup, and
postsolve record creation remain serial. Sorting uses original indices to break
ties, preserving serial candidate and reduction order. Runs that stop on the
wall-clock limit can still differ in how much work they complete.

Set `threads` to `1` for serial execution (the default), to a larger number for
that many worker threads, or to `0` for Rayon's automatic selection, which honors
`RAYON_NUM_THREADS`. `Presolver::new` creates a dedicated pool and returns
`InitError` if initialization fails. With `threads: 1`, it creates no pool.
An internal executor selects ordinary serial iterators or work in the owned
Rayon pool. There is no global pool cache, and Rayon's global pool is unaffected.

The presolver retains its settings and pool across calls; `threads()` reports
the effective thread count. Each call has its own working model and recovery
tape, so calls may run concurrently through `&self`, and results can outlive the
presolver. The time budget and elapsed statistic reset for each call and exclude
pool initialization.

Scans currently stay serial below 1,024 row/column slots or 32,768 estimated
nonzeros, or when the resolved thread count is one. These initial thresholds
limit scan overhead; pool initialization still happens in `Presolver::new`, even
for small problems. Enabling the option does not guarantee a speedup. The scan
estimates use constraint nonzeros for rows and constraint plus Hessian nonzeros
for columns.

The benchmark accepts `--threads N` and `--pool-mode cold|reused` for both time
and size runs. `--preset default|fill|aggressive|unrestricted` selects the baseline,
unrestricted fill with baseline searches, the measured aggressive configuration,
or that configuration without equality length limits. Use `--time-limit-ms N`,
`--equality-row-limit N`, `--equality-column-limit N`, and
`--sparsification all|equalities|off` to override the preset. Pivot experiments
also accept `--equality-pivot-relative R`, `--equality-pivot-attempts N`, and
`--equality-cost-aware true|false`. Metadata records the
complete resulting settings. Propagation experiments accept
`--propagation-relative-gain R` and `--propagation-gain-factor F`; invalid
(nonfinite or negative) values use the corresponding default. Measurements include finite variable-bound
side counts. For example, to compare repeated-call performance:

```sh
cargo run --release -p benchmark -- run time --name serial --problem QAP12 --rule all --trials 5 --pool-mode reused
cargo run --release -p benchmark -- run time --name parallel-4 --problem QAP12 --rule all --trials 5 --threads 4 --pool-mode reused
cargo run --release -p benchmark -- compare time serial parallel-4
```

Timing trials use fresh processes. The default `cold` mode includes presolver
initialization. `reused` initializes the presolver and runs one untimed warm-up
on the same problem before measuring, so it also benefits from warmed caches.
Loading and destruction are outside timing in both modes. Timing comparisons
require matching pool modes, which are recorded in result metadata.

## Results and postsolve

`presolver.presolve(problem)` consumes the problem and returns an outcome plus
size and execution statistics. The one-shot helper `presolve(problem, &settings)`
creates temporary execution resources and returns `Result<PresolveResult,
InitError>`; use `presolve(problem, &settings)?` to propagate initialization
failures. These errors are separate from optimization outcomes:

| Outcome | Meaning |
| --- | --- |
| `Unchanged` | No rule was applied; returns the input problem. |
| `Reduced` | Returns a smaller or transformed problem and its `Postsolve` map. |
| `Solved` | Presolve eliminated the full problem and returns an original-coordinate solution. Remaining opaque cones prevent this outcome. |
| `Infeasible` | Returns original-coordinate Farkas multipliers proving a contradiction. |
| `Unbounded` | Returns both a feasible point and an objective-decreasing recession ray. |

Finding a recession ray alone does not prove primal unboundedness. The entry
point also checks a cheap feasible-point candidate; if that check fails, it
returns a reduced or unchanged problem for a solver to handle.

Applied transformations record recovery data on a tape. Reversing the tape
restores eliminated variables, linear and bound multipliers, conic duals, and
conic slacks. `Postsolve` exposes solution recovery, primal-ray recovery,
infeasibility-certificate recovery, and warm-start reduction. Native dual signs
satisfy `P x + c - Aᵀ y - z + Gᵀ w = 0`: linear and bound multipliers are positive
on lower sides and negative on upper sides. Recovery does not clip values or
adjust them for interiority. Forward warm starts can need solver refinement
after redundant constraints are removed; they are not guaranteed to remain
optimal or stationary.

Reduced constraint storage is exported explicitly with `Problem::into_csc()` or
`Problem::into_conic()`. The latter returns an additional map for translating
conic-form multipliers into native coordinates before postsolve.
See [results](src/result.rs), [postsolve](src/postsolve/mod.rs), and
[conic export](src/problem/conic.rs).

## Development

Requires Rust 1.86 or newer. Run the library tests and documentation examples with:

```sh
cargo test -p presolve
```

Licensed under [Apache-2.0](LICENSE). Attribution is recorded in [NOTICE](NOTICE).
