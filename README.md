# presolve

Sparse presolve and postsolve for convex quadratic and conic optimization in Rust.
The library simplifies a problem before it reaches a solver, then maps solutions
and certificates back to the original variables and constraints.

This README catalogs the implemented rules, their applicability, and how they
interact. All 14 rule families are enabled by default.

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
| Substitution | `short_equalities` | Eliminate a free variable from an equality with 3–8 nonzeros. |
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
nonfinite arithmetic, and forbid a net increase in Hessian nonzeros.
Sources: [substitution.rs](src/core/rules/substitution.rs),
[model.rs](src/core/model.rs), [objective.rs](src/core/objective.rs).

**Singleton columns — `singleton_columns`.** Consider a variable appearing in
exactly one row of the shared constraint matrix, where that row is linear and
has at least two nonzeros.

- In an equality, substitute the variable after discarding bounds implied by
  the equality and the other variables' bounds. If both bound sides remain,
  defer to the doubleton rule. A variable appearing in `P` is also deferred if
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

**Short equalities — `short_equalities`.** Extend equality substitution to rows
with 3–8 nonzeros. The pivot must be free, absent from `P`, occur in 2–8 constraint
rows, and have a coefficient of maximum magnitude in the equality. Among
candidates, prefer the shorter column. The rule removes the pivot and its
equality without introducing a bound row or changing `P`. Its fill allowance is
further capped by the removed row and column, preventing a net increase in
constraint nonzeros, and exploration has a separate work budget and deadline.

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

The [scheduler](src/core/schedule.rs) runs rules in the following order:

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
possible reduction. For this measure, Hessian off-diagonal entries count twice.

| Setting | Default | Meaning |
| --- | --- | --- |
| `time_limit` | 60 seconds | Soft budget including construction and export preparation; an in-progress transaction can finish before the limit is observed. |
| `substitution_fill` | `64` | Maximum new constraint and Hessian coefficients per equality substitution, with a separate prohibition on net Hessian growth. |
| `numerics.feasibility` | `1e-9` | Relative feasibility margin for separation checks; it does not round bounds or cone coordinates to zero. |
| `numerics.parallel` | `1e-12` | Relative coefficient tolerance for parallel candidates; operations requiring exact relations still enforce them. |
| `numerics.huge_bound` | `1e7` | Propagation skips candidate bounds with absolute value at least this large. |

For propagation to improve an already finite bound, the gain must exceed
`max(1e4 * feasibility, 0.01 * abs(old_bound))`, except when the new bound exactly
meets the opposite bound. Singleton-row bound extraction does not use this gain
filter or `huge_bound`. Rule-specific scaling and work limits are currently
internal constants.

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

## Results and postsolve

`presolve(problem, &settings)` consumes the problem and returns an outcome plus
size and execution statistics:

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
adjust them for interiority.

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
