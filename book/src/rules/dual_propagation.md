# Dual propagation

Dual propagation derives bounds on dual multipliers, uses their derivations to propose coordinated movements of primal variables, and verifies each movement before applying a reduction. It can fix a variable at a bound or restrict a one-sided linear inequality to its finite side as an equality.

Unlike [dual fixing](./dual_fixing.md), the movement can involve several variables. Multiplier bounds suggest a reduction; a verified primal direction supplies its justification. The argument does not require an optimal solution to exist before the reduction is applied.

## Coordinated movement example

Consider

\\[
\min_{x,y}2x+y,
\qquad x+y=5,
\qquad x,y\ge0.
\\]

The equality blocks either variable from moving alone. But the transfer

\\[
x^{\prime}=x-t,\qquad y^{\prime}=y+t,\qquad 0\le t\le x
\\]

preserves the equality and nonnegativity while decreasing the objective by \\(t\\). Every feasible point can slide to \\(x=0\\). Presolve can fix that variable and then use ordinary substitution and cleanup.

For the equality multiplier \\(\lambda\\), the nonnegative variables imply \\(\lambda\le2\\) and \\(\lambda\le1\\). Thus the signed bound multiplier of \\(x\\) is \\(2-\lambda\ge1\\), suggesting its lower-bound fixing. The transfer above verifies the suggestion directly.

## Dual restrictions

Use the [shared notation and multiplier signs](../problem_formulation.md#recovery-and-optimality). The rule uses columns without quadratic participation or incident cone rows. For such a variable, stationarity gives

\\[
c_j=\sum_i a_{ij}y_i+z_j.
\\]

The variable's interval supplies these restrictions:

| Variable bounds | Restriction on the column's multiplier sum |
| --- | --- |
| Lower bound only | \\(\sum_i a_{ij}y_i\le c_j\\) |
| Upper bound only | \\(\sum_i a_{ij}y_i\ge c_j\\) |
| No finite bounds | \\(\sum_i a_{ij}y_i=c_j\\) |
| Two finite bounds | No such restriction from this column alone |

The initial multiplier intervals use the signs of the linear row sides. A lower-only row starts with \\(y_i\ge0\\), and an upper-only row with \\(y_i\le0\\). Rows with two finite sides, including equalities, start without a sign restriction. Rows with neither side finite are also left unrestricted by this initialization rather than contributing a tighter interval.

Propagation isolates one multiplier at a time in these column constraints using the other multiplier intervals, much like [bound propagation](./bound_propagation.md) isolates a primal variable. Coefficient signs determine the bound direction. Updated intervals can enable further implications.

## From multiplier bounds to reductions

Two kinds of conclusions are considered:

- A one-sided row multiplier bounded strictly away from zero suggests that its finite side can be made active.
- A signed variable-bound multiplier bounded strictly positive or negative suggests fixing the variable at its lower or upper bound, respectively.

The derivation of each tightened multiplier bound is retained. Combining the column restrictions used in that derivation produces a candidate primal direction \\(d\\). Before using it, presolve checks the direction against the current problem, not just the propagated intervals.

### Direction verification

The direction must have finite components, involve only eligible linear columns, and satisfy

\\[
c^Td\le0.
\\]

Because the moving variables have no quadratic participation, \\(Pd=0\\), so the objective does not increase anywhere along the movement.

Every equality and ranged row must remain unchanged. Other non-target restrictions must remain feasible:

- An equality or ranged row requires \\((Ad)_i=0\\).
- A lower-only row requires \\((Ad)_i\ge0\\).
- An upper-only row requires \\((Ad)_i\le0\\).
- A non-target variable with two finite bounds cannot move; a one-sided variable can move only away from its finite bound.

The designated target is allowed to move toward its finite endpoint. A target row must be a one-sided inequality and move toward its finite side; a target variable must move toward the chosen bound. Movement stops there before violating the target. The other restrictions permit the whole required movement from any feasible starting point.

These sign and structural checks use direct floating-point comparisons, not a tolerance that permits an increasing objective or a nonzero change in an equality. They are numerical checks, not symbolic real-arithmetic proofs.

### Restricting an inequality

For example, minimize \\(x+y\\) subject to \\(x+y\ge5\\), with both variables free. Their dual constraints imply that the row multiplier is \\(1\\). Decreasing \\(x\\) while leaving \\(y\\) fixed lowers the objective until the row reaches \\(x+y=5\\).

The inequality can therefore become an equality. This discards feasible points with larger objective values while retaining enough points to recover an optimum. It can expose equality substitutions for subsequent cleanup.

## Contradictions and numerical safeguards

A contradiction between propagated multiplier intervals does not directly prove primal infeasibility. The pass abandons its proposed reductions when those intervals become inconsistent. A reduction is accepted only after a primal direction is independently verified.

Proposed multiplier bounds must be finite and tighter than the current interval. Conflicts use the configured feasibility margin; a smaller conflict rejects the candidate. Changes to already finite bounds must also exceed the propagation gain threshold described in [Settings](../settings.md#bound-propagation).

If a verified direction does not drive its proposed target but is feasible indefinitely and has strictly negative objective slope, it can supply an improving recession certificate instead. Such a direction establishes unbounded improvement from a feasible point, not feasibility by itself. A zero-cost direction without a reached target does not establish unboundedness.

## Optimality and recovery

For an accepted fixing or row restriction, any original feasible point can slide to the selected bound or row side without increasing its objective. In exact arithmetic, the reduction preserves the optimal value and permits recovery of an original optimum when one exists. It need not preserve every original feasible point.

Finite fixings use [fixed-variable recovery](./fixed_variables.md#recovery). Restricting a row to a side keeps the same variable coordinates and row coefficients. The reduced equality's multiplier already has a sign valid for that original side under the verified direction argument; no separate coordinate transformation is required for that restriction. Earlier reductions are still reversed during postsolve.

## Configuration and scheduling

`rules.dual_propagation` defaults to `false`; the aggressive preset enables it. It is separate from `rules.dual_fixing`. The rule runs once among the final passes, after the redundant-bound-removal stage, when those enabled stages are reached within the soft time budget. Bound removal can expose free or one-sided variables that provide stronger dual restrictions.

`dual_propagation.work_limit` defaults to four times the current constraint nonzero count. Propagation and direction extraction each use an allowance of that size, limiting both implication work and long proof reconstruction. A work-limited pass may still use partial conclusions if they remain consistent, but every applied direction must pass verification on the current model.

Successful reductions trigger local cleanup and equality substitutions rather than restarting the full exploration cycle. See [Settings](../settings.md) and [Rule Scheduling](../scheduling.md) for the surrounding controls and order.
