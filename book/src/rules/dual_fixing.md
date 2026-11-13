# Dual fixing

Dual fixing places a variable at a bound when movement toward that bound cannot violate a constraint and cannot increase the objective. The chosen value is then substituted as in [fixed-variable elimination](./fixed_variables.md).

The variable need not already have equal bounds, and it can appear in several constraints. The proof combines a safe direction of movement with the sign of the objective derivative. It does not require first solving a separate dual problem.

## Example

Consider

\\[
\min_{x,y}3x+2y,
\qquad x+y\le10,
\qquad x,y\ge0.
\\]

Decreasing \\(x\\) lowers the objective and makes the upper-bound constraint easier to satisfy. Every feasible point can therefore be moved to \\(x=0\\) without increasing its cost. Fixing \\(x=0\\) leaves the objective \\(2y\\) and the constraints \\(0\le y\le10\\).

By contrast, a row \\(x+y\ge10\\) can block decreasing \\(x\\). A favorable objective coefficient alone is not enough to justify fixing.

For an upper-bound example, minimizing \\(-3x+2y\\) subject to \\(x+y\ge10\\), \\(0\le x\le8\\), and \\(y\ge0\\) allows increasing \\(x\\) without worsening the objective or violating the row. Presolve can fix \\(x=8\\).

## Safe directions and constraint locks

Use the [shared problem notation](../problem_formulation.md). A constraint locks a direction if moving a variable in that direction could violate one of its finite row bounds. A lock is conservative: the existence of a lock does not mean every particular move is infeasible, but the rule cannot assume unrestricted movement in that direction.

For a nonzero coefficient \\(a_{ij}\\), the locks are:

| Row restriction | Positive coefficient | Negative coefficient |
| --- | --- | --- |
| Finite upper side | Locks increasing the variable | Locks decreasing the variable |
| Finite lower side | Locks decreasing the variable | Locks increasing the variable |
| Two finite sides, including an equality | Locks both directions | Locks both directions |

A row with neither side finite imposes no lock. An incident cone row conservatively locks both directions. A variable's own bounds are handled separately: they determine the endpoints toward which it may move.

If no row locks decreasing \\(x_j\\), lowering it toward its lower bound preserves the constraints. The analogous statement holds for increasing toward the upper bound. Other variables stay fixed during this argument.

## Independent objective terms

First suppose \\(x_j\\) has no quadratic cross terms. Its objective contribution and derivative are

\\[
q_j(t)=\tfrac12 p t^2+c_jt,
\qquad q_j^{\prime}(t)=pt+c_j,
\qquad p=P_{jj}\ge0.
\\]

Because the derivative is nondecreasing, two endpoint tests suffice:

- If decreasing is unlocked, \\(L_j\\) is finite, and \\(pL_j+c_j\ge0\\), fix \\(x_j=L_j\\). The objective is nondecreasing throughout its interval, so decreasing the variable cannot hurt.
- If increasing is unlocked, \\(U_j\\) is finite, and \\(pU_j+c_j\le0\\), fix \\(x_j=U_j\\). The objective is nonincreasing throughout its interval, so increasing the variable cannot hurt.

The endpoint derivative must also be finite. Negative diagonal curvature is skipped by this rule. For a linear term, \\(p=0\\), so these tests reduce to checking the sign of \\(c_j\\).

Unlike [empty-column elimination](./empty_columns.md), this rule does not generally choose an interior quadratic minimizer. Its finite-value reductions use safe bound endpoints.

## Coupled quadratic terms

When cross terms are present, the derivative depends on other variables:

\\[
g_j(x)=P_{jj}x_j+c_j+\sum_{k\ne j}P_{jk}x_k.
\\]

To justify a lower-bound fixing, bound this derivative from below at \\(x_j=L_j\\), using every other variable's interval. A positive coefficient uses the other variable's lower bound, and a negative coefficient uses its upper bound:

\\[
\underline g_j(L_j)=P_{jj}L_j+c_j+\sum_{k\ne j}\min(P_{jk}L_k,P_{jk}U_k).
\\]

If this lower estimate has a reliably positive sign, decreasing is unlocked, and \\(P_{jj}>0\\), then decreasing to \\(L_j\\) cannot increase the objective for any feasible values of the other variables. For an upper-bound fixing, use the upper derivative estimate

\\[
\overline g_j(U_j)=P_{jj}U_j+c_j+\sum_{k\ne j}\max(P_{jk}L_k,P_{jk}U_k),
\\]

and require a reliably negative sign with increasing unlocked. Sums range over nonzero cross coefficients.

For example, with objective \\(x^2+xy+y^2+3x\\), bounds \\(x\ge0\\), \\(0\le y\le2\\), and row \\(x+y\le10\\), decreasing \\(x\\) is unlocked. At \\(x=0\\), its derivative is \\(y+3\ge3\\), so \\(x=0\\) is an optimal choice for every feasible fixed value of \\(y\\).

The coupled test is restricted to Hessian rows with 2–33 stored entries and positive diagonal curvature. It requires a finite derivative estimate separated from zero by more than \\(64\\) times machine epsilon times the sum of absolute contributing terms, including the linear coefficient. Unbounded or cancellation-sensitive estimates can therefore prevent a fixing even when a mathematical argument exists.

## Missing endpoints and zero cost

For a linear term with \\(c_j>0\\), no downward locks and no lower bound give the improving direction \\(-e_j\\). For \\(c_j<0\\), no upward locks and no upper bound give \\(e_j\\). These are recession directions: they establish unbounded improvement from a feasible point, not feasibility of the rest of the problem.

When \\(p=c_j=0\\), movement has no objective cost. A finite safe endpoint can still be chosen. If the safe direction has no finite endpoint, the variable and its incident rows can instead be removed by projecting out the variable. The missing endpoint allows it to satisfy those one-sided rows for any retained point satisfying the remaining restrictions. This is not a fixing to infinity and does not establish unbounded improvement.

## Optimality and recovery

For an accepted finite fixing, every original feasible point can be moved to the chosen endpoint without increasing the objective. Substitution preserves the optimal value and permits recovery of an original optimum when one exists. Other feasible points can be discarded, so the correspondence need not preserve every original point.

Postsolve restores the chosen finite value and recovers its signed bound multiplier from stationarity, as in [fixed-variable recovery](./fixed_variables.md#recovery).

For a zero-cost variable projected out through an infinite safe direction, postsolve instead intersects its original interval with the restrictions imposed by the saved rows at the recovered remaining variables. It chooses a finite value in that interval by clipping zero to it. The corresponding row and variable-bound multipliers are recovered as zero. Recession directions are also mapped back through earlier reductions.

## Configuration

The `rules.dual_fixing` switch enables both the simple and coupled tests and defaults to `true`. The simple test runs during cleanup; the coupled test runs in medium phases. The separate `rules.dual_propagation` switch controls a different rule. See [Settings](../settings.md) and [Rule Scheduling](../scheduling.md) for configuration and ordering.

Numerical rejection of a fixing leaves the variable available to the solver or other rules; it does not establish infeasibility.
