# Dominated-column elimination

Dominated-column elimination fixes a variable at a bound when transferring its contribution to another variable preserves feasibility and does not increase the objective. It removes the fixed variable using the shared fixed-value substitution machinery.

Unlike parallel-column aggregation, this rule does not require proportional constraint columns. The permitted transfer depends on each constraint's direction. This implementation considers variables with no quadratic objective terms, including cross terms.

## Example

Consider

\\[
\min_{x,y} 2x+y,
\qquad x+y\ge5,
\qquad x,y\ge0.
\\]

Increasing \\(y\\) while decreasing \\(x\\) by the same amount leaves the constraint unchanged and reduces cost. Since \\(y\\) has no upper bound, this transfer can continue until \\(x = 0\\). Presolve fixes \\(x\\) there, leaving the objective \\(y\\) and the row \\(y \ge  5\\).

For example, \\(x = 2\\), \\(y = 3\\) has cost \\(7\\); replacing it by \\(x = 0\\), \\(y = 5\\) has cost \\(5\\) and the same constraint contribution.

If \\(y \le  3\\) is added, this argument no longer permits fixing \\(x = 0\\): satisfying the row requires \\(x \ge  2\\). The receiving variable needs enough room to accept the transfer.

## Mathematical formulation

Let \\(x_j\\) be the proposed dominant variable and \\(x_k\\) the proposed dominated variable. Consider the transfer

\\[
x_j^{\prime}=x_j+t,
\qquad x_k^{\prime}=x_k-t,
\qquad t\ge0.
\\]

Other variables are unchanged. With no quadratic terms involving either variable, the objective change is

\\[
f(x^{\prime})-f(x)=t(c_j-c_k).
\\]

Thus \\(c_j\le c_k\\) ensures that the transfer does not increase the objective. Equal costs are allowed for a bound fixing; strict improvement is required for an improving recession certificate.

### Constraint directions

For a linear row \\(i\\), its expression changes by \\(t d_i\\), where

\\[
d_i=A_{ij}-A_{ik}.
\\]

A missing coefficient is treated as zero. For the transfer to preserve every feasible row value, the implementation requires:

| Row bounds | Required coefficient difference | Reason |
| --- | --- | --- |
| Finite lower and upper bounds | \\(d\_i = 0\\) | Preserve the row expression exactly. |
| Upper bound only | \\(d\_i \le  0\\) | Do not increase the row expression. |
| Lower bound only | \\(d\_i \ge  0\\) | Do not decrease the row expression. |
| Neither side finite | Any finite \\(d\_i\\) | No row restriction to violate. |

The first case includes equalities and ranged rows. Non-linear row domains encountered in the pair's support reject the relation.

For example, transferring from \\(x\\) to \\(y\\) preserves \\(x + y \ge  5\\) and increases the left-hand side of \\(x + 2y \ge  6\\), so both lower-bound constraints remain satisfied. The same transfer could violate \\(x + 2y \le  6\\) and would fail that row's check.

## Open bound sides and fixing targets

The transfer increases the dominant variable and decreases the dominated variable. Their opposite bound sides cannot be violated by this movement, but the dominant upper side and dominated lower side can stop it.

A bound side is treated as open when it is infinite or a retained linear row already implies that variable bound using the other variables' current bounds. An implied side remains satisfied because the transfer preserves the rows and the other necessary variable bounds; the proving row is retained for recovery.

After proving dominance, the implementation selects a target in this order:

1. If the dominated variable has a finite lower bound and the dominant upper side is open, fix the dominated variable at its lower bound.
2. Otherwise, if the dominant variable has a finite upper bound and the dominated lower side is open, fix the dominant variable at its upper bound.
3. If both directional endpoints are infinite and costs differ, return an improving recession certificate.

If neither fixing is justified, the pair is left unchanged. When both directional endpoints are infinite but costs are equal, the code also leaves the pair unchanged rather than claiming unbounded improvement.

Thus the rule does not always remove the variable called dominated: it can fix the dominant variable at its upper bound instead.

### Optimality and recession direction

For an accepted fixing, every original feasible point can be transferred to a feasible point with that fixed value and no larger objective. Substituting the value therefore preserves the optimal value and permits recovery of an original optimum when one exists. It can discard other feasible points, so the reduction is not a one-to-one feasible-set correspondence.

When the dominant upper bound and dominated lower bound are infinite and \\(c_j<c_k\\), the direction is

\\[
d=e_j-e_k,\qquad c^Td=c_j-c_k<0.
\\]

The row-direction checks make this a feasible recession direction for the constraints, and neither variable has a bound blocking it. Since neither variable participates in the quadratic objective, \\(Pd=0\\). This establishes unbounded improvement from a feasible point; it does not prove that the rest of the problem is feasible.

## Applicability

Dominance uses direct cost and coefficient-difference comparisons, without a parallel tolerance. Nonfinite differences reject the relation. Both variables must have no quadratic objective terms or incident cone rows. The search is bounded rather than exhaustive; [Settings](../settings.md#dominated-columns) describes enabling it and widening its scope.

## Recovery

Postsolve restores the chosen fixed value and reconstructs its signed bound multiplier from stationarity. If a finite implied bound on the other variable justified the transfer, its multiplier can first move to the proving row. This accounts for that implied restriction before recovering the eliminated variable's multiplier.

Infinite open sides need no implied-bound proof. A returned recession direction is mapped through earlier reductions to the original coordinates.

See [Settings](../settings.md) for rule switches and tuning, and [Rule scheduling](../scheduling.md) for ordering and stopping conditions.
