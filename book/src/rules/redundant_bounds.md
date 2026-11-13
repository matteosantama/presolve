# Redundant-bound removal

Redundant-bound removal deletes an explicit variable bound when a retained linear constraint and the other variables' current bounds already guarantee it. It keeps the variables, constraint coefficients, and objective unchanged.

Unlike [bound propagation](./bound_propagation.md), which tightens bounds, this rule removes bound sides that add no restriction. Unlike [implied-free equality elimination](./implied_free_equalities.md), it does not substitute a variable or delete its defining row.

## Example

Suppose

\\[
x+y\le5,\qquad y\ge2,\qquad x\le10.
\\]

The row and the lower bound on \\(y\\) imply

\\[
x\le5-y\le3.
\\]

Thus \\(x\le10\\) is redundant. Removing that explicit upper bound does not allow larger feasible values of \\(x\\): the retained row and \\(y\ge2\\) still guarantee \\(x\le3\\).

The rule removes the upper bound rather than replacing it by \\(3\\). A lower bound on \\(x\\), if present, is checked separately.

## Mathematical criterion

Use the [shared problem notation](../problem_formulation.md). For a row containing \\(x_k\\), compute the range of the remaining terms from the other variables' bounds:

\\[
\alpha_{i,-k}=\sum_{j\ne k}\min(a_{ij}L_j,a_{ij}U_j),
\qquad
\beta_{i,-k}=\sum_{j\ne k}\max(a_{ij}L_j,a_{ij}U_j).
\\]

Sums range over nonzero coefficients. For \\(a_{ik}>0\\), the row implies

\\[
\widehat L_k=\frac{\ell_i-\beta_{i,-k}}{a_{ik}},
\qquad
\widehat U_k=\frac{u_i-\alpha_{i,-k}}{a_{ik}}.
\\]

For \\(a_{ik}<0\\), division reverses the inequalities:

\\[
\widehat L_k=\frac{u_i-\alpha_{i,-k}}{a_{ik}},
\qquad
\widehat U_k=\frac{\ell_i-\beta_{i,-k}}{a_{ik}}.
\\]

A finite explicit lower bound is redundant if a finite implied lower bound satisfies \\(\widehat L_k\ge L_k\\). An explicit upper bound is redundant if a finite implied upper bound satisfies \\(\widehat U_k\le U_k\\). A single retained linear row must supply the proof; this pass does not search arbitrary combinations of constraints.

## Avoiding circular proofs

Consider

\\[
x=y,\qquad x\ge0,\qquad y\ge0.
\\]

Either lower bound follows from the equality and the other lower bound. But removing both would admit \\(x=y=-1\\).

The implementation checks bounds sequentially using the current problem. After removing \\(x\ge0\\), it updates the affected row activities. The remaining bound \\(y\ge0\\) can no longer be proved from that equality and an unbounded \\(x\\), so it stays. Which bound survives can depend on the visit order.

## Applicability and numerical checks

The pass visits live variables with a nonempty constraint column and at least one finite bound, ordered by column length. For a variable with a finite upper bound, it tries that side first. Variables with only a finite lower bound are checked during the same sweep; lower sides of initially two-sided variables are deferred to a second sweep. Every attempt uses the bounds remaining at that point.

Only linear rows provide implication proofs. Cone rows are skipped, but quadratic participation does not itself exclude a variable: the proof concerns feasibility and leaves the objective unchanged.

A candidate implication must be finite. The redundancy comparison uses direct containment, without a feasibility tolerance that would permit weakening the bound. Cancellation-sensitive residual activities can be recomputed from the other row terms. If no usable proof is found, the explicit bound remains; this does not establish infeasibility or prove that the bound is mathematically necessary.

## Equivalence and recovery

In exact arithmetic, each removed bound is already guaranteed by restrictions that remain. Sequential removal therefore preserves the full feasible set, objective values, and optimal solutions.

No variable coordinates change, so this reduction needs no primal reconstruction. It adds no separate recovery record. A removed bound can have zero multiplier; its former restriction is enforced by the retained constraints and bounds. Earlier reductions still undergo their ordinary postsolve recovery.

## Configuration and scheduling

`rules.redundant_bounds` defaults to `true` and remains enabled in the aggressive preset. The pass runs once after the main exploration phases and the equality-dependency and sparsification stages, if the scheduling state has not marked the time limit reached.

Keeping bounds during exploration lets them assist other reductions. Removing redundant sides afterward avoids unnecessary explicit bounds in the solver's input and can expose free or one-sided variables for later rules, including implied-free equality elimination and dual propagation. This sweep has no separate configurable work allowance or internal deadline check; the surrounding time budget is soft.

See [Settings](../settings.md) and [Rule Scheduling](../scheduling.md) for the surrounding controls and order.
