# Implied-free equality elimination

Implied-free equality elimination substitutes a variable out of a linear equality when that same equality and the other variables' bounds already guarantee its explicit bounds. It removes both the variable and the equality without adding a constraint to preserve those bounds.

The variable need not be unbounded. “Implied-free” means that its explicit bounds impose no additional restriction on the substitution.

## Example

Suppose

\\[
x+y=5,\qquad 1\le y\le3,\qquad 0\le x\le10.
\\]

Solving the equality gives \\(x=5-y\\). The bounds on \\(y\\) imply

\\[
2\le x\le4,
\\]

which already satisfies \\(x\\)'s explicit bounds. Presolve can substitute \\(5-y\\) for \\(x\\) in the objective and every other constraint, then remove \\(x\\) and the defining equality if the substitution checks pass.

For example, an additional row \\(x+2y\le8\\) becomes \\(y\le3\\), and an objective \\(2x+y\\) becomes \\(10-y\\). Postsolve reconstructs \\(x=5-y\\).

If instead \\(x\ge3\\), the implied interval \\([2,4]\\) does not guarantee that lower bound. This rule rejects that pivot. Ordinary equality substitution could preserve the necessary restriction as \\(y\le2\\), but this rule requires every explicit bound to be redundant.

## Mathematical formulation

Write the defining equality as

\\[
a_{ik}x_k+\sum_{j\ne k}a_{ij}x_j=b_i,
\qquad a_{ik}\ne0.
\\]

Using the other variables' bounds, let the residual sum lie in \\([r_{\min},r_{\max}]\\). For a positive pivot coefficient, the implied interval is

\\[
\frac{b_i-r_{\max}}{a_{ik}}
\le x_k\le
\frac{b_i-r_{\min}}{a_{ik}}.
\\]

For a negative pivot coefficient, the interval endpoints reverse. Infinite bounds may leave one or both endpoints infinite.

Every finite explicit lower bound must be no larger than the finite implied lower endpoint, and every finite explicit upper bound must be no smaller than the finite implied upper endpoint. The implementation uses direct comparisons for this redundancy check. An absent explicit bound needs no implication proof.

Once both sides are redundant, substitute

\\[
x_k=\frac{b_i}{a_{ik}}-\sum_{j\ne k}\frac{a_{ij}}{a_{ik}}x_j.
\\]

The objective and other rows are transformed as in [short equality elimination](./short_equalities.md#mathematical-formulation), with no remaining bound row for the eliminated variable. In exact arithmetic, original and reduced feasible points correspond through this reconstruction and have equal objective values.

## Applicability and safeguards

The rule considers linear equalities with at least two nonzero coefficients. It checks the implied interval from the defining equality itself, rather than relying on an implication that could disappear when that equality is removed.

The shared [equality settings](../settings.md#equality-substitution) control pivot stability, candidate attempts, column degree, quadratic participation, cost ranking, and fill. The default settings require a candidate column with at least two entries, no quadratic participation, and a pivot of largest magnitude in its row. The implication proof replaces the explicit-freedom requirement, independently of `equalities.require_free_variable`.

Equality rows containing a column with at most two entries are considered before other equality rows. Candidates of degree two bypass the configured maximum row and column lengths, so a long equality can still be useful when eliminating a sparse column. Higher-degree candidates obey those size limits. The default settings consider at most one ranked candidate per visit to an equality. Affected rows can be revisited during the same pass.

A successful implication check does not guarantee substitution: the shared substitution routine must also accept the numerical data, fill, and Hessian-growth limits. Rejected candidates leave the problem unchanged and do not establish infeasibility. Affected equalities can be reconsidered after successful substitutions.

## Configuration and scheduling

`rules.implied_free_equalities` defaults to `true` and remains enabled in the aggressive preset. The rule participates in the fast and medium exploration phases and in selected cleanup passes described in [Rule Scheduling](../scheduling.md).

Each invocation uses `equalities.work_limit`, with a default allowance of sixteen times the current constraint nonzero count. Work accounting and the time deadline limit scanning, implication calculations, and substitution attempts. Work and fill limits are separate controls.

## Recovery

Postsolve uses the saved defining equality to reconstruct the eliminated variable from the surviving variables. The shared equality-substitution recovery also maps row and signed variable-bound multipliers back to the original problem. No extra constraint is needed for the eliminated bounds because the defining equality and surviving bounds already guarantee them.
