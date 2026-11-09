# Empty-row elimination

An empty row is a linear constraint with no nonzero coefficients. Its left-hand side is always zero, so presolve either removes the row when zero satisfies its bounds or reports infeasibility when it does not.

This rule removes constraints, not variables. It does not change the objective or any variable bounds.

## Examples

The row \\(-2 \le  0 \le  7\\) is satisfied for every variable vector and can be removed.

The row \\(2 \le  0 \le  7\\) is impossible because its lower bound exceeds zero. Likewise, \\(-7 \le  0 \le  -2\\) is impossible because its upper bound is below zero. Either row proves that the entire problem is infeasible.

For an empty equality, \\(0 = 0\\) is redundant, while \\(0 = 3\\) proves infeasibility in exact arithmetic.

## Mathematical formulation

Use the [shared problem notation](../problem_formulation.md). If row \\(i\\) has no nonzero coefficients, then \\(A_{i:}=0\\) and its constraint reduces to

\\[
\ell_i \le 0 \le u_i.
\\]

In exact arithmetic, the row is redundant if \\(\ell_i\le0\\) and \\(u_i\ge0\\). Otherwise, no variable vector can satisfy it and the problem is infeasible. Infinite bounds representing an absent lower or upper restriction do not prevent removal when the other side is satisfied.

### Reduced problem and equivalence

Let \\(R\\) contain all row indices except \\(i\\). Removing a satisfied empty row leaves

\\[
\ell_R \le A_{R:}x \le u_R.
\\]

All variables, variable bounds, and objective data remain unchanged. In exact arithmetic, the original and reduced problems have the same feasible variable vectors and objective values, and therefore the same optimal solutions.

An empty row can be present initially or arise during another reduction. For example, fixing \\(x = 3\\) in \\(2x \le  6\\) leaves \\(0 \le  0\\), which this rule can remove. Fixing \\(x = 3\\) in \\(2x \ge  7\\) instead leaves \\(0 \ge  1\\), exposing infeasibility.

## Numerical tolerance

Infeasibility is tested with the feasibility margin described in [Settings](../settings.md#numerical-tolerances). A small violation within tolerance can be treated as satisfied. For example, \\(0 \le  -7e-18\\) can be removed when the configured tolerance covers that residual.

This differs from exact feasible-set equivalence: a recovered point can retain that small violation in the original problem.

## Recovery

No primal variable needs reconstruction. The removed row's multiplier is restored as zero; it contributes no coefficients to stationarity and zero satisfies complementarity for a satisfied row.

See [Settings](../settings.md) for rule switches and tuning, and [Rule scheduling](../scheduling.md) for ordering and stopping conditions.
