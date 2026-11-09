# Parallel-column elimination

Parallel-column elimination combines variables whose constraint columns are proportional when their objective terms also depend only on the same combined variable. It removes one variable and replaces the retained variable's interval with the interval of possible combined values.

Proportional constraint columns alone are not sufficient. The implementation requires compatible quadratic terms and, for aggregation, compatible linear objective coefficients. When linear costs differ, a separate fallback can fix a variable or identify an improving recession direction.

## Example

Suppose every constraint uses \\(x\\) and \\(y\\) through \\(x + 2y\\), the objective is \\(3x + 6y + w\\), and their bounds are \\(1 \le  x \le  4\\) and \\(0 \le  y \le  3\\).

Introduce

\\[
z=x+2y,\qquad 1\le z\le10.
\\]

The row \\(x + 2y + w \le  12\\) becomes \\(z + w \le  12\\), and the objective becomes \\(3z + w\\). Two variables have become one.

If the reduced solution has \\(z = 7\\), one valid recovery is \\(x = 1\\), \\(y = 3\\). Recovery need not reproduce a particular original pair: several pairs can represent the same aggregate.

## Mathematical formulation

Use the [shared problem notation](../problem_formulation.md).

For a nonzero scale \\(\rho\\), suppose

\\[
A_{:k}=\rho A_{:j},
\qquad P_{:k}=\rho P_{:j},
\qquad c_k=\rho c_j.
\\]

These relations make the constraints and objective depend on the pair only through

\\[
z=x_j+\rho x_k.
\\]

The retained coordinate \\(j\\) represents \\(z\\). Remove column \\(k\\) from the constraint matrix, remove row and column \\(k\\) from the quadratic matrix, and drop \\(c_k\\). All remaining coefficients, row bounds, and the objective constant stay unchanged; coefficients associated with coordinate \\(j\\) now multiply the aggregate.

### Aggregate bounds

For \\(\rho>0\\),

\\[
L_z=L_j+\rho L_k,\qquad U_z=U_j+\rho U_k.
\\]

For \\(\rho<0\\),

\\[
L_z=L_j+\rho U_k,\qquad U_z=U_j+\rho L_k.
\\]

All other variable bounds are unchanged. For continuous variables, every value in this interval can be split into an original pair within its bounds.

In exact arithmetic, aggregation maps every original feasible point to a reduced feasible point with the same objective value, and every reduced feasible point admits such an original pair. This preserves the optimal value and permits recovery of optimal solutions, but the correspondence need not be one-to-one.

## Objective-dominance fallback

If the constraint and quadratic relations hold but linear costs differ, consider

\\[
\Delta=c_k-\rho c_j,
\qquad d_k=-\operatorname{sign}(\Delta),
\qquad d_j=-\rho d_k.
\\]

All other direction components are zero. The direction satisfies \\(Ad=0\\) and \\(Pd=0\\), while \\(c^Td=-|\Delta|<0\\). It preserves the constraints and quadratic contribution while decreasing the objective linearly.

The code examines the variable endpoints in these improving directions:

- If neither endpoint is finite, it returns an improving recession certificate. The full problem is unbounded below if its remaining constraints are feasible.
- If only the retained variable's endpoint is infinite, it fixes the other variable at its finite improving endpoint.
- If only the other variable's endpoint is infinite, it fixes the retained variable at its finite improving endpoint.
- If both endpoints are finite, this fallback does not fix either variable.

For example, with identical constraint columns, costs \\(2x + y\\), and both variables nonnegative without upper bounds, transferring from \\(x\\) to \\(y\\) lowers the cost. Since \\(y\\) can increase indefinitely, the fallback can fix \\(x = 0\\) without changing the attainable optimum.

## Numerical safeguards

This rule requires exact floating-point scaled agreement for the constraint and quadratic columns, and compatible linear costs for aggregation. Tolerance-only matches are rejected. The scale and compared coefficients must be finite, and the scale nonzero.

Aggregate bounds reject NaNs and overflow from finite endpoints. A failed aggregation does not establish infeasibility. Exact floating-point comparisons are not symbolic proofs; the equivalence argument assumes exact proportionality and arithmetic.

## Recovery

Postsolve splits the aggregate into a pair satisfying the saved original bound intervals and \\(z=x_j+\rho x_k\\). The split need not be unique. Signed variable-bound multipliers satisfy \\(z_k=\rho z_j\\) under the proportional gradient and constraint relations; here the subscripted quantities are multipliers, not the aggregate variable.

A dominance-based fixing instead uses ordinary [fixed-variable recovery](./fixed_variables.md). A broader [dominated-column rule](./dominated_columns.md) can consider nonproportional columns when enabled.

See [Settings](../settings.md) for rule switches and tuning, and [Rule scheduling](../scheduling.md) for ordering and stopping conditions.
