# Parallel-row elimination

Parallel-row elimination combines linear constraints whose coefficient rows are proportional. It rescales their bounds to a common row expression, intersects the intervals, and removes one row when the checks permit.

This rule changes constraint rows and their bounds. It does not change the objective or variable bounds, although a resulting equality can enable another reduction.

## Example

Consider \\(x + 2y \le  5\\) and \\(2x + 4y \le  8\\). Multiplying the second row by \\(\tfrac12\\) gives \\(x + 2y \le  4\\).

The two constraints are therefore equivalent to the tighter row \\(x + 2y \le  4\\); the original upper bound of \\(5\\) is redundant.

A negative scale reverses the row sides. For example, \\(-2x - 4y \ge  -8\\) also becomes \\(x + 2y \le  4\\) after multiplication by \\(-\tfrac12\\).

Parallel rows need not merely duplicate an upper bound. Combining \\(x + 2y \ge  4\\) with \\(2x + 4y \le  8\\) gives the equality \\(x + 2y = 4\\). Replacing the first lower bound by \\(5\\) instead makes the pair infeasible in exact arithmetic.

## Mathematical formulation

Suppose rows \\(i\\) and \\(r\\) satisfy

\\[
A_{i:}=\rho A_{r:},\qquad \rho\ne0,
\\]

with row bounds \\(\ell_i\le A_{i:}x\le u_i\\) and \\(\ell_r\le A_{r:}x\le u_r\\). Row \\(i\\) is retained as the base row.

For \\(\rho>0\\), the second row implies bounds on the base expression

\\[
\widehat\ell=\rho\ell_r,\qquad \widehat u=\rho u_r.
\\]

For \\(\rho<0\\), the sides exchange roles:

\\[
\widehat\ell=\rho u_r,\qquad \widehat u=\rho\ell_r.
\\]

Intersect these with the base row's interval:

\\[
\ell_i^{\prime}=\max(\ell_i,\widehat\ell),
\qquad u_i^{\prime}=\min(u_i,\widehat u).
\\]

If the intersection is nonempty, row \\(r\\) can be removed and the base row receives the intersected bounds. Equal finite endpoints make the retained row an equality. An empty intersection proves infeasibility in exact arithmetic.

For exactly proportional rows, this operation preserves all feasible variable vectors and objective values. Approximate coefficient matches require the numerical qualifications below.

## Numerical safeguards

Matching variable support alone is insufficient: coefficients must agree after scaling. The scale must be finite and nonzero. Each coefficient pair is checked using the parallel tolerance \\(\tau\\):

\\[
|a-\rho b|\le\tau\max(|a|,|\rho b|).
\\]

The comparison also distinguishes exact equality under floating-point scaling from tolerance-only agreement. A finite bound that overflows during rescaling rejects the merge.

An empty intersection reports infeasibility only for an exact coefficient match and a gap exceeding the scaled feasibility margin. Otherwise the merge is rejected without changing the rows. A tolerance-only match cannot turn a previously non-equality retained row into an equality.

Accepted approximate merges keep the retained coefficients and discard the other row. Exact feasible-set equivalence requires truly proportional rows; tolerance-only agreement can leave a residual in the original problem. See [Settings](../settings.md#numerical-tolerances).

## Recovery

Primal values do not change. The removed row's multiplier starts at zero. If the retained row's multiplier uses a side tightened by the other row, that contribution transfers back to its source, scaled by \\(\rho\\). For exactly proportional rows, this preserves stationarity; approximate merges can leave a corresponding residual.

See [Settings](../settings.md) for rule switches and tuning, and [Rule scheduling](../scheduling.md) for ordering and stopping conditions.
