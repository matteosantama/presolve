# Singleton-row elimination

Singleton-row elimination replaces a linear constraint containing exactly one nonzero coefficient with bounds on its variable, then removes the constraint once those bounds capture its restriction.

The variable remains in the problem. Its objective coefficients and its coefficients in other rows are unchanged by this rule.

## Example

Suppose a problem contains the row \\(6 \le  2x \le  14\\), and the variable already has bounds \\(0 \le  x \le  5\\).

Dividing the row by \\(2\\) gives \\(3 \le  x \le  7\\). Intersecting these bounds with the existing bounds gives \\(3 \le  x \le  5\\).

The row can now be removed: its restriction is fully represented by the variable bounds.

## Mathematical formulation

Using the [shared problem notation](../problem_formulation.md), suppose row \\(i\\) contains exactly one nonzero coefficient, in column \\(k\\):

\\[
\ell_i \le a_{ik}x_k \le u_i,
\qquad a_{ik} \ne 0.
\\]

### Implied variable bounds

For a positive coefficient, divide both row bounds by \\(a_{ik}\\). For a negative coefficient, division reverses the inequalities, so the bounds exchange roles:

For \\(a_{ik}>0\\),

\\[
\widehat L_k=\ell_i/a_{ik},
\qquad \widehat U_k=u_i/a_{ik}.
\\]

For \\(a_{ik}<0\\),

\\[
\widehat L_k=u_i/a_{ik},
\qquad \widehat U_k=\ell_i/a_{ik}.
\\]

Intersect the implied interval with the existing variable bounds:

\\[
L_k^{\prime} = \max(L_k,\widehat L_k),
\qquad
U_k^{\prime} = \min(U_k,\widehat U_k).
\\]

For example, the row \\(6 \le -2x_k \le 14\\) implies \\(-7 \le x_k \le -3\\).

### Infeasibility, redundancy, and fixed variables

In exact arithmetic, \\(L_k^{\prime} > U_k^{\prime}\\) proves infeasibility. For example, \\(6 \le  2x \le  14\\) conflicts with an existing upper bound \\(x \le  2\\).

If the existing interval already lies inside the implied interval, the row is redundant and can be removed without changing the bounds. Thus \\(3 \le  x \le  7\\) adds nothing when \\(4 \le  x \le  6\\) already holds.

If the updated bounds coincide at a finite value, the variable becomes eligible for [fixed-variable elimination](./fixed_variables.md). Singleton-row elimination itself does not substitute or remove the variable.

### Reduced problem and equivalence

Let \\(R\\) contain all row indices except \\(i\\). The reduced linear constraints are

\\[
\ell_R \le A_{R:}x \le u_R,
\\]

with the updated bounds on \\(x_k\\); all other variable bounds are unchanged. Here, \\(A_{R:}\\) selects the remaining rows and all columns.

In exact arithmetic, the reduced and original problems have the same feasible variable vectors. The objective is unchanged, so objective values and optimal solutions are preserved. No primal variable needs to be reconstructed for this reduction alone.

## Numerical tolerance and recovery

A conflicting implied bound establishes infeasibility only when the gap exceeds the configured feasibility margin. Nonfinite candidates are not applied. The row is removed only when the stored variable bounds actually contain its restriction, using direct comparisons rather than tolerance alone.

Primal values do not change. During dual recovery, a bound multiplier attributable to a transferred row restriction can move back to that row. For coefficient \\(a_{ik}\\), transferring \\(z_k\\) contributes \\(z_k/a_{ik}\\) to the row multiplier, preserving stationarity in exact arithmetic.

See [Settings](../settings.md) for rule switches and tuning, and [Rule scheduling](../scheduling.md) for ordering and stopping conditions.
