# Bound propagation

Bound propagation uses linear constraints and current variable bounds to derive tighter variable bounds. It can also detect infeasibility, remove redundant row sides, or expose a fixed variable for another reduction.

Unlike singleton-row elimination, this rule works on rows containing several variables. It does not change the objective.

## Example

Suppose a row is \\(x + y \le  10\\), with \\(y \ge  3\\) and \\(0 \le  x \le  20\\).

Even at its smallest allowed value, \\(y\\) contributes at least \\(3\\), so \\(x\\) cannot exceed \\(7\\). The updated bound is \\(0 \le  x \le  7\\).

Negative coefficients reverse the direction. For \\(-2x + y \le  10\\) and \\(y \ge  4\\), the row implies \\(-2x \le  6\\), hence \\(x \ge  -3\\).

The original row generally remains necessary: knowing \\(x \le  7\\) alone does not enforce \\(x + y \le  10\\) when \\(y\\) can exceed \\(3\\).

## Mathematical formulation

Use the [shared problem notation](../problem_formulation.md).

For each nonzero coefficient in row \\(i\\), define the smallest and largest possible contributions of its variable:

\\[
m_{ij}=\min(a_{ij}L_j,a_{ij}U_j),
\qquad M_{ij}=\max(a_{ij}L_j,a_{ij}U_j).
\\]

The row's minimum and maximum activities over the variable-bound box are

\\[
\alpha_i=\sum_j m_{ij},
\qquad \beta_i=\sum_j M_{ij}.
\\]

These are bounds on the row expression, not a guarantee that all values in the interval satisfy the other constraints. Sums here and below range over nonzero coefficients.

### Infeasibility and redundant sides

In exact arithmetic, \\(\alpha_i>u_i\\) or \\(\beta_i<\ell_i\\) proves infeasibility: even the most favorable variable values cannot satisfy that row side.

If \\(\alpha_i\ge\ell_i\\), the lower row side is already guaranteed by the variable bounds. If \\(\beta_i\le u_i\\), the upper row side is guaranteed. Either redundant side can be replaced by an infinite bound; if both sides are redundant, the row can be deleted.

### Isolate one variable

For a variable \\(x_k\\) with \\(a_{ik}\ne0\\), compute the activity of the remaining terms:

\\[
\alpha_{i,-k}=\sum_{j\ne k}m_{ij},
\qquad \beta_{i,-k}=\sum_{j\ne k}M_{ij}.
\\]

For \\(a_{ik}>0\\), the implied variable bounds are

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

Each implication requires a finite row side and a finite corresponding residual activity. An absent or unbounded side supplies no finite candidate through this formula.

Intersect each accepted candidate with the current variable interval:

\\[
L_k^{\prime}=\max(L_k,\widehat L_k),
\qquad U_k^{\prime}=\min(U_k,\widehat U_k).
\\]

A conflicting candidate can establish infeasibility. Equal finite updated bounds expose a variable for [fixed-variable elimination](./fixed_variables.md).

### Equivalence

In exact arithmetic, the derived bounds hold for every original feasible point, and removed row sides are already implied by retained variable bounds. Consequently, these updates preserve feasible variable vectors, objective values, and optimal solutions. Tightening a variable alone does not justify deleting its source row.

## Numerical safeguards

Activity contradictions are confirmed by direct recomputation before a certificate is returned. Contradictions use the configured feasibility margin; a smaller conflict rejects the candidate rather than storing an inconsistent interval. Redundant-side removal requires direct containment without a tolerance.

Unbounded residuals supply no finite candidate. Cancellation-sensitive residual calculations can be recomputed from the other terms. Nonfinite, excessively large, or insignificant updates can be skipped; these limits are described in [Settings](../settings.md#bound-propagation).

Tighter bounds can trigger more propagation or structural reductions, but work and time limits mean presolve need not reach a mathematical fixed point. See [Rule scheduling](../scheduling.md).

## Recovery

Primal variable values are unchanged. During reverse recovery, multipliers associated with tightened bounds can move back to their source rows, with corresponding adjustments to bound multipliers to preserve stationarity. The source equations remain available even if later reductions remove the rows.

See [Settings](../settings.md) for rule switches and tuning, and [Rule scheduling](../scheduling.md) for ordering and stopping conditions.
