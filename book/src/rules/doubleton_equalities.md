# Doubleton equality elimination

A doubleton equality is a linear equality with exactly two nonzero coefficients. Presolve can solve it for one variable, substitute that expression into the objective and other constraints, and remove the variable when numerical and sparsity checks permit.

Unlike singleton-column elimination, the eliminated variable can appear in several other constraint rows. Its bounds must also survive the substitution.

## Example

Suppose the problem contains \\(2x + 3y = 12\\), another row \\(x + z \le  8\\), and the objective \\(4x + y + z\\).

Solving the equality for \\(x\\) gives

\\[
x=6-\tfrac32 y.
\\]

Substitution changes the other row to \\(-\tfrac32y + z \le  2\\) and the objective to \\(24 - 5y + z\\).

If \\(x\\) is unbounded, the defining equality can be removed. If \\(0 \le  x \le  6\\), those bounds instead become \\(0 \le  y \le  4\\). These restrictions can be retained as a singleton row; singleton-row elimination can convert it into variable bounds.

## Mathematical formulation

Consider the equality

\\[
a x_j+b x_k=d,
\qquad a\ne0,\quad b\ne0.
\\]

Eliminating \\(x_j\\) expresses it in terms of the surviving variable:

\\[
x_j=\delta+s x_k,
\qquad \delta=d/a,\quad s=-b/a.
\\]

### Other constraint rows

For another row \\(r\\), substitution gives

\\[
\ell_r-A_{rj}\delta
\le \sum_{t\ne j,k}A_{rt}x_t+(A_{rk}+A_{rj}s)x_k
\le u_r-A_{rj}\delta.
\\]

Thus a row containing the eliminated variable changes its bounds and its coefficient of the survivor. A survivor coefficient can be created, combined with an existing coefficient, or canceled entirely. Infinite row bounds remain infinite.

### Preserve variable bounds

Bounds \\(L_j\le x_j\le U_j\\) become

\\[
L_j\le\delta+s x_k\le U_j.
\\]

For a positive pivot \\(a\\), these bounds are equivalent to the singleton row

\\[
d-aU_j\le b x_k\le d-aL_j.
\\]

For a negative pivot, the lower and upper variable bounds exchange roles. These restrictions must be combined with the survivor's existing bounds; they can expose infeasibility or a fixed variable. If the eliminated variable has no bounds, the defining row can be deleted instead.

### Objective update and equivalence

Let \\(z\\) contain all variables except \\(x_j\\). Write the reconstruction as \\(x=Tz+h\\), where the surviving coordinates are copied unchanged and coordinate \\(j\\) is \\(\delta+s x_k\\). Here \\(h_j=\delta\\) and the other entries of \\(h\\) are zero.

For the [shared quadratic objective](../problem_formulation.md), substitution gives

\\[
P_{\mathrm{red}}=T^TPT,
\qquad c_{\mathrm{red}}=T^T(c+Ph),
\qquad c_{0,\mathrm{red}}=c_0+c^Th+\tfrac12 h^TPh.
\\]

Cross terms can therefore change the survivor's quadratic coefficients as well as linear coefficients and the constant. For a linear objective, the eliminated term contributes \\(c_j\delta\\) to the constant and \\(c_js\\) to the survivor's coefficient.

In exact arithmetic, substitution with all necessary variable bounds retained gives a one-to-one correspondence between feasible points of the reduced and original problems, with equal objective values. Optimal solutions are preserved.

## Applicability

Elimination requires a usable pivot and finite transformed data. The selected coefficient ratio must have magnitude between \\(10^{-7}\\) and \\(10^7\\). For quadratic variables with unequal coefficient magnitudes, the larger pivot is preferred to limit amplification. Otherwise singleton columns, approximately integral ratios, and shorter columns guide the choice.

Necessary bounds on the eliminated variable survive as a singleton constraint. When enabled, singleton-row elimination converts it into bounds on the survivor. Unbounded variables need no retained bound row.

Suitable trees of linear equalities can be reduced together when coefficients have equal magnitudes and all eliminated variables are free. This changes how much work is needed, not the underlying substitution argument. Numerical checks and limits on added coefficients can reject individual or grouped reductions.

## Recovery

Postsolve uses the saved defining equality to reconstruct the eliminated variable after recovering the surviving variables. It also maps row and signed variable-bound multipliers back to the original problem. With necessary bounds retained, the recovered point has the same objective value in exact arithmetic.

See [Settings](../settings.md) for rule switches and tuning, and [Rule scheduling](../scheduling.md) for ordering and stopping conditions.
