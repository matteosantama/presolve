# Short equality elimination

Short equality elimination solves a linear equality containing several variables for one selected variable, substitutes that expression into the objective and other constraints, and removes the variable when numerical and sparsity checks permit.

It generalizes [doubleton equality elimination](./doubleton_equalities.md) to rows with at least three nonzero coefficients. The implementation limits row size, candidate variables, substitution work, and new coefficients to avoid expensive or unstable reductions.

## Example

Suppose a problem contains \\(3x + 2y - z = 12\\), another row \\(x + w \le  7\\), and the objective \\(6x + y + w\\). Assume \\(x\\) has no bounds or quadratic objective terms.

Solving for \\(x\\) gives

\\[
x=4-\tfrac23y+\tfrac13z.
\\]

The other row becomes \\(-\tfrac23y + \tfrac13z + w \le  3\\), and the objective becomes \\(24 - 3y + 2z + w\\). The defining equality and variable \\(x\\) can be removed if the substitution passes the implementation's checks.

This example illustrates why fill matters: replacing \\(x\\) can introduce both \\(y\\) and \\(z\\) into rows that previously contained neither.

## Mathematical formulation

Consider an equality in row \\(i\\), with a nonzero pivot coefficient \\(a_{ik}\\):

\\[
\sum_j a_{ij}x_j=b_i.
\\]

Let \\(R\\) contain all variable indices except \\(k\\), and set

\\[
\delta=b_i/a_{ik},
\qquad s_j=-a_{ij}/a_{ik}\quad(j\in R).
\\]

Coefficients absent from the row have \\(s_j=0\\). The substitution is

\\[
x_k=\delta+s^Tx_R.
\\]

For another linear row \\(r\\), substitution gives

\\[
\ell_r-A_{rk}\delta
\le (A_{rR}+A_{rk}s^T)x_R
\le u_r-A_{rk}\delta.
\\]

Rows containing the eliminated variable can gain, combine, or lose coefficients. Bounds on the eliminated variable become

\\[
L_k\le\delta+s^Tx_R\le U_k.
\\]

Any sides already implied by the equality and the other variables' bounds can be omitted. Necessary sides must remain as a constraint on the surviving variables. If no sides remain, the defining equality can be deleted entirely.

### Objective update and equivalence

Write the reconstruction as \\(x=Tx_R+h\\): surviving coordinates are copied unchanged, coordinate \\(k\\) uses the slopes \\(s\\), and \\(h_k=\delta\\) with all other entries zero.

For the [shared quadratic objective](../problem_formulation.md), the reduced data are

\\[
P_{\mathrm{red}}=T^TPT,
\qquad c_{\mathrm{red}}=T^T(c+Ph),
\qquad c_{0,\mathrm{red}}=c_0+c^Th+\tfrac12h^TPh.
\\]

In exact arithmetic, retaining all necessary bound restrictions gives a one-to-one correspondence between feasible original and reduced points, with equal objective values. The reduction preserves optimal solutions. The formulas also apply to quadratic candidates when the settings permit them; the defaults restrict candidates to variables without quadratic participation.

## Applicability and tuning

The default configuration considers equalities with 3–8 nonzero coefficients, eliminating only unbounded variables without quadratic participation whose constraint columns have 2–8 entries. It requires a pivot with the largest coefficient magnitude in the row and retains one candidate per row per pass.

The [equality settings](../settings.md#equality-substitution) allow bounded or quadratic candidates, alternative pivots, and larger searches. Necessary variable bounds remain as a constraint on the survivors. A retained row can still contain several variables.

Substitution is accepted only when transformed data pass numerical checks and the work, fill, and Hessian-growth policies allow it. These are independent controls; increasing a work budget does not remove a fill limit. A rejected candidate does not establish infeasibility.

## Recovery

Postsolve uses the saved defining equality to reconstruct the eliminated variable after recovering the surviving variables. It also maps row and signed variable-bound multipliers back to the original problem. With necessary bounds retained, the recovered point has the same objective value in exact arithmetic.

See [Settings](../settings.md) for rule switches and tuning, and [Rule scheduling](../scheduling.md) for ordering and stopping conditions.
