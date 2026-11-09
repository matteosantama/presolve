# Fixed-variable elimination

Fixed-variable elimination substitutes a variable's known value into the objective and constraints, then removes the variable. A variable qualifies when its lower and upper bounds are exactly equal and finite; nearly equal bounds do not qualify.

## Quadratic example

Consider

\\[
f(x,z)=x^2+3xz+2z^2+5z+10,
\\]

with \\(z=7\\). Substitution gives

\\[
f(x,7)=x^2+21x+143.
\\]

The cross term becomes \\(21x\\), while the terms involving only \\(z\\), together with the original constant, contribute \\(143\\). Likewise, the constraint becomes

\\[
x+z\le10 \quad\Longrightarrow\quad x\le3.
\\]

## Mathematical formulation

Using the [shared notation](../problem_formulation.md), suppose \\(x_k=v\\) and let \\(R\\) contain the remaining variable indices.

### Objective update

Substituting \\(x_k = v\\) gives the reduced objective

\\[
f_{\mathrm{red}}(x_R) = \tfrac12 x_R^T P_{RR} x_R + (c_R + P_{Rk}v)^T x_R + c_0 + c_kv + \tfrac12 P_{kk}v^2.
\\]

Thus, the reduced objective data are

\\[
P_{\mathrm{red}} = P_{RR},
\qquad
c_{\mathrm{red}} = c_R + P_{Rk}v,
\qquad
c_{0,\mathrm{red}} = c_0 + c_kv + \tfrac12 P_{kk}v^2.
\\]

Quadratic coefficients between remaining variables are unchanged. Cross terms involving the fixed variable contribute to the remaining linear coefficients, while terms involving only the fixed variable contribute to the constant.

### Linear constraint update

For the linear constraints, substitution gives

\\[
\ell - A_{:k}v
\le A_{:R}x_R
\le u - A_{:k}v.
\\]

Infinite row bounds remain unchanged. Equality constraints are the special case \\(\ell = u = b\\), giving

\\[
A_{:R}x_R = b - A_{:k}v.
\\]

The bounds on the remaining variables are unchanged by this substitution. For a cone row, subtract the fixed contribution from its right-hand side while retaining the cone block.

### Primal recovery and equivalence

Postsolve reconstructs the variable vector by inserting \\(x_k = v\\) at its saved index and retaining the recovered values of \\(x_R\\).

In exact arithmetic, this gives a one-to-one correspondence between feasible points of the reduced and original problems, with

\\[
f_{\mathrm{red}}(x_R) = f(x).
\\]

Consequently, the reduction preserves objective values and optimal solutions.

## Recovery

Postsolve restores \\(x_k=v\\). For linear constraints, its signed bound multiplier follows from stationarity:

\\[
z\_k=(Px+c)\_k-A\_{:k}^Ty.
\\]

Because the original variable is fixed, this combined multiplier can have either sign. Cone-row contributions are included when present, using the [general stationarity convention](../problem_formulation.md#recovery-and-optimality).

See [Settings](../settings.md) for rule switches and tuning, and [Rule scheduling](../scheduling.md) for ordering and stopping conditions.
