# Empty-column elimination

An empty column is a variable with no nonzero coefficients in the constraint matrix. Its objective terms and variable bounds determine whether presolve can choose an optimal value, remove it by quadratic elimination, or identify an improving recession direction.

An empty constraint column does not necessarily mean an independent variable: quadratic cross terms can still couple it to other variables in the objective.

## Linear examples

Suppose \\(x\\) appears in no constraints and has bounds \\(1 \le  x \le  5\\).

- With objective \\(3x + 2y\\), the best choice is \\(x = 1\\). Removing \\(x\\) leaves the objective \\(3 + 2y\\).
- With objective \\(-3x + 2y\\), the best choice is \\(x = 5\\), leaving \\(-15 + 2y\\).
- With no objective term involving \\(x\\), any value in its interval works. The implementation chooses zero clipped to the interval, giving \\(x = 1\\) here.

If the objective coefficient is negative and there is no upper bound, increasing \\(x\\) decreases the objective without limit. This establishes an improving direction; the full problem is unbounded below if the remaining constraints are feasible.

## Mathematical formulation

Using the [shared problem notation](../problem_formulation.md), suppose column \\(k\\) is empty, so \\(A_{:k}=0\\). First assume there are no quadratic cross terms involving \\(x_k\\): \\(P_{jk}=0\\) for every \\(j\ne k\\).

Let \\(R\\) denote the remaining variable indices and \\(p=P_{kk}\\). The variable's independent objective contribution is

\\[
q_k(t)=\tfrac12 p t^2+c_k t,
\qquad L_k\le t\le U_k.
\\]

### Choosing the value

For \\(p>0\\), the unconstrained minimizer is \\(-c_k/p\\). Project it onto the variable's interval:

\\[
v=\min\bigl(U_k,\max(L_k,-c_k/p)\bigr).
\\]

For \\(p=0\\), choose the lower bound when \\(c_k>0\\), or the upper bound when \\(c_k<0\\). When both \\(p\\) and \\(c_k\\) are zero, the implementation chooses

\\[
v=\min\bigl(U_k,\max(L_k,0)\bigr).
\\]

A finite chosen value can be substituted and the variable removed. A negative diagonal coefficient is skipped by this rule.

### Reduced problem and optimality

With a finite minimizing value \\(v\\), the reduced objective is

\\[
\tfrac12 x_R^T P_{RR}x_R+c_R^Tx_R+c_0+q_k(v).
\\]

The remaining constraints use \\(A_{:R}\\) with unchanged row bounds. The remaining variable bounds are also unchanged.

Every feasible reduced point can be extended with \\(x_k=v\\), and any original feasible point can have its \\(x_k\\) replaced by \\(v\\) without increasing the objective. Thus the reduction preserves the optimal value and permits recovery of an original optimum when one exists. It does not preserve every original feasible point or its objective value.

### Quadratic example

For the independent objective contribution \\(x^2 - 6x\\) and bounds \\(0 \le  x \le  2\\), the quadratic convention gives \\(p = 2\\) and \\(c = -6\\). The unconstrained minimizer is \\(3\\), so the chosen value is \\(2\\). Substitution contributes \\(-8\\) to the objective constant.

### Improving recession direction

For \\(p=0\\) and \\(c_k\ne0\\), a nonfinite improving endpoint gives the direction

\\[
d_k=-\operatorname{sign}(c_k),
\qquad d_j=0\quad(j\ne k).
\\]

Because the constraint column is empty, \\(Ad=0\\). There are no quadratic terms involving the variable, and \\(c^Td=-|c_k|<0\\). The absent bound permits movement indefinitely in this direction. This direction does not itself establish that the remaining constraints are feasible.

## Coupled quadratic variables

An unbounded variable with an empty constraint column but quadratic cross terms may instead be removed using its stationarity condition, provided \\(P_{kk}>0\\):

\\[
x_k=-\frac{c_k+\sum_{j\ne k}P_{kj}x_j}{P_{kk}}.
\\]

This is an affine substitution, not a constant fixing. It requires finite transformed data and respects fill limits without permitting net Hessian growth. Bounded coupled variables are retained by this rule because their minimizing expression can require clipping. See [Settings](../settings.md) for enabling quadratic elimination.

## Recovery

For a constant fixing, postsolve restores the chosen value and its signed bound multiplier \\(pv+c_k\\), since the eliminated column had no constraint entries. For a coupled quadratic elimination it restores the variable using the affine expression above, after recovering the variables on which it depends. Numerical rejection of an elimination does not establish infeasibility.

See [Settings](../settings.md) for rule switches and tuning, and [Rule scheduling](../scheduling.md) for ordering and stopping conditions.
