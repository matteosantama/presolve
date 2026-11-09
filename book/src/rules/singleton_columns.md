# Singleton-column elimination

A singleton column is a variable with exactly one nonzero coefficient in the constraint matrix. Singleton-column elimination uses that constraint to express the variable in terms of the remaining variables, substitutes the expression into the objective, and removes the variable when the required checks pass.

The variable can still appear in the objective. Its bounds must also be preserved: they can prevent substitution or become a constraint on the remaining variables.

## Equality example

Suppose \\(x\\) appears only in the equality \\(2x + 3y = 10\\) and has no variable bounds. Solve for it:

\\[
x = 5 - \tfrac32 y.
\\]

For the objective \\(x + 2y\\), substitution gives \\(5 + \tfrac12y\\). The reduced problem removes \\(x\\) and this equality. After solving, postsolve restores \\(x\\) using the saved expression.

## Mathematical formulation

Let column \\(k\\) have its only nonzero coefficient in linear row \\(i\\). Write that row as

\\[
\ell_i \le a_{ik}x_k + r(x_R) \le u_i,
\qquad
r(x_R) = \sum_{j\ne k} a_{ij}x_j,
\\]

where \\(R\\) denotes the remaining variables and \\(a_{ik}\ne0\\).

For an equality with \\(\ell_i=u_i=b\\), the substitution is

\\[
x_k = \frac{b-r(x_R)}{a_{ik}}.
\\]

The objective is evaluated with this expression in place of \\(x_k\\). For a quadratic objective this can change both its linear and quadratic coefficients. No other constraint row contains \\(x_k\\).

### Preserve necessary variable bounds

An existing bound on \\(x_k\\) can be omitted if the substitution and the other variables' bounds already guarantee it. Otherwise it must survive as a constraint on \\(r(x_R)\\).

For a positive pivot and retained bounds \\(L_k\le x_k\le U_k\\), this constraint is

\\[
b-a_{ik}U_k \le r(x_R) \le b-a_{ik}L_k.
\\]

For a negative pivot the lower and upper variable bounds exchange roles. If neither bound needs to survive, the row can be removed entirely.

For example, substituting \\(x = 5 - y\\) with a necessary bound \\(x \le  4\\) leaves \\(y \ge  1\\).

In exact arithmetic, equality substitution with all necessary bounds retained gives a one-to-one correspondence between original and reduced feasible points and preserves their objective values.

## Inequality case

For an inequality, choosing a row boundary requires an optimality argument. The project applies this case only when the singleton variable has no quadratic objective terms, including cross terms. Let its linear objective coefficient be \\(c_k\\).

For a minimization problem, the preferred row side is:

- The lower side if \\(c\_k/a\_{ik}>0\\).
- The upper side if \\(c\_k/a\_{ik}<0\\).
- Either finite side if \\(c\_k = 0\\); the implementation prefers the lower side when finite.

The chosen side must be finite and reachable without violating the variable's bounds. Substitution then uses that side as an equality, retaining any necessary bound on the other side of the variable's interval.

This restriction can discard feasible points with worse or equal objectives. Unlike equality substitution, it need not preserve every feasible point; it preserves the optimal value and permits recovery of an original optimum when one exists.

### Reachable boundary example

Consider

\\[
\min_{x,y} x+2y
\quad\text{subject to}\quad
x+y\ge5,\qquad 0\le y\le3,
\\]

with \\(x\\) unbounded. For every fixed \\(y\\), the smallest feasible \\(x\\) is \\(5 - y\\). Substitution leaves the objective \\(5 + y\\), so the reduced optimum is \\(y = 0\\); postsolve recovers \\(x = 5\\).

If \\(1 \le  x \le  4\\) instead, the lower bound does not block reaching \\(x = 5 - y\\) over the original interval for \\(y\\). The upper bound becomes \\(y \ge  1\\). The reduced problem has \\(1 \le  y \le  3\\), giving \\(y = 1\\) and \\(x = 4\\).

### Boundary blocked by a bound

If \\(x \ge  3\\) and \\(0 \le  y \le  3\\), the boundary expression \\(x = 5 - y\\) gives \\(x = 2\\) at \\(y = 3\\), violating the lower bound. The best feasible choice for a given \\(y\\) is now \\(\max(3,5-y)\\).

The singleton-column rule skips this substitution. This does not prove infeasibility or prevent another presolve rule from simplifying the problem.

## Applicability and recovery

For equalities, unnecessary variable bounds can be dropped while necessary bounds survive as a row on the remaining variables. Numerical stability and limits on added coefficients can prevent an otherwise valid elimination.

For inequalities, the chosen boundary must be reachable without violating the variable's bounds. Quadratic participation prevents this linear-cost boundary argument. A rejected substitution does not establish infeasibility.

Postsolve restores the variable using the saved equation and reconstructs the corresponding row and bound multipliers.

See [Settings](../settings.md) for rule switches and tuning, and [Rule scheduling](../scheduling.md) for ordering and stopping conditions.
