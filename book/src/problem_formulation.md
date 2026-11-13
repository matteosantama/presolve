# Problem formulation

The presolver works with continuous variables and convex quadratic objectives. The linear part of the problem is written as

\\[
\min_x f(x)=\tfrac12 x^T P x+c^T x+c_0,
\qquad \ell\le Ax\le u,
\qquad L\le x\le U.
\\]

Here, the symmetric positive semidefinite matrix \\(P\\) defines the quadratic terms, \\(c\\) the linear objective, and \\(c_0\\) the constant. A linear objective is the special case \\(P=0\\). Each row of \\(A\\) is a linear constraint; each column corresponds to a variable.

The vectors \\(\ell,u\\) contain row bounds, while \\(L,U\\) contain variable bounds. A lower bound of negative infinity or an upper bound of positive infinity represents an absent restriction. An equality has \\(\ell_i=u_i=b_i\\); \\(Ax=b\\) denotes a collection of such equalities. Equal finite variable bounds fix a variable.

## Indices and submatrices

Subscripts select entries or submatrices. In particular:

- \\(A_{i:}\\) is row \\(i\\), and \\(A_{:k}\\) is column \\(k\\).
- \\(R\\) denotes a set of retained indices, specified by each reduction.
- \\(x_R\\) selects retained variables; \\(P_{RR}\\) selects their quadratic submatrix.
- \\(A_{:R}\\) retains columns and \\(A_{R:}\\) retains rows.

A nonzero count refers to stored nonzero coefficients. A singleton row or column contains exactly one nonzero coefficient; a doubleton contains exactly two.

## Cone constraints

The same problem can also contain cone blocks. Writing \\(G\\) for the matrix formed by the cone rows, these constraints have the form

\\[
Gx+s=h,\qquad s\in K,
\\]

where \\(K\\) is a product of convex cones. In the native API, linear and cone rows share `Problem::a`; the tags in `Problem::rows` distinguish their meanings. The rule pages use \\(A\\) for the linear rows unless otherwise stated. A rule requiring linear rows does not automatically apply to cone rows.

## Input representation

`Problem::p` contains the upper triangle of \\(P\\); `None` represents a linear objective. `Problem::c` and `Problem::c0` hold the remaining objective data. `Problem::variable_bounds` holds variable intervals; an empty vector means every variable is free.

Input validity is the caller's responsibility: dimensions, sorted unique sparse row indices, finite coefficients, convexity, and cone geometry must be consistent. Infinite endpoints represent missing bounds, not infinite matrix coefficients.

## Recovery and optimality

A reduced problem is solved in its new coordinates, then postsolve maps the result to the original coordinates. Some rules preserve every feasible point through an invertible substitution; others retain a subset sufficient to recover an optimum or aggregate several original points into one reduced point. Each rule states its mathematical guarantee.

For linear rows, signed multipliers satisfy

\\[
Px+c=A^Ty+z.
\\]

Here \\(y\\) contains row multipliers and \\(z\\) variable-bound multipliers. Both are positive on lower sides and negative on upper sides. An equality or fixed variable can have a multiplier of either sign. Cone multipliers \\(w\\) use the dual-cone convention, giving

\\[
Px+c-A^Ty-z+G^Tw=0.
\\]

Equivalence arguments in the rule pages use exact arithmetic. Floating-point safeguards and feasibility tolerances can affect which reductions are accepted. An improving recession direction alone does not prove feasibility; an unboundedness result requires a feasible point as well.
