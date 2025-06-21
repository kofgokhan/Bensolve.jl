Following is taken from [Bensolve](http://bensolve.org/) documentation:

# Bensolve.jl

[![Build Status](https://github.com/kofgokhan/Bensolve.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/kofgokhan/Bensolve.jl/actions/workflows/CI.yml?query=branch%3Amain)

Bensolve is a solver for vector linear programs (VLP), in particular, for the subclass of
multiple objective linear programs (MOLP). The present version is based on Benson’s
algorithm and its extensions, see e.g. [3, 4, 1, 2, 10, 8, 5, 6, 9] and the references therein.
For the theoretical background of this program, the reader is referred to [8, 6].
The present version utilizes the GNU Linear Programming Kit (GLPK). 

Bensolve assumes the following formulation of a MOLP:
```math
\begin{align}

\begin{split}
\text{minimize (or maximize)} \quad & 

\begin{pmatrix}
P_{11} x_1 + P_{12} x_2 + \dots + P_{1n} x_n \\
P_{21} x_1 + P_{22} x_2 + \dots + P_{2n} x_n \\
\dots \\
P_{q1} x_1 + P_{q2} x_2 + \dots + P_{qn} x_n \\
\end{pmatrix}

\end{split} \\

\begin{split}

\text{subject to} \quad
a_1 & \le B_{11} x_1 + B_{12} x_2 + \dots + B_{1n} x_n \le b_1 \\
a_2 & \le B_{21} x_1 + B_{12} x_2 + \dots + B_{1n} x_n \le b_2 \\
& \dots \\
a_m & \le B_{m1} x_1 + B_{m2} x_2 + \dots + B_{mn} x_n \le b_m
\end{split} \\

\begin{split}
l_1 & \le x_1 \le s_1 \\
l_2 & \le x_2 \le s_2 \\
& \dots \\
l_n & \le x_n \le s_n
\end{split}

\end{align}
```

# Vector Linear Program
In this more general setting, the minimization (or maximization) in (1) is defined with
respect to a partial ordering $\le_C$ induced by a polyhedral cone $C \subseteq \mathbb{R}^q$, that is

```math
\begin{align}

\begin{pmatrix}
y_1 \\ y_2 \\ \dots \\ y_q
\end{pmatrix}  \le_C 
\begin{pmatrix}
z_1 \\ z_2 \\ \dots \\ z_q
\end{pmatrix} \iff
\begin{pmatrix}
z_1 - y_1 \\ z_2 - y_2 \\ \dots \\ z_q - y_q
\end{pmatrix} \in C.

\end{align}
```

Bensolve assumes that $C$ has a non-empty interior and contains no lines. The polyhedral cone $C$ is assumed to be given by one of the following two representations:

The `CONE` representation is given by a matrix $Y$ with $q$ rows and $o$ columns. A vector $(y_1, \dots, y_q)$ belongs to $C$ if and only if there are nonnegative real numbers $v_1, v_2, \dots, v_o \ge 0$ such that

```math
\begin{align}
\begin{split}
y_1 & = Y_{11} v_1 + Y_{12} v_2 + \dots + Y_{1o} v_o \\
y_2 & = Y_{21} v_1 + Y_{22} v_2 + \dots + Y_{2o} v_o \\
& \dots \\
y_q & = Y_{q1} v_1 + Y_{q2} v_2 + \dots + Y_{qo} v_o.
\end{split}
\end{align}
```

The columns of the matrix $Y$ are generating vectors of the polyhedral cone $C$.
The `DUALCONE` representation is given by a matrix $Z$ with $q$ rows and $p$ columns. A vector $(y_1, \dots, y_q)$ belongs to $C$ if and only if the following inequalities are satisfied:

```math
\begin{align}
\begin{split}
Z_{11} y_1 + Z_{21} y_2 + \dots + Z_{q1} y_q & \ge 0 \\ 
Z_{12} y_1 + Z_{22} y_2 + \dots + Z_{q2} y_q & \ge 0 \\ 
& \dots \\
Z_{1p} y_1 + Z_{2p} y_2 + \dots + Z_{qp} y_q & \ge 0. 
\end{split}
\end{align}
```

The columns of the matrix $Z$ are generating vectors of the dual cone of the polyhedral cone $C$.

# Using the package

You can install the package from the package manager.

```
] add https://github.com/kofgokhan/Bensolve.jl
```

You can then test the package with one of the example problems. Let us take a look at one of the examples:

```math
\begin{align*}
\text{minimize} \quad & 
\begin{pmatrix}
x_1 - x_2 \\
x_1 + x_2
\end{pmatrix} \\
\text{subject to} \quad 
& 2 x_1 + x_2 \ge 6 \\
& x_1 + 2 x_2 \ge 6 \\
& x_1 \ge 0 \\
& x_2 \ge 0
\end{align*}
```

```
using Bensolve

status, upper_img, lower_img, solve_time = solve("ex01.vlp")
```

Alternatively, you can pass in the parameters of the problem yourself.

```
P = [
    1 -1
    1  1
]

B = [
    2 1
    1 2
]

a = [6, 6]
b = [Inf, Inf]

l = [0, 0]
s = [Inf, Inf]

molp_solve(P, B, a, b, l, s)
```

Using JuMP:

```
using JuMP, Bensolve

model = Model(Bensolve.Optimizer)
@variable(model, l[i] <= x[i = 1:2] <= s[i])
@constraint(model, a .<= B * x .<= b)
@objective(model, Min, P * x)
optimize!(model)

status = termination_status(model)
N = result_count(model) # number of vertices and extreme directions of the upper image
objective_value(model, result = 1) # first vertex / extreme direction
```

# Citation policy
If you are using this version of Bensolve for scientific papers, please cite it as:

[A] Löhne, A, Weißing, B.: Bensolve - VLP solver, version 2.1.x, www.bensolve.org

```
@article{,
   author = {Andreas Löhne and Benjamin Weißing},
   doi = {10.1016/j.ejor.2016.02.039},
   issn = {03772217},
   issue = {3},
   journal = {European Journal of Operational Research},
   month = {8},
   pages = {807-813},
   title = {The vector linear program solver Bensolve – notes on theoretical background},
   volume = {260},
   year = {2017},
}
```
[B] Löhne, A, Weißing, B.: The vector linear program solver Bensolve – notes on theoretical background, European J. Oper. Res., 260(3):807–813, 2017

```
@misc{,
   author = {Andreas Löhne and Benjamin Weißing},
   publisher = {www.bensolve.org},
   title = {Bensolve - VLP solver},
}
```


# References

[1] H. Benson. Further analysis of an outcome set-based algorithm for multiple-objective
linear programming. Journal of Optimization Theory and Applications, 97(1):1–10,
1998.

[2] H. Benson. An outer approximation algorithm for generating all efficient extreme
points in the outcome set of a multiple objective linear programming problem. Journal
of Global Optimization, 13:1–24, 1998.

[3] J. P. Dauer. Analysis of the objective space in multiple objective linear programming.
J. Math. Anal. Appl., 126(2):579–593, 1987.

[4] J. P. Dauer and Y.-H. Liu. Solving multiple objective linear programs in objective
space. European J. Oper. Res., 46(3):350–357, 1990.

[5] M. Ehrgott, A. Löhne, and L. Shao. A dual variant of Benson’s outer approximation
algorithm. Journal of Global Optimization, 52(4):757–778, 2012.

[6] A. H. Hamel, A. Löhne, and B. Rudloff. A Benson type algorithm for linear vector
optimization and applications. Journal of Global Optimization, 59(4):811–836, 2014.

[7] F. Heyde and A. Löhne. Geometric duality in multiple objective linear programming.
SIAM J. Optim., 19(2):836–845, 2008.

[8] A. Löhne. Vector Optimization with Infimum and Supremum. Vector Optimization.
Springer, Berlin, 2011.

[9] A. Löhne and B. Weißing. The vector linear program solver Bensolve – notes on
theoretical background. European J. Oper. Res., 260(3):807–813, 2017.

[10] L. Shao and M. Ehrgott. Approximately solving multiobjective linear programmes in
objective space and an application in radiotherapy treatment planning. Math. Methods
Oper. Res., 68(2):257–276, 2008.
