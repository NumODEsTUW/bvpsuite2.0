# BVPSUITE2.0

## Scope of bvpsuite2.0

The package **bvpsuite2.0** has been developed at the Institute for Analysis and Scientific Computing, Vienna University of Technology, and can be used for the numerical solution of implicit boundary value problems (BVPs) in singular and regular ordinary differential equations (ODEs) as well as eigenvalue problems (EVPs), and Index-1 Differential-Algebraic Equations (DAEs). The ODE system can have a general implicit form and be of mixed order [\[1\]](#fn1) subject to multi-point boundary conditions (BCs). Furthermore the interval can be finite or semi-infinite [\[2\]](#fn2). The BVP may also include unknown parameters to be calculated together with the unknown solution, in which case additional boundary conditions are required. Moreover, in the scope of the code are parameter-dependent problems, where a path in the solution-parameter space is to be followed. This path may feature turning points
**bvpsuite2.0** is a product of years of research at the Institute for Analysis and Scientific Computing on the analysis and numerical solution of ODEs with time and space singularities. Through the years, the code evolved from **sbvp**, a collocation code for explicit singular ODEs of the first order, to **bvpsuite1.1**, which can also handle explicit and implicit problems of arbitrary order, to **bvpsuite2.0**, which improves on usability and readability of the code. The code in the package is now structured in a simpler modular form.

* * *
1.  *The highest involved derivative may vary with the solution component and it can also be zero, which means that algebraic constraint which do not involve derivatives are also admitted*. [↩︎](#fnref1)
2.  *\[a, ∞), where a ≥ 0 and Dirichlet-type boundary conditions are posed at infinity*. [↩︎](#fnref2)
* * *


## Solver

As in the previous versions of the code, collocation is used for the numerical solution of the underlying boundary value problems. A collocation solution is a piecewise polynomial function which satisfies the given ODE at a finite number of nodes (collocation points). This approach shows advantageous convergence properties compared to other direct higher order methods (see [link 1](https://pdfs.semanticscholar.org/0fc4/408259d7358323a9722f960b9209eab06d8d.pdf)). However, the superconvergence in the context of collocation applied to singular problems can not be guaranteed in general, see for example [link 2](https://www.jstor.org/stable/2156572?seq=1) and [link 3](https://www.jstor.org/stable/2157516?seq=1). Our estimate for the global error of the collocation solution is a classical error estimate based on mesh halving. To make the computations more efficient, an adaptive mesh selection strategy based on the control of residual and the global error estimate has been implemented. A detailed description of this mesh selection algorithm is given in [link 4](https://link.springer.com/article/10.1007/s11075-010-9374-0).  
A flowchart of the code:

![bvpsuite2.0_flowchart.png](https://github.com/NumODEsTUW/bvpsuite2.0/blob/main/bvpsuite2.0Flowcharts/bvpsuite2.0_flowchart.png)

and of the pathfollowing case:

![bvpsuite2.0pathfollowing_flowchart.png](https://github.com/NumODEsTUW/bvpsuite2.0/blob/main/bvpsuite2.0Flowcharts/bvpsuite2.0pathfollowing_flowchart.png)


## How to get started using bvpsuite2.0

The [manual](https://github.com/NumODEsTUW/bvpsuite2.0/blob/main/manual_bvpsuite2.0.pdf) aims to help users to get started with the computations of their own problems with **bvpsuite2.0**.  
At the beginning of the second section of the manual, **bvpsuite2.0** is used to compute the solution of a simple example. Then the model problem in a general form, for which **bvpsuite2.0** can find approximate solutions, is presented. Afterwards, the input arguments and also the output of the function call are briefly discussed.  
Finally, in the third and final section, we present a few examples of BVPs to highlight some of the features of the seven modules of the code.

* * *
* * *

