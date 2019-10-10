BPALM package
=============


This matlab package provides generic solvers for multi-block Bregman proximal alternating linearized minimization (BPALM) and its adaprive variants (A-BPALM1 and A-BPALM2) for solving structured nonsmooth nonconvex problems of the form

  min_x f(x)+sum_{i=1}^N g_i(x_i),

where f is relatively smooth and g_i (i=1,...,n) are proper and lower semicontinuous.

We provide three solvers for solving the above-mentioned problem:
- BPALM: Bregman Proximal Alternating Linearized Minimization;
- A-BPALM1: an adaptive BPALM (see Algorithm 2 in [1]);
- A-BPALM2: an adaptive BPALM (see Remark 3.13 in [1]);

# Using the package

## Installation

All the functions in the package are given in the folder /MatlabCodes.

In order to use the function we recommend to execute the following command

```Matlab
addpath(genpath('.'))
```

if you are not working in the root folder of the package or replacing '.' by the location of the folder on your machine.


## Example: orthogonal nonnegative matrix factorization (ONMF)

We recommend to look at the following files to see how to use the package:
* Demo/demo_ONMF.m: contains an example for comparing a fixed penalty version of BPALM, A-BPALM1, and A-BPALM2 for ONMF with synthetic data.
* Demo/demo_continuation_ONMF.m: contains an example for comparing a continuation version of BPALM, A-BPALM1, and A-BPALM2 for ONMF with synthetic data.

## Solving your own optimization problem

You need to to write:
- a function for providing the function value and the gradient of the objective at point x (add the file to folder MatlabCodes/Test_functions);
- a function for providing the function value and the gradient of the kernel at point x (add the file to folder MatlabCodes/Test_functions);
- a function for providing a solution of the subproblem of BPALM (add the file to folder MatlabCodes/Subproblems);
- a demo file like demo_ONMF.m or demo_continuation_ONMF.m for calling the algorithms.

# References

[1] M. AHOOKHOSH, L.T.K. HIEN, N. GILLIS, AND P. PATRINOS, Multi-block Bregman proximal alternating linearized minimization and its application to orthogonal nonnegative matrix factorization, (2019) https://arxiv.org/pdf/1908.01402.pdf

