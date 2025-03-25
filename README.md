# MultiPolyEig
Global Solver for Polynomial Multiparameter Eigenvalue Problems

Implements the hidden variable Dixon resultant method from
[1] E. Graf and A. Townsend, "A Hidden Variable Resultant Method for the
Polynomial Multiparameter Eigenvalue Problem," 2025.

For d=1,2,3, multipolyeig(F,options) returns the common eigenvalues of 
F = cell(F1,...,Fd),where Fi are d-parameter polynomial matrices, given 
as d+2 dimensional tensors where the last d dimensions correspond to 
variables and the first two to matrix dimensions. The output is given as a 
matrix with rows that give eigentuples (x_1,...,x_d) of F.

Options in options are res, the final tolerance for the residual check,
and evblock, which determines the number of entries used to compute 
the first d-1 coordinates of solutions.

Test files implement the examples given in section 6 of [1].

The example in testLeakyWaves requires RandomJointEig, 
https://github.com/borplestenjak/RandomJointEig,
Bor Plestenjak, Hauke Gravenkamp, Daniel Kiefer 2024.


