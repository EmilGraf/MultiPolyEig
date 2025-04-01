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

Emil Graf and Alex Townsend, 2025.

# Test Files

Test files implement the examples given in section 6 of [1].

# Aeroelastic Flutter

The file testAEF implements the example in section 6.1 of [1], which is a 
2PEP related to aeroelastic flutter from A. Pons, S. Gutschmidt, 
Multiparameter spectral analysis for aeroelastic instability problems, 
J. Appl. Mech. 85 (6) (2018) 061011. doi:10.1115/ 1.4039671.

The example in testAEF requires files from MultiParEig, 
https://www.mathworks.com/matlabcentral/fileexchange/47844-multipareig, Bor Plestenjak, 2025.
To run this example, first download and install MultiParEig.

# Leaky Waves

The file testLeakyWaves implements the example in section 6.2 of [1],
which is a 3PEP that computes wavenumbers of leaky waves from
H. Gravenkamp, B. Plestenjak, D. A. Kiefer, E. Jarlebring, 
Computation of leaky waves in layered structures coupled to unbounded 
media by exploiting multiparameter eigenvalue problems, 
J. Sound Vib. 596 (2025) 118716. https://doi.org/10.1016/j.jsv.2024.118716.

The example in testLeakyWaves requires files from RandomJointEig, 
https://github.com/borplestenjak/RandomJointEig,
Bor Plestenjak, Haoze He, Hauke Gravenkamp, Daniel Kiefer 2024. 
To run this example, first download RandomJointEig and ensure that it is on the MATLAB path.




