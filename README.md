# Sparse linear systems solvers
Implementation of different iterative methods to solve large sparse linear systems of equations written in matrix form, i.e., Ax=b.

### Sequential Monte Carlo solvers
Iterative methods that uses Monte Carlo to estimate the inner product of the solution with a weighting vector.
- Monte Carlo Almost optimal (MAO) algorithm introduced by Dimov.
<br />Reference: Dimov, I., & Karaivanova, A. (1997). Iterative Monte Carlo algorithms for linear algebra problems. Numerical Analysis and Its Applications, 150-160.
- A new faster algorithm I developed during my research.

### Basic iterative solvers
Basic iterative methods studied in 'Advanced Matrix Algorithims' graduate course.
<br />Reference: Saad, Y. (2003). Iterative methods for sparse linear systems (Vol. 82). siam.
- Block Gauss Seidel method
- Block Jacobi method

### Generalized Minimum Residual (GMRES) solver 
Iterative projection method studied in 'Advanced Matrix Algorithims' graduate course.
<br />References: Saad, Y. (2003). Iterative methods for sparse linear systems (Vol. 82). siam.

### Permutation functions
Permutation (reordering) methods studied in 'Advanced Matrix Algorithims' graduate course to improve the convergence of iterative solvers.
<br />Reference: Saad, Y. (2003). Iterative methods for sparse linear systems (Vol. 82). siam.
- Reverse Cuthill-McKee method
- Breadth First Search method

