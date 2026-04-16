# Introduction

This document records project implemented changes, planned changes, and potential directions. 

The version history follows [Semantic Versioning](https://semver.org/). 

Release dates are mm/dd/yyyy format.

## Ideas for the future

### Adjustments
- Adaptive time discretization in Runge-Kutta methods.
- Adaptive integration in methods that require it.
- QR algorithm should be made implicitly shifted.
- relative and abaolute tolerances.

### Compute
- Adjustable CPU/GPU multiprocessing routines/wrappers.
- pretty much everything will need to be adjusted to ensure competive speed. 
- "!" function versions, mutate state rather than make copies.
- Support for sparse matrices and arrays.

### Features
- Adaptive integration.
- Graph sorts.
- Optimal path solvers.
- Complete eigenvalue implementation (similar to LAPACK,ARPACK functionality).
- Reimplement some LAPACK,ARPACK functions natively, see if native can compete with wrapper implementations.
- Tensor network solvers.


---
Begin Version History
---

# v0.2.0 (Unreleased)
This update focuses on expanding core functionality; keeping scalability in mind. Future updates will adress optomizations to make the package competetive. This has not been tested for breaks with v0.1.0. Breaks should be expected.

## Added

### Examples
- OEO_rk4.jl: optoelectronic oscillator example using 4th order Runge-Kutta. 
- SHO_glrk.jl: Simple harmonic oscillator example using Gauss-Legendre Runge-Kutta.
- SHO_rk4.jl: Simple harmonic oscillator example using 4th order Runge-Kutta.

### Source
- array_sorts.jl
    - bubble_sort: allow descending order sorting (TODO), added partial sorting.
    - insertion_sort
    - merge_sort
- diffeq_solvers.jl
    - rk1
    - rk2
    - rk23 (TODO)
    - rk3
    - rk4
    - rk45 (TODO)
    - glrk: Gauss-Legendre Runge Kutta. (TODO)
- matrices.jl
    - diag
    - dot
    - I
    - lpnorm: 
        - L1 norm (TODO)
        - L2 norm
        - L3 norm (TODO)
    - normalize
    - tr
- matrix_solvers.jl
    - eig_vals: Returns eigenvalues of a matrix using QR.   
        - add other methods, arnoldi, lanczos (TODO)
    - qr_decomp: Accepts rectangular matrices.
    - round_number!: Rounds elements of a matrix to 0.
- polynomial.jl
    - lagrange_basis
    - lagrange_interpolation
    - lagrange_poly
    - normalized_legendre
    - polyroots

### Test
- test.sh: test modularization with flags
    - ability to run all tests
    - run all tests with code coverage
    - add modular testing (TODO)

### Website (TODO)
A new website was made using Documenter.jl for code documentation.
    - link CHANGELOG.md to website (TODO).

### Tests
All functions were tested and tested again. (TODO)

## Changed

- matrix_solvers: 
    - arnoldi: Defined expected types.
    - gs: Defined expected types.
    - lanczos: Defined expected types.
    - power_iteration: Defined expected types. Changed convergence type from iteration count to tolerance check.
    - qr_decomp: Defined expected types. Changed convergence type from iteration count to tolerance check.

## Fixed
- matrices.jl
    - Improper indexing in MBO. 

---

# v0.1.0 (4/14/2026)
Initial feature release, introduces ideas of this 
package.

## Added
### Examples
- ising_gse.jl: Shows the usage of MBO and ge functions.

### Source
- array_sorts.jl
    - bubble_sort: Sorts an array in O(n^2), fast  partial sorts. Only in descending order at the moment.
- compare.jl
    - vec_compare: Compares time taken for 2 functions that operate on vectors.
- matrices.jl 
    - MBO: Computes the many body operator given a set of interaction matrices.
- matrix_solvers.jl
    - arnoldi: Approximates the upper-Hessenberg form of a general matrix.
    - lanczos: Approximates the tridiagonal form of a hermitian matrix.
    - gs: uses Arnoldi to compute ground state energy and corresponding eigenvector of a given Hamiltonian matrix.
    - power_iteration: Computes maximum amplitude eigenvalue and corresponding eigenvector.
    - qr_decomp: QR decomposition for square matrices.
- Simularity.jl

### Test
 - runtests.jl
    - Full program testing suite.


