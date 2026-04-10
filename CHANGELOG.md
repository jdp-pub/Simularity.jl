# v0.1.0 ()
Initial feature release, introduces ideas of this 
package. No optimizations.

## Added
### Examples
- ising_gje.jl: shows the usage of MBO and ge functions.

### Source
- runtests.jl: for testing functions. 
- matrix_solvers.jl
- array_sorts.jl
- compare.jl
- matrices.jl 
    - MBO: Computes the many body operator given a set of interaction matrices.

# v0.2.0 (Unreleased)

## Added
### Examples
- OEO_rk4.jl: optoelectronic oscillator example using 4th order Runge-Kutta. 
- SHO_glrk.jl: Simple harmonic oscillator example using Gauss-Legendre Runge-Kutta.
- SHO_rk4.jl: Simple harmonic oscillator example using 4th order Runge-Kutta.

### Source
- diffeq_solvers.jl
    - rk1 (TODO)
    - rk2
    - rk3
    - rk4
    - rk45 (TODO)
    - glrk (TODO)

### Website (TODO)
A new website was made using Documenter.jl for code documentation.

### Tests
All functions were tested and tested again. (TODO)

## Changed
### Optimizations
All functions were verified in optimization tests and comparisons with modern standard. (TODO)


## Fixed
- Improper indexing in MBO. 
