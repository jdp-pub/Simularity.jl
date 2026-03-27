module Simularity

#=
~~ matrix solvers
arnoldi(square matrix, length of side of output, matrix)
    returns upper hessenberg matrix, krylov vectors, gs eigenvalue 

gs(matrix, guess vector, variational tolerance)
    returns smallest real eigenvalue and corresponding eigencevtor to accuracy vtol

lanczos(symmetric matrix, length of side of output matrix)
    returns tridiagonal matrix, krylov vectors, gs eigenvalue 

power_iteration(matrix, complex initial guess vector, number of iterations)
    returns dominant eigenvalue, eigenvector
=#
include("matrix_solvers.jl")

#=
~~ array sorts
bubble_sort(vector)
    returns sorted vector O(n^2)

quick_sort(vector)
    returns sorted vector O(nlog(n))

=#
include("array_sorts.jl")

#=
~~ compare
vec_compare(vector of functions, length of random vectors as input, number of runs to average per function)
    returns vector of average times

=#
include("compare.jl")

include("matrices.jl")
end
