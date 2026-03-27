include("Simularity.jl")
using Test

@testset "Simularity.jl" begin
    #### Matrix solvers
    
    #### Array sorts
    
    # small predetermined array
    A = [3,6,2,4]
    @test alg_tk.bubble_sort(A) == [2,3,4,6]

    A = [8.45,9,78.6,34,5.76,4.39,8.57,6]
    @test alg_tk.bubble_sort(A) == [4.39,5.76,6,8.45,8.57,9,34,78.6]


    # timing large rangom arrays
    f = [alg_tk.bubble_sort]
    n = 10 # length of random array  
    @test all((alg_tk.vec_compare(f,n) .< [n^2/1E6]) .== 1)
end
