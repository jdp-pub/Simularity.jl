include("Simularity.jl")
using Test

@testset "Simularity.jl" begin
    #### Matrix solvers
    A = [1.5 0.5; 0.5 -1.5]
    out = Simularity.power_iteration(A)
    @test isapprox(out[1],1) || isapprox(out[1],-1) 
    #[-1 0]))
    
    #### Array sorts
    
    # small predetermined array
    A = [3,6,2,4]
    @test Simularity.bubble_sort(A) == [2,3,4,6]

    A = [8.45,9,78.6,34,5.76,4.39,8.57,6]
    @test Simularity.bubble_sort(A) == [4.39,5.76,6,8.45,8.57,9,34,78.6]


    # timing large rangom arrays
    f = [Simularity.bubble_sort]
    n = 3 # length of random array  
    @test all((Simularity.vec_compare(f,n) .< [n^2/1E6]) .== 1)
end
