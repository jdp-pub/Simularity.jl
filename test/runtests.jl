include("Simularity.jl")
using Test

@testset "Simularity.jl" begin

    # Matrices
    
    ## diag()
    ### Square
    A = [1 0 0 0; 3 4 3 3; 4 4 5 4; 6 6 6 7]
    @test Simularity.diag(A) == [1, 4, 5, 7]
    
    ### Rectangular
    A = [1 0 0 0; 3 4 3 3; 4 4 5 4]
    @test Simularity.diag(A) == [1, 4, 5]

    A = [0 0 0; 4 3 3; 4 5 4; 6 6 7]
    @test Simularity.diag(A) == [0, 3, 4]

    ## dot()
    ### Vector dot
    x = [1,2,3]
    y = [4,5,6]
    @test Simularity.dot(x,y) == 32

    ### Matrix dot
    A = [1 1 1; 2 2 2; 3 3 3]
    B = [2 2 2; 3 3 3; 4 4 4]
    @test Simularity.dot(A,B) == 60

    # I()
    @test Simularity.I() == [1 0; 0 1]
    
    n = 5
    @test Simularity.I(n) == [1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1]
    
    # tr()
    ### Square
    A = [1 0 0 0; 3 4 3 3; 4 4 5 4; 6 6 6 7]
    @test Simularity.tr(A) == 17
    
    ### Rectangular
    A = [1 0 0 0; 3 4 3 3; 4 4 5 4]
    @test Simularity.tr(A) == 10

    A = [0 0 0; 4 3 3; 4 5 4; 6 6 7]
    @test Simularity.tr(A) == 7 
    
    
    # Matrix solvers
   
    ## eig_vals
    ### Diagonal
    A = [1 0; 0 2]
    #@test Simularity.eig_vals(A) == [1, 2]

    ### Nondiagonal
    A = [1 5; 3 2]
    #@test Simularity.eig_vals(A) == [1/2*(3+sqrt(61)), 1/2*(3-sqrt(61))]
    
    ## Power Iteration
    #A = [1.5 0.5; 0.5 -1.5]
    #out = Simularity.power_iteration(A)
    #@test isapprox(out[1],1) || isapprox(out[1],-1)
    
    ## QR decomposition
    ### Square
    A = [1 2; 3 4]
    Q,R = Simularity.qr_decomp(A)
    @test isapprox(Q,[sqrt(10)/10 3*sqrt(10)/10; 3*sqrt(10)/10 -sqrt(10)/10],atol=1E-10)
    @test isapprox(R,[sqrt(10) 7*sqrt(10)/5; 0 sqrt(10)/5],atol=1E-10)
    
    ### Rectangular
    A = [1 2 4; 3 4 5]
    Q,R = Simularity.qr_decomp(A)
    @test isapprox(Q,[sqrt(10)/10 3*sqrt(10)/10; 3*sqrt(10)/10 -sqrt(10)/10],atol=1E-10)
    @test isapprox(R,[sqrt(10) 7*sqrt(10)/5 19*sqrt(10)/10; 0 sqrt(10)/5 7*sqrt(10)/10],atol=1E-10)

    A = [1 2; 3 4; 5 6]
    Q,R = Simularity.qr_decomp(A)
    @test isapprox(Q,[sqrt(35)/35 13*sqrt(210)/210; 3*sqrt(35)/35 2*sqrt(210)/105; sqrt(35)/7 -sqrt(210)/42],atol=1E-10)
    @test isapprox(R,[sqrt(35) 44*sqrt(35)/35; 0 2*sqrt(210)/35],atol=1E-10)
    
    # Array sorts
    
    # small predetermined array
    A = [3,6,2,4]
    @test Simularity.bubble_sort(A) == [2,3,4,6]

    A = [8.45,9,78.6,34,5.76,4.39,8.57,6]
    @test Simularity.bubble_sort(A) == [4.39,5.76,6,8.45,8.57,9,34,78.6]

    A = [2,8,5,3,9,4]
    @test Simularity.insertion_sort(A) == [2,3,4,5,8,9]

    A = [2, 1, 4, 3]
    @test Simularity.merge_sort(A) == [1, 2, 3, 4]

    # timing large rangom arrays
    f = [Simularity.bubble_sort]
    n = 3 # length of random array  
    @test all((Simularity.vec_compare(f,n) .< [n^2/1E6]) .== 1)
end
