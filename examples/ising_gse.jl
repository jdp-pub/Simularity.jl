include("../src/Simularity.jl")
using .Simularity

function main()
    n = 6
    #n = 14 # takes a little over 1m 

    h = zeros(2^n,2^n)
    Z = [1 0; 0 -1]
    J = 0.6435

    # set up the system, hamiltonian for nearest neighbor model
    for nx in 1: n-1
        h = h + J*Simularity.MBO([Z,Z],[nx,nx+1],n)
    end

    # conduct computation with variational method in this case
    E, E_psi = Simularity.gs(h)

 
    println("Ground state: ", E_psi)
    println("Ground state energy: ", E)
end

main()
