module alg_tk

########
#=
METHODS   
~~ matrix solvers
power_iteration(matrix, complex initial guess vector, number of iterations)
    returns dominant eigenvalue, eigenvector

gs(matrix, guess vector, variational tolerance)
    returns smallest real eigenvalue and corresponding eigencevtor to accuracy vtol

lanczos(symmetric matrix, length of side of output matrix)
    returns tridiagonal matrix, krylov vectors, gs eigenvalue 

arnoldi(square matrix, length of side of output, matrix)
    returns upper hessenberg matrix, krylov vectors, gs eigenvalue 

=#
#######


using LinearAlgebra
using SparseArrays

function power_iteration(A,x::Array{ComplexF64}=complex.(rand(Float64,size(A,1))),k::Int=11)
    # returns dominant unit eigenvector and eigenvalue of A
    # A can be complex
    # x initial is a real vector
    # returns complex eigenvalue and eigenvector
    
    # get the action of A
    for kx in 1:k
        # eigenvector
        x = normalize(A*x)
        display(x)
    end

    # rayleigh quotient
    # eigenevalue
    l = dot(A*x,x)/dot(x,x)

    return l, x
end

function gs(A,x::Array{ComplexF64}=complex.(rand(Float64,size(A,1))),vtol=1E-8)
    #findes the smallest eigenvalue and eigenvector 
    #using the variational method

    #variational method to find lowest eigenvector
    x = normalize(x)
    E = x'*A*x/(x'*x)
    p = (A*x -E*x)/(x'*x)
    p0 = zeros(size(p))
    E0 = 0.

    while !isapprox(real(E0),real(E),atol=vtol)
        p0 = p
        E0 = E
        E = (p0'*A*p0/(p0'*p0))
        p = p0 - (A*p0 - E0*p0)/(p0'*p0)
    end

    # ground state 
    return E, normalize(p)
end

function lanczos(A,k::Int=size(A,2))
    # approximates the tridiagonal form and 
    # eigenvectors of a hermitian matrix A. 

    # ground state
    E,p = gs(A)

    # krylov space for symmextric matrix
    # q columns are krylov subspace vectors
    q = [p]
    v = []
    h = spzeros(k,k)
    for kx in 1:k
        v = A*q[kx]
        h[kx,kx] = real(q[kx]'*v)
        if kx == 1
            v = v - h[kx,kx]*q[kx]
        else
            v = v - h[kx,kx]*q[kx] - h[kx,kx-1]*q[kx-1]
        end

        if kx != k
            h[kx,kx+1] = norm(v)
            h[kx+1,kx] = norm(v)
            push!(q,normalize(v))
        end
    end
    #eigenvectors
    q = reduce(hcat,q)

    return h,q,E
end

function arnoldi(A,k::Int=size(A,2))
    # works with nonsymmetric matrices, 
    # approximates the upper hessenberg form
    # of A

    # ground state eigenvalue, eigenvector
    E,p = gs(A)
    v = []
    h = spzeros(k,k)
    q = [p]

    # build nonsymmetric krylov space
    for m in 1:k
        v = A*q[m]
        for j in 1:m
            h[j,m] = q[j]'*v
            v = v - h[j,m]*q[j]
        end
        if m != k
            h[m+1,m] = norm(v)
            push!(q,normalize(v))
        end
    end
    # eigenvectors
    q = reduce(hcat,q)

    return h,q,E
end

function qr_decomp(A,k::Int=size(A,2))
    # A is invertible square matrix
    # returns Q, a unitary matrix, 
    # and R, an upper triangular matrix
    Q = [normalize(A[:,1])]
    R = spzeros(k,k)
    R[1,1] = norm(A[:,1])
    for kx in 2:k
        v = A[:,kx] 
        for m in 1:kx-1

            R[kx,m] = A[:,kx]'*Q[m]
            v = v - R[kx,m]*Q[m]
        end
        R[kx,kx] = norm(v)
        push!(Q,normalize(v))
    end

   return Q,R
end

end
