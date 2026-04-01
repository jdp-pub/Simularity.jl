using LinearAlgebra
using SparseArrays

"""

# Arguments
- ``
- ``
- ``
- ``


# Return 

# Notes

# Examples
```julia


```
"""
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

"""

# Arguments
- ``
- ``
- ``
- ``


# Return 

# Notes

# Examples
```julia


```
"""
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

"""

# Arguments
- ``
- ``
- ``
- ``


# Return 

# Notes

# Examples
```julia


```
"""
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

"""
    power_iteration(A::AbstractMatrix{<:Number},x::AbstractVector{<:Number}=complex.(rand(Float64,size(A,1))),k::Int64=100,vtol::Float64=1E-6)

# Arguments
- `A::AbstractMatrix{<:Number}`: The matrix to be examined.
- `x::AbstractVector{<:Number}`: Guess for eigenvector, helps with convergence.
- `k::Int64`: Maximum iterations for convergence
- `vtol::Float64`: Convergence tolerance. 


# Return 
Largest magnitude complex eigenvalue with its unit eigenvector of A.

# Notes
This is a stochastic method unless an initial guess is supplied. 
Best suited for a getting single eigen value/vector quickly.

# Examples
```julia


```
"""
function power_iteration(A::AbstractMatrix{<:Number},x::AbstractVector{<:Number}=complex.(rand(Float64,size(A,1))),k::Int64=100,vtol::Float64=1E-6)
    kx = 0
    xtemp = x
    dottemp = NaN 
    while !(isapprox(dot(x,xtemp),dottemp,atol=vtol) || kx==k)
        dottemp = dot(x,xtemp)
        xtemp = x
        x = normalize(A*x) # action of A on x
        kx = kx+1
    end

    # rayleigh quotient
    # eigenevalue
    l = dot(A*x,x)/dot(x,x)
    return l, x
end

"""

# Arguments
- ``
- ``
- ``
- ``


# Return 

# Notes

# Examples
```julia


```
"""
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


