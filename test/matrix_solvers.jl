
"""
    arnoldi(A::AbstractMatrix,k::Int=size(A,2))

# Arguments
- `A::AbstractMatrix`: The matrix of interest.
- `k::Int`: The number of eigenvectors to consider in the subspace. 

# Return 
Upper Hessenberg form of A and the eivenvectors. 

# Description
Useful for non-Hermitian and Hermitian matrices.

"""
function arnoldi(A::AbstractMatrix,k::Int=size(A,2))
    # begin with ground state eigenvalue, eigenvector
    E,p = gs(A)
    v = []
    h = zeros(k,k)
    q = [p]

    # build nonsymmetric krylov space
    for m in 1:k
        v = A*q[m]
        for j in 1:m
            h[j,m] = q[j]'*v
            v = v - h[j,m]*q[j]
        end
        if m != k
            h[m+1,m] = lpnorm(v)
            push!(q,normalize(v)) # should be removed in favor of writing to preallocated space
        end
    end

    # eigenvectors
    q = reduce(hcat,q)

    return h,q
end

function eig_vals(A::AbstractMatrix;mode::String="qr",k::Int=1000,vtol::Number=1E-8)


    if mode == "qr"

        # should be changed to implicit qr algorihtm during optimization phase

        Ao = A
        Atemp = A*2
        kx = 0


        Q,R = qr_decomp(A)

        while !(isapprox(diag(Ao),diag(Atemp),atol=vtol) || kx==k)
            Atemp = Ao 
            Q,R = qr_decomp(Atemp)
            Ao = R*Q
            kx = kx+1
        end

        return diag(Ao)
    end
end

"""

# Arguments
- ``
- ``
- ``
- ``


# Return 

# Description


"""
function gs(A::AbstractMatrix,x::AbstractVector=complex.(rand(Float64,size(A,1))),vtol::Number=1E-8)
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
    lanczos(A::AbstractMatrix,k::Int=size(A,2))

# Arguments
- `A::AbstractMatrix`: The matrix of interest.
- `k::Int`: The number of eigenvectors to consider in the subspace. 

# Return 
Tridiagonal form of A and the eigenvectors.

# Description
Useful for hermitian matrices, faster than other methods in the valid case. 

"""
function lanczos(A::AbstractMatrix,k::Int=size(A,2))
    # begin with ground state
    E,p = gs(A)

    # krylov space for symmextric matrix
    # q columns are krylov subspace vectors
    q = [p]
    v = []
    h = zeros(k,k)
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

    return h,q
end

"""
    power_iteration(A::AbstractMatrix,x::AbstractVector=complex.(rand(Float64,size(A,1))),k::Int=100,vtol::Number=1E-6)

# Arguments
- `A::AbstractMatrix`: The matrix to be examined.
- `x::AbstractVector`: Guess for eigenvector, helps with convergence.
- `k::Int`: Maximum iterations for convergence
- `vtol::Float`: Convergence tolerance. 


# Return 
Largest magnitude complex eigenvalue with its unit eigenvector of A.

# Description
This is a stochastic method unless an initial guess is supplied. 
Best suited for a getting single eigen value/vector quickly.


"""
function power_iteration(A::AbstractMatrix,x::AbstractVector=complex.(rand(Float64,size(A,1))),k::Int=100,vtol::Number=1E-6)
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
    qr_decomp(A::AbstractMatrix,k::Int=size(A,2))

# Arguments
- `A::AbstractMatrix`: An invertible square matrix.

# Return 
A unitary matrix and an upper triangular matrix.

# Description
QR decomposition. Useful for performing higher level operations.

"""
function qr_decomp(A::AbstractMatrix)
    krow,kcol = size(A)

    # integers are cast to float, complex operation retained
    (eltype(A) isa Integer) ? type = typeof(float.(A)[1,1]) : type = typeof(float.(A)[1,1])

    krow < kcol ? Q = zeros(type, krow, krow) :  Q = zeros(type, krow, kcol)
    Q[:,1] = normalize(A[:,1])

    krow < kcol ? R = zeros(type, krow, kcol) : R =zeros(type, kcol, kcol)
    R[1,1] = lpnorm(A[:,1])

    v = Vector{Float64}(undef, krow)

    for kx in 2:kcol
        v = A[:,kx] 
        for m in 1:kx-1
            R[m,kx] = A[:,kx]'*Q[:,m]
            v = v - R[m,kx]*Q[:,m]
        end
        if kx <= kcol && kx <= krow
            R[kx,kx] = lpnorm(v)
            Q[:,kx] = normalize(v)
        end
    end
   return Q,R
end

function round_number!(A::AbstractMatrix;atol::Number=1E-18)
    x,y = size(A)
    for ix in 1:x
        for iy in 1:y
            if isnan() abs(A[x,y]) < atol 
                A[x,y] = 0
            end
        end
    end
end
