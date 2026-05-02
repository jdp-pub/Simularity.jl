"""
    diag(A::AbstractMatrix)

# Arguments
- `A::AbstractMatrix`: The matrix of interest.

# Return 
The diagonal elements of a matrix as an array.

"""
diag(A::AbstractMatrix) = [A[xn,xn] for xn in 1:minimum(size(A))]

"""
    dot(a::AbstractArray,b::AbstractArray)

# Arguments
- `a::AbstractArray`: Array on left side of dot product.
- `b::AbstractArray`: Array on left side of dot product.

# Return 
The dot product of two arrays (vectors).

"""
dot(a::AbstractArray,b::AbstractArray) = sum(a.*b)

"""
    I(n::Int=2)


# Arguments
- `n::Int`: The length of one side of the matrix.

# Return 
nxn identity matrix.

"""
function I(n::Int=2)
    In = zeros(n,n)
    for ix in 1:n
        In[ix,ix] = 1
    end
    return In
end

"""
    dot(a::AbstractArray,b::AbstractArray)

# Arguments
- `a::AbstractArray`: Array on left side of dot product.
- `b::AbstractArray`: Array on left side of dot product.

# Return 
The dot product of two arrays (vectors).

# Description

# References
"""
function MBO(Ol,pos,l)
    # many body operator
    # pos and Ol should be sorted by index previously
    I = [1 0; 0 1]
    O = I
    k = 1
    if in(1,pos)
        O = Ol[1]
        k = k+1
    end
    
    for i in 2:l
        if in(i,pos)
            O = kron(O,Ol[k])
            k = k+1
        else
            O = kron(O,I)
        end
    end
    return O
end

"""
    dot(a::AbstractArray,b::AbstractArray)

# Arguments
- `a::AbstractArray`: Array on left side of dot product.
- `b::AbstractArray`: Array on left side of dot product.

# Return 
The dot product of two arrays (vectors).

# Description

# References

https://en.wikipedia.org/wiki/Norm_(mathematics)

"""
lpnorm(A,p::Number=2) = (sum(abs.(A).^p))^(1. /p)

"""
    dot(a::AbstractArray,b::AbstractArray)

# Arguments
- `a::AbstractArray`: Array on left side of dot product.
- `b::AbstractArray`: Array on left side of dot product.

# Return 
The dot product of two arrays (vectors).

# Description

# References
"""
normalize(x::AbstractVector,L::Int=2) = x/lpnorm(x,L)

"""
    dot(a::AbstractArray,b::AbstractArray)

# Arguments
- `a::AbstractArray`: Array on left side of dot product.
- `b::AbstractArray`: Array on left side of dot product.

# Return 
The dot product of two arrays (vectors).

# Description

# References
"""
tr(A::AbstractMatrix) = sum(diag(A))
