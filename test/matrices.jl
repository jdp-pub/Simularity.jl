"""
    diag(A::AbstractMatrix)

# Arguments
- `A::AbstractMatrix`: The matrix of interest.

# Return 
The diagonal elements of a matrix as an array.

"""
diag(A::AbstractMatrix) = [A[xn,xn] for xn in 1:minimum(size(A))]

"""
    dot(a::AbstractArray{<:Number},b::AbstractArray{<:Number})

# Arguments
- `a::AbstractArray{<:Number}`: Array on left side of dot product.
- `b::AbstractArray{<:Number}`: Array on left side of dot product.

# Return 
The dot product of two arrays (vectors).

"""
dot(a::AbstractArray{<:Number},b::AbstractArray{<:Number}) = sum(a.*b)

"""
    I(n::Int=2)


# Arguments
- `a::AbstractArray{<:Number}`: Array on left side of dot product.
- `b::AbstractArray{<:Number}`: Array on left side of dot product.

# Return 
The dot product of two arrays (vectors).

# Description

# References
"""
function I(n::Int=2)
    In = zeros(n,n)
    for ix in 1:n
        In[ix,ix] = 1
    end
    return In
end

"""
    dot(a::AbstractArray{<:Number},b::AbstractArray{<:Number})

# Arguments
- `a::AbstractArray{<:Number}`: Array on left side of dot product.
- `b::AbstractArray{<:Number}`: Array on left side of dot product.

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
    dot(a::AbstractArray{<:Number},b::AbstractArray{<:Number})

# Arguments
- `a::AbstractArray{<:Number}`: Array on left side of dot product.
- `b::AbstractArray{<:Number}`: Array on left side of dot product.

# Return 
The dot product of two arrays (vectors).

# Description

# References
"""
function lpnorm(A,L::Int=2)
    if L == 1

    elseif L == 2
        return sqrt(sum(A.*A))


    elseif L == 3

    end
end

"""
    dot(a::AbstractArray{<:Number},b::AbstractArray{<:Number})

# Arguments
- `a::AbstractArray{<:Number}`: Array on left side of dot product.
- `b::AbstractArray{<:Number}`: Array on left side of dot product.

# Return 
The dot product of two arrays (vectors).

# Description

# References
"""
normalize(x::AbstractVector{<:Number},L::Int=2) = x/lpnorm(x,L)

"""
    dot(a::AbstractArray{<:Number},b::AbstractArray{<:Number})

# Arguments
- `a::AbstractArray{<:Number}`: Array on left side of dot product.
- `b::AbstractArray{<:Number}`: Array on left side of dot product.

# Return 
The dot product of two arrays (vectors).

# Description

# References
"""
tr(A::AbstractMatrix{<:Number}) = sum(diag(A))
