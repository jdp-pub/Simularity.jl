"""
    diag(A::AbstractMatrix)

# Arguments
- `A::AbstractMatrix`: The matrix of interest.

# Return 
The diagonal elements of a matrix as an array[^Matrix_(mathematics)].

# References
[^Matrix_(mathematics)]: [Matrix (Mathematics), https://en.wikipedia.org/wiki/Matrix_(mathematics) (May 05, 2026).](Matrix_(mathematics))

"""
diag(A::AbstractMatrix) = [A[xn,xn] for xn in 1:minimum(size(A))]

"""
    dot(a::AbstractArray,b::AbstractArray)

# Arguments
- `a::AbstractArray`: Array on left side of dot product.
- `b::AbstractArray`: Array on left side of dot product.

# Return 
The dot product of two arrays[^Dot_product].

# References
[^Dot_product]: [Dot Product, https://en.wikipedia.org/wiki/Dot_product (May 05, 2026).](https://en.wikipedia.org/wiki/Dot_product)

"""
dot(a::AbstractArray,b::AbstractArray) = sum(a.*b)

"""
    I(n::Int=2)


# Arguments
- `n::Int`: The length of one side of the matrix.

# Return 
nxn identity matrix[^Matrix_(mathematics)].

# References
[^Matrix_(mathematics)]: [Matrix (Mathematics), https://en.wikipedia.org/wiki/Matrix_(mathematics) (May 05, 2026).](Matrix_(mathematics))

"""
function I(n::Int=2)
    In = zeros(n,n)
    for ix in 1:n
        In[ix,ix] = 1
    end
    return In
end

"""
    MBO(Ol::AbstractArray,pos::AbstractArray,l::Int)

# Arguments
- `Ol::AbstractArray`: List of operators in the coupling. Should be sorted by index 
- `pos::AbstractArray`: The positions of items corresponting to the operators in Ol. Should be sorted by index 
- `l::Int`: The total number of items in the system the operator acts on.

# Return 
Many body operator as a matrix[^Second_quantization].

# References
[^Second_quantization]: [Second Quantization, https://en.wikipedia.org/wiki/Second_quantization (May 05, 2026).](https://en.wikipedia.org/wiki/Second_quantization)

"""
function MBO(Ol::AbstractArray,pos::AbstractArray,l::Int)
    # many body operator
    # pos and Ol previously
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
    lpnorm(A::AbstractArray,p::Number=2)

# Arguments
- `A::AbstractArray`: The array that produces the norm.
- `p::Number=2`: Degree of the norm. 

# Return 
The norm of the array[^Norm_(mathematics)].


# References
[^Norm_(mathematics)]: [Norm (Mathematics), https://en.wikipedia.org/wiki/Norm_(mathematics) (May 05, 2026).](https://en.wikipedia.org/wiki/Norm_(mathematics))

"""
lpnorm(A::AbstractArray,p::Number=2) = (sum(abs.(A).^p))^(1. /p)

"""
    normalize(x::AbstractVector,p::Int=2)

# Arguments
- `x::AbstractArray`: The array to normalize.
- `p::Number`: The degree of the norm to use.

# Return 
A normalized array[^Matrix_(mathematics)].

# References
[^Matrix_(mathematics)]: [Matrix (Mathematics), https://en.wikipedia.org/wiki/Matrix_(mathematics) (May 05, 2026).](Matrix_(mathematics))

"""
normalize(x::AbstractVector,p::Int=2) = x./lpnorm(x,p)

"""
    tr(A::AbstractMatrix)

# Arguments
- `A::AbstractMatrix`: Array on left side of dot product.

# Return 
The trace of matrix A[^Matrix_(mathematics)].

# References
[^Matrix_(mathematics)]: [Matrix (Mathematics), https://en.wikipedia.org/wiki/Matrix_(mathematics) (May 05, 2026).](Matrix_(mathematics))

"""
tr(A::AbstractMatrix) = sum(diag(A))
