
diag(A::AbstractMatrix{<:Number}) = [A[xn,xn] for xn in 1:minimum(size(A))]
dot(a::AbstractArray{<:Number},b::AbstractArray{<:Number}) = sum(a.*b)

function I(n::Int=2)
    In = zeros(n,n)
    for ix in 1:n
        In[ix,ix] = 1
    end
    return In
end

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

function lpnorm(A,L::Int=2)
    if L == 1

    elseif L == 2
        return sqrt(sum(A.*A))


    elseif L == 3

    end
end


normalize(x::AbstractVector{<:Number},L::Int=2) = x/lpnorm(x,L)
tr(A::AbstractMatrix{<:Number}) = sum(diag(A))
