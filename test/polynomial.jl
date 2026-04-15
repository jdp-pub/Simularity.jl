lagrange_basis(xl::AbstractArray,x::Number) = [lagrange_poly(xl,x,j) for j in 1:length(xl)]


lagrange_interpolation(xl::AbstractArray,yl::AbstractArray,x::Number) = sum(yl .* lagrange_basis(xl,x))


function lagrange_poly(xl::AbstractArray,x::Number,j::Int)
    k = length(xl)
    l = 0
    for m in 0:k
        if m != j
            l = l * (x-xl[m])/(xl[j]-xl[m])
        end
    end
    return l
end


normalized_legendre_poly(s::Int) = [factorial(s)^2/factorial(2*s)*(-1)^(s-k)*binomial(s,k)*binomial(s+k,k) for k in 0:s]


polyroots(p::Vector) = eig_vals(hcat(I(length(p)-1)[:,2:length(p)-1],-p[1:length(p)-1]))

