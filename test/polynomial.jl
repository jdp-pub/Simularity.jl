
"""
    lagrange_basis(xl::AbstractArray,x::Number) 

# Arguments
- `xl`: Discrete points in function domain. The set of nodes for the space.
- `x`: The position in the domain.

# Return 
The Lagrange basis [l_0(x),l_1(x),... l_k(x)] where k is the number of nodes, l_j is the jth lagrange polynomial.

# Description
Builds the basis needed for interpolation[^Lagrange_polynomial].

# References
[^Lagrange_polynomial]: [Lagrange Polynomial, https://en.wikipedia.org/wiki/Lagrange_polynomial (accessed April 16, 2026).](https://en.wikipedia.org/wiki/Lagrange_polynomial)

"""
lagrange_basis(xl::AbstractArray,x::Number) = [lagrange_poly(xl,x,j) for j in 1:length(xl)]

"""
    lagrange_interpolation(xl::AbstractArray,yl::AbstractArray,x::Number)

# Arguments
- `xl`: Discrete points in function domain. The set of nodes for the space.
- `yl`: The output of some function acting on the nodes; y_i = f(x_i).
- `x`: The position in the domain.

# Return 
The output of f(x_m) for some x_m that lies between nodes in xl.

# Description
Lagrange interpolation method[^Lagrange_polynomial]. The output of f(x_m) is from a
minimum order polynomial needed to build the function f.

# References
[^Lagrange_polynomial]: [Lagrange Polynomial, https://en.wikipedia.org/wiki/Lagrange_polynomial (accessed April 16, 2026).](https://en.wikipedia.org/wiki/Lagrange_polynomial)

"""
lagrange_interpolation(xl::AbstractArray,yl::AbstractArray,x::Number) = sum(yl .* lagrange_basis(xl,x))

"""
    lagrange_poly(xl::AbstractArray,x::Number,j::Int)

# Arguments
- `xl`: Discrete points in function domain. The set of nodes for the space.
- `x`: The position in the domain.
- `j`: The order of the polynomial.

# Return 
The Lagrange polynomial of order j at position x over domain xl.

# Description
The Lagrange polynomial[^Lagrange_polynomial]. Used in many important methods[^Lagrange_polynomial].

# References
[^Lagrange_polynomial]: [Lagrange Polynomial, https://en.wikipedia.org/wiki/Lagrange_polynomial (accessed April 16, 2026).](https://en.wikipedia.org/wiki/Lagrange_polynomial)

"""
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

"""
    normalized_legendre_poly(s::Int)

# Arguments
- `s`: The order of the polynomial.

# Return 
The Legendre polynomial coefficients in a vector in ascending order. 
Symbollicly; [c0_x^0 c_1x^1 ... c_sx^s] with coefficient c for variable x up to order s.

# Description

# References

"""
normalized_legendre_poly(s::Int) = [factorial(s)^2/factorial(2*s)*(-1)^(s-k)*binomial(s,k)*binomial(s+k,k) for k in 0:s]


polyroots(p::Vector) = eig_vals(hcat(I(length(p)-1)[:,2:length(p)-1],-p[1:length(p)-1]))

