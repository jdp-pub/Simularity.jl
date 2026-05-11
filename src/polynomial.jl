
"""
    lagrange_basis(xl::AbstractArray,x::Number) 

# Arguments
- `xl::AbstractArray`: Discrete points in function domain. The set of nodes for the space.
- `x::Number`: The position in the domain.

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
- `xl::AbstractArray`: Discrete points in function domain. The set of nodes for the space.
- `yl::AbstractArray`: The output of some function acting on the nodes; y_i = f(x_i).
- `x::Number`: The position in the domain.

# Return 
The output of f(x_m) for some x_m that lies between nodes in xl.

# Description
Lagrange interpolation method[^Lagrange_polynomial]. The output of f(x_m) is from a
minimum order polynomial needed to build the function f at the point of interpolation.

# References
[^Lagrange_polynomial]: [Lagrange Polynomial, https://en.wikipedia.org/wiki/Lagrange_polynomial (accessed April 16, 2026).](https://en.wikipedia.org/wiki/Lagrange_polynomial)

"""
lagrange_interpolation(xl::AbstractArray,yl::AbstractArray,x::Number) = sum(yl .* lagrange_basis(xl,x))

"""
    lagrange_poly(xl::AbstractArray,x::Number,j::Int)

# Arguments
- `xl::AbstractArray`: Discrete points in function domain. The set of nodes for the space.
- `x::Number`: The position in the domain.
- `j::Int`: The order of the polynomial.

# Return 
The Lagrange polynomial of order j at position x over domain xl.

# Description
The Lagrange polynomial used for building a basis which allows 
    efficient solution to a wide class of problems[^Lagrange_polynomial].

# References
[^Lagrange_polynomial]: [Lagrange Polynomial, https://en.wikipedia.org/wiki/Lagrange_polynomial (accessed April 16, 2026).](https://en.wikipedia.org/wiki/Lagrange_polynomial)

"""
function lagrange_poly(xl::AbstractArray,x::Number,j::Int)
    k = length(xl)
    l = 1
    for m in 1:k
        if m != j
            l = l * (x-xl[m])/(xl[j]-xl[m])
        end
    end
    return l
end

"""
    normalized_legendre_poly(s::Int,x::Number)

# Arguments
- `s::Int`: The order of the polynomial.
- `x::Number`: The order of the polynomial.

# Return 
The Legendre polynomial of order s evaluated at x.

# Description
The Legendre polynomials are complete and orthogonal[^Legendre_polynomials]. 
The Legendre polynomials can be expressed in many ways[^Legendre_polynomials] but this specific
form was acquired from Melvin Leok[^Gauss-Legendre-method-as-a-collocation-method], although.

# References
[^Legendre_polynomials]: [Legendre Polynomials, https://en.wikipedia.org/wiki/Legendre_polynomials (accessed April 17, 2026).](https://en.wikipedia.org/wiki/Legendre_polynomials)
[^Gauss-Legendre-method-as-a-collocation-method]: [Gauss-Legendre method as a collocation method,  https://www.youtube.com/watch?v=o_vtNxWNnUQ&t=13s (accessed April 17, 2026).](https://www.youtube.com/watch?v=o_vtNxWNnUQ&t=13s)
"""
normalized_legendre_poly(s::Int,x::Number) = sum([factorial(s)^2/factorial(2*s)*(-1)^(s-k)*binomial(s,k)*binomial(s+k,k)*x^k for k in 0:s])

"""
    normalized_legendre_poly_coef(s::Int)

# Arguments
- `s::Int`: The order of the polynomial.

# Return 
The Legendre polynomial coefficients in a vector in ascending order. 
Symbollicly; [c0_x^0 c_1x^1 ... c_sx^s] with coefficient c for variable x up to order s.

# Description
The Legendre polynomial coefficients which are coefficients of complete and orthogonal polynomials[^Legendre_polynomials]. 
The Legendre polynomials can be expressed in many ways[^Legendre_polynomials] but this specific
form was acquired from Melvin Leok[^Gauss-Legendre-method-as-a-collocation-method].


# References
[^Legendre_polynomials]: [Legendre Polynomials, https://en.wikipedia.org/wiki/Legendre_polynomials (accessed April 17, 2026).](https://en.wikipedia.org/wiki/Legendre_polynomials)
[^Gauss-Legendre-method-as-a-collocation-method]: [Gauss-Legendre method as a collocation method,  https://www.youtube.com/watch?v=o_vtNxWNnUQ&t=13s (accessed April 17, 2026).](https://www.youtube.com/watch?v=o_vtNxWNnUQ&t=13s)
"""
normalized_legendre_poly_coef(s::Int) = [factorial(s)^2/factorial(2*s)*(-1)^(s-k)*binomial(s,k)*binomial(s+k,k) for k in 0:s]

"""
    polyroots(p::AbstractVector)

# Arguments
- `p::AbstractVector`: The coefficients of the polynomial in ascending order.

# Return 
The roots of a polynomial as a vector.

# Description
Finds the roots of a polynomial given the coefficients by 
finding the eigenvalues of the companion matrix of the
polynomial[^Companion_matrix].

# References
[^Companion_matrix]: [Companion Matrix, https://en.wikipedia.org/wiki/Companion_matrix (accessed April 17, 2026).](https://en.wikipedia.org/wiki/Companion_matrix)

"""
polyroots(p::AbstractVector) = eig_vals(hcat(I(length(p)-1)[:,2:length(p)-1],-p[1:length(p)-1]))

