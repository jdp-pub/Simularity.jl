"""

# References
https://sites.math.rutgers.edu/~falk/math573/lecture13.pdf
"""
function adaptive()

end

"""
    riemann_integration(x::AbstractArray,y::AbstractArray,mode::String="mid")

# Arguments
- `x::AbstractArray`: The horizontal discretization. 
- `y::AbstractArray`: The output of the map f at position x_i; y_i = f(x_i). 
- `mode::String`: Options: left, right, mid. Determines rectangle placement.

# Returns
A numerical approximation of the integration of a function y on some interval given in x.

# Description
Integration using Riemann sums[^riemann-sums]. In practice, the midpoint rule becomes the trapezoidal rule[^the-midpoint-and-trapezoidal-rules] since this is a discrete space. This is apparent when one uses linear interpolation in the midpoint rule to approximate the midpoint f((x_1+x_2)/2)->(f(x_1)+f(x_2))/2. Left and right rules included.

# References
[^riemann-sums]: [Riemann Sums, https://math.libretexts.org/Bookshelves/Calculus/Calculus_3e_(Apex)/05%3A_Integration/5.03%3A_Riemann_Sums (accessed April 15, 2026).](https://math.libretexts.org/Bookshelves/Calculus/Calculus_3e_(Apex)/05%3A_Integration/5.03%3A_Riemann_Sums)
[^the-midpoint-and-trapezoidal-rules]: [The Midpoint and Trapezoidal Rules,  https://courses.lumenlearning.com/calculus2/chapter/the-midpoint-and-trapezoidal-rules/ (accessed April 15, 2026).](https://courses.lumenlearning.com/calculus2/chapter/the-midpoint-and-trapezoidal-rules/)
"""
function riemann_integration(x::AbstractArray,y::AbstractArray,mode::String="mid")
    if mode == "mid"
        return sum([(x[nx+1]-x[nx])*(y[nx]+y[nx+1])/2 for nx in 1:length(x)-1])
    elseif mode == "left"
        return sum([(x[nx+1]-x[nx])*y[nx] for nx in 1:length(x)-1])
    elseif mode == "right"
        return sum([(x[nx+1]-x[nx])*y[nx+1] for nx in 1:length(x)-1])
    end
end

function monte_carlo()

end




