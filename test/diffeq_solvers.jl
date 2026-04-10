
"""
    glrk(f,y0::Vector{<:Number},ti::Number=0,tf::Number=10,n::Int=1000,kn::Int=10,fargs::Vector{Any}=[])

# Arguments
- `f`: Function that describes dynamical system. 
- `y0::Vector{<:Number}`: Initial conditions.
- `ti::Number`: Initial time.
- `tf::Number`: End time.
- `n::Int`: The number of timesteps.
- `s::Int`: Runge-Kutta stages.
- `fargs::Vector{Any}`: Additional parameters to pass to f.


# Return 
The time evolution of supplied parameters and the corresponding time series.

# Description
Gauss-Legendre Runge-Kutta ODE solver. Good for nth order implicit and explicit 
Runge-Kutta, dynamical evolution of systems that can be cast as an array of ODEs.

# References
"""
function glrk(f,y0::Vector{<:Number},ti::Number=0,tf::Number=10,n::Int=1000,s::Int=10,fargs::Vector=[])

    y = y0
    yt = y
    
    yl = Vector{typeof(y0)}(undef,n)
    yl[1] = y0

    t = range(ti,tf,n)
    dt = t[2]-t[1]
    # set up polynomial in ascending order, index is polynomial power
    p = normalized_legendre(s) 
    
    # find roots of polynomial
    roots = polyroots(p)

    # integrate roots

    stop

    # compute runge kutta
    for nx in 2:n


        yl[nx] = y
    end

    ytl = [[row[i] for row in yl] for i in 1:length(yl[1])]
    return ytl,t
end

"""
    rk2(f,y0::Vector{<:Number},ti::Number=0,tf::Number=10,n::Int=1000,fargs::Vector{Any}=[])


# Arguments
- `f`: Function that describes dynamical system. 
- `y0::Vector{<:Number}`: Initial conditions.
- `ti::Number`: Initial time.
- `tf::Number`: End time. 
- `n::Int`: The number of timesteps.
- `fargs::Vector{Any}`: Additional parameters to pass to f.

# Return 
The time evolution of supplied parameters and the corresponding time series.

# Description
Second order Runge-Kutta ODE solver, dynamical evolution of systems that can be cast 
as an array of ODEs.

"""
function rk2(f,y0::Vector{<:Number},ti::Number=0,tf::Number=10,n::Int=1000,fargs::Vector{Any}=[])

    y = y0
    yt = y
    
    yl = Vector{typeof(y0)}(undef,n)
    yl[1] = y0

    t = range(ti,tf,n)
    dt = t[2]-t[1]

    k1 = typeof(y0)(undef,n)
    k2 = typeof(y0)(undef,n)

    for nx in 2:n
        k1 = f(y,t[nx],fargs).*dt

        yt = y + k1/2
        k2 = f(yt,t[nx]+dt/2,fargs).*dt

        y = y + (k1 + k2)/2

        yl[nx] = y
    end

    ytl = [[row[i] for row in yl] for i in 1:length(yl[1])]
    return ytl,t
end


"""
    rk3(f,y0::Vector{<:Number},ti::Number=0,tf::Number=10,n::Int=1000,fargs::Vector{Any}=[])


# Arguments
- `f`: Function that describes dynamical system. 
- `y0::Vector{<:Number}`: Initial conditions.
- `ti::Number`: Initial time.
- `tf::Number`: End time. 
- `n::Int`: The number of timesteps.
- `fargs::Vector{Any}`: Additional parameters to pass to f.

# Return 
The time evolution of supplied parameters and the corresponding time series.

# Description
Third order Runge-Kutta ODE solver, dynamical evolution of systems that can be cast 
as an array of ODEs.

"""
function rk3(f,y0::Vector{<:Number},ti::Number=0,tf::Number=10,n::Int=1000,fargs::Vector{Any}=[])

    y = y0
    yt = y
    
    yl = Vector{typeof(y0)}(undef,n)
    yl[1] = y0

    t = range(ti,tf,n)
    dt = t[2]-t[1]

    k1 = typeof(y0)(undef,n)
    k2 = typeof(y0)(undef,n)
    k3 = typeof(y0)(undef,n)

    for nx in 2:n
        k1 = f(y,t[nx],fargs).*dt

        yt = y + k1/2
        k2 = f(yt,t[nx]+dt/2,fargs).*dt

        yt = y + 2*k2 - k1
        k3 = f(yt,t[nx]+dt,fargs).*dt


        y = y + (k1 + 4k2 + k3)/6

        yl[nx] = y
    end

    ytl = [[row[i] for row in yl] for i in 1:length(yl[1])]
    return ytl,t
end

"""
    rk4(f,y0::Vector{<:Number},ti::Number=0,tf::Number=10,n::Int=1000,fargs::Vector{Any}=[])


# Arguments
- `f`: Function that describes dynamical system. 
- `y0::Vector{<:Number}`: Initial conditions.
- `ti::Number`: Initial time.
- `tf::Number`: End time. 
- `n::Int`: The number of timesteps.
- `fargs::Vector{Any}`: Additional parameters to pass to f.

# Return 
The time evolution of supplied parameters and the corresponding time series.

# Description
Fourth order Runge-Kutta ODE solver, dynamical evolution of systems that can be cast 
as an array of ODEs.

"""
function rk4(f,y0::Vector{<:Number},ti::Number=0,tf::Number=10,n::Int=1000,fargs::Vector{Any}=[])

    y = y0
    yt = y
    
    yl = Vector{typeof(y0)}(undef,n)
    yl[1] = y0

    t = range(ti,tf,n)
    dt = t[2]-t[1]

    k1 = typeof(y0)(undef,n)
    k2 = typeof(y0)(undef,n)
    k3 = typeof(y0)(undef,n)
    k4 = typeof(y0)(undef,n)

    for nx in 2:n
        k1 = f(y,t[nx],fargs).*dt

        yt = y + k1/2
        k2 = f(yt,t[nx]+dt/2,fargs).*dt

        yt = y + k2/2
        k3 = f(yt,t[nx]+dt/2,fargs).*dt

        yt = y + k3
        k4 = f(yt,t[nx]+dt,fargs).*dt

        y = y + (k1 + 2*k2 + 2*k3 + k4)/6

        yl[nx] = y
    end

    ytl = [[row[i] for row in yl] for i in 1:length(yl[1])]
    return ytl,t
end


