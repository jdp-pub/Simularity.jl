"""
    glrk(f,y0::AbstractVector,ti::Number=0,tf::Number=10,n::Int=1000,kn::Int=10,fargs::AbstractVector=[])

# Arguments
- `f`: Function that describes dynamical system. 
- `y0::AbstractVector`: Initial conditions.
- `ti::Number`: Initial time.
- `tf::Number`: End time.
- `n::Int`: The number of timesteps.
- `s::Int`: Runge-Kutta stages.
- `fargs::AbstractVector`: Additional parameters to pass to f.


# Return 
The time evolution of supplied parameters and the corresponding time series.

# Description
Gauss-Legendre Runge-Kutta ODE solver[^list-of-runge-kutt-methods]. Good for nth order implicit and explicit 
Runge-Kutta, dynamical evolution of systems that can be cast as an array of ODEs.

# References
[^list-of-runge-kutt-methods]: [List of Runge-Kutta Methods, https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods (accessed April 14, 2026).](https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods)

"""
function glrk(f,y0::AbstractVector,ti::Number=0,tf::Number=10,n::Int=1000,s::Int=10,fargs::AbstractVector=[],xn::Int=100000)

    # needs stability tests and optomization
    

    y = y0
    yt = y
    
    yl = Vector{typeof(y0)}(undef,n)
    yl[1] = y0

    t = range(ti,tf,n)
    dt = t[2]-t[1]

    # calculating the butcher-tableau may need to be moved to its own function
    # set up polynomial in ascending order, index is polynomial power
    p = normalized_legendre_poly_coef(s) 

    
    # find roots of polynomial
    roots = polyroots(p)

    aij = zeros(length(roots),length(roots))
    bi = zeros(length(roots))
    lj = zeros(xn)
    li = zeros(xn)

    for i in 1:length(roots)
        for j in 1:length(roots)
            lj = [lagrange_poly(roots,x,j) for x in range(0,roots[i],xn)]
            
            # try higher accuracy integrators to reduce iteration count
            aij[i,j] = riemann_integration(range(0,roots[i],xn),lj)
        end
        li = [lagrange_poly(roots,x,i) for x in range(0,1,xn)]
        bi[i] = riemann_integration(range(0,1,xn),li)
    end


    # compute runge kutta
    for nx in 2:n
        yi = [y + sum([aij[i,j]*f(y,t[nx],fargs) for j in 1:length(roots)])*dt for i in 1:length(roots)]
        y = y + sum([bi[i]*f(y,t[nx],fargs) for i in 1:length(roots)])*dt

        yl[nx] = y
    end

    ytl = [[row[i] for row in yl] for i in 1:length(yl[1])]
    return ytl,t
end

"""
    glrka(f,y0::AbstractVector,ti::Number=0,tf::Number=10,n::Int=1000,kn::Int=10,fargs::AbstractVector=[])

# Arguments
- `f`: Function that describes dynamical system. 
- `y0::AbstractVector`: Initial conditions.
- `ti::Number`: Initial time.
- `tf::Number`: End time.
- `n::Int`: The number of timesteps.
- `s::Int`: Runge-Kutta stages.
- `fargs::AbstractVector`: Additional parameters to pass to f.


# Return 
The time evolution of supplied parameters and the corresponding time series.
Each index of the position series is a time evolution of one of the input parameters.
For example, if an initial condition is specified [y1 y2], then the output is [y1t y2t] 
where y1t = [yt1 yt2 ... ytn].

# Description
Gauss-Legendre Runge-Kutta Adaptive ODE solver[^list-of-runge-kutt-methods]. Good for nth order implicit and explicit 
Runge-Kutta, dynamical evolution of systems that can be cast as an array of ODEs. This method has an adaptive time step.

# References
[^list-of-runge-kutt-methods]: [List of Runge-Kutta Methods, https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods (accessed April 14, 2026).](https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods)

"""

function glrka(f,y0::AbstractVector,ti::Number=0,tf::Number=10,n::Int=1000,s::Int=10,fargs::AbstractVector=[],xn::Int=100000)
    y = y0
    yt = y
    
    yl = Vector{typeof(y0)}(undef,n)
    yl[1] = y0

    t = range(ti,tf,n)
    dt = t[2]-t[1]

    # calculating the butcher-tableau may need to be moved to its own function
    # set up polynomial in ascending order, index is polynomial power
    p1 = normalized_legendre_poly_coef(s) 
    p2 = normalized_legendre_poly_coef(s+1) 

    
    # find roots of polynomial
    roots1 = polyroots(p)
    roots2 = polyroots(p)

    aij1 = zeros(length(roots),length(roots))
    bi1 = zeros(length(roots))
    lj1 = zeros(xn)
    li1 = zeros(xn)

    aij2 = zeros(length(roots),length(roots))
    bi2 = zeros(length(roots))
    lj2 = zeros(xn)
    li2 = zeros(xn)

    for i in 1:length(roots1)
        for j in 1:length(roots1)
            lj1 = [lagrange_poly(roots1,x,j) for x in range(0,roots1[i],xn)]
            
            # try higher accuracy integrators to reduce iteration count
            aij1[i,j] = riemann_integration(range(0,roots1[i],xn),lj1)
        end
        li1 = [lagrange_poly(roots1,x,i) for x in range(0,1,xn)]
        bi1[i] = riemann_integration(range(0,1,xn),li1)
    end

    for i in 1:length(roots2)
        for j in 1:length(roots2)
            lj1 = [lagrange_poly(roots2,x,j) for x in range(0,roots2[i],xn)]
            
            # try higher accuracy integrators to reduce iteration count
            aij2[i,j] = riemann_integration(range(0,roots2[i],xn),lj2)
        end
        li2 = [lagrange_poly(roots2,x,i) for x in range(0,1,xn)]
        bi2[i] = riemann_integration(range(0,1,xn),li2)
    end

    # compute runge kutta
    # TODO: add adaptive check
    for nx in 2:n
        yi = [y + sum([aij[i,j]*f(y,t[nx],fargs) for j in 1:length(roots)])*dt for i in 1:length(roots)]
        y = y + sum([bi[i]*f(y,t[nx],fargs) for i in 1:length(roots)])*dt

        yl[nx] = y
    end

    ytl = [[row[i] for row in yl] for i in 1:length(yl[1])]
    return ytl,t
end

"""
    rk1(f,y0::AbstractVector,ti::Number=0,tf::Number=10,n::Int=1000,fargs::AbstractVector=[])


# Arguments
- `f`: Function that describes dynamical system. 
- `y0::AbstractVector`: Initial conditions.
- `ti::Number`: Initial time.
- `tf::Number`: End time. 
- `n::Int`: The number of timesteps.
- `fargs::AbstractVector`: Additional parameters to pass to f.

# Return 
The time evolution of supplied parameters and the corresponding time series.

# Description
First order Runge-Kutta ODE solver or identically Euler's method[^list-of-runge-kutt-methods], dynamical evolution 
of systems that can be cast as an array of ODEs.

# References
[^list-of-runge-kutt-methods]: [List of Runge-Kutta Methods, https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods (accessed April 14, 2026).](https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods)

"""
function rk1(f,y0::AbstractVector,ti::Number=0,tf::Number=10,n::Int=1000,fargs::AbstractVector=[])
    y = y0
    yt = y
    
    yl = Vector{typeof(y0)}(undef,n)
    yl[1] = y0

    t = range(ti,tf,n)
    dt = t[2]-t[1]

    k1 = typeof(y0)(undef,n)
    k2 = typeof(y0)(undef,n)

    for nx in 2:n
        y = y + f(yt,t[nx],fargs).*dt
        yl[nx] = y
    end

    ytl = [[row[i] for row in yl] for i in 1:length(yl[1])]
    return ytl,t
end


"""
    rk2(f,y0::AbstractVector,ti::Number=0,tf::Number=10,n::Int=1000,fargs::AbstractVector=[])


# Arguments
- `f`: Function that describes dynamical system. 
- `y0::AbstractVector`: Initial conditions.
- `ti::Number`: Initial time.
- `tf::Number`: End time. 
- `n::Int`: The number of timesteps.
- `fargs::AbstractVector`: Additional parameters to pass to f.

# Return 
The time evolution of supplied parameters and the corresponding time series.

# Description
Second order Runge-Kutta ODE solver[^list-of-runge-kutt-methods], dynamical evolution of systems that can be cast 
as an array of ODEs.

# References
[^list-of-runge-kutt-methods]: [List of Runge-Kutta Methods, https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods (accessed April 14, 2026).](https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods)

"""
function rk2(f,y0::AbstractVector,ti::Number=0,tf::Number=10,n::Int=1000,fargs::AbstractVector=[])

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
    rk3(f,y0::AbstractVector,ti::Number=0,tf::Number=10,n::Int=1000,fargs::AbstractVector=[])


# Arguments
- `f`: Function that describes dynamical system. 
- `y0::AbstractVector`: Initial conditions.
- `ti::Number`: Initial time.
- `tf::Number`: End time. 
- `n::Int`: The number of timesteps.
- `fargs::AbstractVector`: Additional parameters to pass to f.

# Return 
The time evolution of supplied parameters and the corresponding time series.

# Description
Third order Runge-Kutta ODE solver[^list-of-runge-kutt-methods], dynamical evolution of systems that can be cast 
as an array of ODEs.

# References
[^list-of-runge-kutt-methods]: [List of Runge-Kutta Methods, https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods (accessed April 14, 2026).](https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods)


"""
function rk3(f,y0::AbstractVector,ti::Number=0,tf::Number=10,n::Int=1000,fargs::AbstractVector=[])

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
    rk4(f,y0::AbstractVector,ti::Number=0,tf::Number=10,n::Int=1000,fargs::AbstractVector=[])


# Arguments
- `f`: Function that describes dynamical system. 
- `y0::AbstractVector`: Initial conditions.
- `ti::Number`: Initial time.
- `tf::Number`: End time. 
- `n::Int`: The number of timesteps.
- `fargs::AbstractVector`: Additional parameters to pass to f.

# Return 
The time evolution of supplied parameters and the corresponding time series.

# Description
Fourth order Runge-Kutta ODE solver[^list-of-runge-kutt-methods], dynamical evolution of systems that can be cast 
as an array of ODEs.

# References
[^list-of-runge-kutt-methods]: [List of Runge-Kutta Methods, https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods (accessed April 14, 2026).](https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods)

"""
function rk4(f,y0::AbstractVector,ti::Number=0,tf::Number=10,n::Int=1000,fargs::AbstractVector=[])

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

"""
    rk45(f,y0::AbstractVector,ti::Number=0,tf::Number=10,n::Int=1000,fargs::AbstractVector=[])


# Arguments
- `f`: Function that describes dynamical system. 
- `y0::AbstractVector`: Initial conditions.
- `ti::Number`: Initial time.
- `tf::Number`: End time. 
- `n::Int`: The number of timesteps.
- `fargs::AbstractVector`: Additional parameters to pass to f.

# Return 
The time evolution of supplied parameters and the corresponding time series. 
Each index of the position series is a time evolution of one of the input parameters.
For example, if an initial condition is specified [y1 y2], then the output is [y1t y2t] 
where y1t = [yt1 yt2 ... ytn].

# Description
Dormand-Price Runge-Kutta ODE solver[^list-of-runge-kutt-methods], dynamical evolution of systems that can be cast 
as an array of ODEs. This method has an adaptive time step.

# References
[^list-of-runge-kutt-methods]: [List of Runge-Kutta Methods, https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods (accessed April 14, 2026).](https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods)

"""
function rk45(f,y0::AbstractVector,ti::Number=0,tf::Number=10,n::Int=1000,fargs::AbstractVector=[])

    y = y0
    yt4 = y
    yt5 = y
    
    yl = Vector{typeof(y0)}(undef,n)
    yl[1] = y0

    t = range(ti,tf,n)
    dt = t[2]-t[1]

    k14 = typeof(y0)(undef,size(y0))
    k24 = typeof(y0)(undef,size(y0))
    k34 = typeof(y0)(undef,size(y0))
    k44 = typeof(y0)(undef,size(y0))


    k15 = typeof(y0)(undef,size(y0))
    k25 = typeof(y0)(undef,size(y0))
    k35 = typeof(y0)(undef,size(y0))
    k45 = typeof(y0)(undef,size(y0))

    # TODO 
    # - match dormand price coefficients
    # - add adaptive step
    for nx in 2:n
        k14 = f(y,t[nx],fargs).*dt

        yt4 = y + k1/2
        k24 = f(yt,t[nx]+dt/2,fargs).*dt

        yt4 = y + k2/2
        k34 = f(yt,t[nx]+dt/2,fargs).*dt

        yt4 = y + k3
        k44 = f(yt,t[nx]+dt,fargs).*dt

        k15 = f(y,t[nx],fargs).*dt

        yt5 = y + k1/2
        k25 = f(yt,t[nx]+dt/2,fargs).*dt

        yt5 = y + k2/2
        k35 = f(yt,t[nx]+dt/2,fargs).*dt

        yt5 = y + k3
        k45 = f(yt,t[nx]+dt,fargs).*dt

        

        y = y + (k1 + 2*k2 + 2*k3 + k4)/6

        yl[nx] = y
    end

    ytl = [[row[i] for row in yl] for i in 1:length(yl[1])]
    return ytl,t
end


