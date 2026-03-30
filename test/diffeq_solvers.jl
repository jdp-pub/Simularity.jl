function rk4(f,y0,ti,tf,n,fargs)

    y = y0
    yt = y
    
    yl = Vector{Vector{Float64}}(undef,n)
    yl[1] = y0

    tl = zeros(n)
    tl[1] = ti


    dt = (ti-tf)/n
    t = range(ti,tf,n)

    for nx in 2:n
        k1 = f(y,t[nx],fargs)*dt

        yt = y + k1/2
        k2 = f(y,t[nx],fargs)*dt

        yt = y + k2/2
        k3 = f(y,t[nx],fargs)*dt

        yt = y + k3
        k4 = f(y,t[nx],fargs)*dt

        y = y + (k1 + 2*k2 + 2*k3 + k4)/6

        yl[nx] = y
        tl[nx] = t[nx]
    end

    ytl = [[row[i] for row in yl] for i in 1:length(yl[1])]
    return ytl,tl
end
