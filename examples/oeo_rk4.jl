include("../test/Simularity.jl")
using Plots



function Vx(y,t,args)

    o = args[1]
    kx = args[2] 
    w1 = args[3]
    T = args[4] 
    theta = args[5] 
    delta = args[6] 
    a = args[7] 
    phi = args[8] 
    xl = args[9] 
    zl = args[10]
    gx = args[11]
    dt = args[12]

    x = y[1]
    z = y[2]

    if mod(kx,4)==0
        # keeps track of position in time
        xl[gx]=y[1]
        zl[gx]=y[2]
        args[11] = gx+1
    end

    # total function calls, grows 4x faster than gx for rk4
    args[2]=kx+1

    # time delay parameters
    ztd = 0
    xtd = 0
    

    # when current time is after the time delay reaches the mzm, get x,z at time t-T
    if gx > convert(Int,floor(T/dt))
        ztd = zl[gx-convert(Int,floor(T/dt))]
        xtd = xl[gx-convert(Int,floor(T/dt))]
    end

    # principle equations of motion
    dx = z
    dz = - o^2*x - delta*z - delta*a*ztd*sin(2*(xtd+theta+phi))

    return [dx,dz]
end

function main()

    # a more complicated example
    # 2nd order diff eq with time delay
    # Voltage of an optoelectronic oscillator in time

    n = 1000000

    fargs = Vector{Any}(undef,12)

    # a dictionary might work better here, ill have to experiment
    fargs[1] = 1.0*pi*2.0# o, bandpass filter central frequency
    fargs[2] = 1 # kx, total function calls
    fargs[3] = fargs[1] # w1, frequency examined
    fargs[4] = 200.5 # T, time delay
    fargs[5] = 0.0 # theta 
    fargs[6] = 5.0/1000.0*2.0*pi # delta, bandwidth of filter
    fargs[7] = 1.0 # a, round trip gain, <1 osccilates w/ gain, >1 decays, 1 is critical
    fargs[8] = pi/4.0 # phi
    tf = 50.0*fargs[4]
    fargs[9] = Vector{Float64}(undef,n) # xl
    fargs[10] = Vector{Float64}(undef,n) # zl
    fargs[11] = 1
    fargs[12] = tf/n # dt

    y = [0.1, 0.0]
    x,t = Simularity.rk4(Vx,y,0.0,tf,n,fargs)
    p = plot(t,x[1],
            label="OEO",
            xlabel="Time (s)",
            ylabel="Amplitude"
            )
    
    display(p)
    readline() 
end

main()
