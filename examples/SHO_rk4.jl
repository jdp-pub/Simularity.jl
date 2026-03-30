include("../test/Simularity.jl")
#using .Simularity
using Plots

# function pattern to pass to rk4
function SHO(y,t,args)

    # simple harmonic oscillator model
    k = args[1]
    m = args[2]

    # cast second order diffeq as 2 first order diffeqs
    x = y[1]
    z = y[2]

    # set function
    dx = z
    dz = -k/m*x
    
    return [dx,dz]
end


function main()

    # number of steps
    n = 1000
    
    # initial conditions
    y = [1,0]
    
    # any arguments to pass to the function
    fargs = [1,1]
    
    # compute with rk4
    x,t = Simularity.rk4(SHO,y,0,10,n,fargs)

    p = plot(t,x,
            label=["Position" "Velocity"],
            xlabel="Time (s)",
            ylabel="Amplitude"
            )
    
    display(p)
    readline()
end

main()
