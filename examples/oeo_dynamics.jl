include("../test/Simularity.jl")
#using .Simularity
using Plots

function Vx(y,t,args)

    o = args[1]
    A = args[2]
    w1 = args[3]
    T = args[4] 
    theta = args[5] 
    delta = args[6] 
    ea = args[7] 
    a = args[8] 
    em = args[9] 
    phi = args[10] 
    epc = args[11] 
    ep = args[12] 
    dt = args[13] 
    xl = args[14] 
    zl = args[15] 


    ztd = 0
    xtd = 0
    
    x = y[1]
    z = y[2]

    if t > T
        ztd = zl[-int(T/dt)]
        xtd = xl[-int(T/dt)]
    end

    dx = - o^2*x - delta*z - delta*a*ztd*sin(2*(xtd+theta)+phi)

    return [dx,dz]
end

function main()

    fargs[1] # o
    fargs[2] # A
    fargs[3] # w1
    fargs[4] # T
    fargs[5] # theta 
    fargs[6] # delta
    fargs[7] # ea
    fargs[8] # a
    fargs[9] # em 
    fargs[10] # phi
    fargs[11] # epc
    fargs[12] # ep
    fargs[13] # dt
    fargs[14] # xl
    fargs[15] # zl

    n = 1000
    x,t = Simularity.rk4(Vx,y,0,10,n,fargs)

    p = plot(t,x,
            label=["Position" "Velocity"],
            xlabel="Time (s)",
            ylabel="Amplitude"
            )
    
    display(p)
    readline() 
end

main()
