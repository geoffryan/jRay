include("ode.jl")

using ODE

f(x, y) = [y[2],-y[1]]

xa = 0.0
ya = [0.0, 1.0]
xb = 2
dx = 0.1

y_fe = ODE.evolve(xa, ya, xb, dx, f, ODE.fe)
y_rk2 = ODE.evolve(xa, ya, xb, dx, f, ODE.rk2_mp)
y_rk4 = ODE.evolve(xa, ya, xb, dx, f, ODE.rk4)

println(y_fe)
println(y_rk2)
println(y_rk4)
println([sin(xb), cos(xb)])
