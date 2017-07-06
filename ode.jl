module ODE

function fe(x::AbstractFloat, y, dx::AbstractFloat, f, args...)
    y + dx*f(x, y, args...)
end

function rk2_mp(x::AbstractFloat, y, dx::AbstractFloat, f, args...)
    y1 = y + 0.5*dx*f(x, y, args...)
    y + dx*f(x+0.5*dx, y1, args...)
end

function rk4(x::AbstractFloat, y, dx::AbstractFloat, f, args...)
    k1 = f(x, y, args...)
    k2 = f(x+0.5*dx, y+0.5*dx*k1, args...)
    k3 = f(x+0.5*dx, y+0.5*dx*k2, args...)
    k4 = f(x+dx, y+dx*k3, args...)
    y + dx*(k1+2k2+2k3+k4)/6.
end

function evolve(xa, ya, xb, dx, f, scheme, args...; callback=nothing, maxiter=100000)
    x = xa
    y = ya
    last = false
    for i in 1:maxiter
        if x+dx > xb
            dx = xb-x
            last = true
        end
        y = scheme(x, y, dx, f, args...)
        x += dx
        if callback != nothing && callback(x, y)
            break
        end
        if last || x >= xb
            break
        end
    end

    y
end

end
