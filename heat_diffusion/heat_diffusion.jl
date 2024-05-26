using Plots

steps = 100
L = 50
dy = 1.0 / L
# dt = 10^(-4) 
dt = 0.01

function init_step(x::Float64)::Float64

    if x < 0.4 || x > 0.6
        return 0.
    else
        return 1.
    end

end

function init_gauss(x::Float64, sigma::Float64 = 0.1)::Float64

    return exp(-(x - 0.5)^2 / (sigma^2))

end

function f(y::Array{Float64, 1}, dx::Float64)::Array{Float64, 1}

    # add zero padding to the array
    y_pad = cat([0.], y, dims=1)
    y_pad = cat(y_pad, [0.], dims=1)

    return (y_pad[1:end-2] .- 2*y_pad[2:end-1] .+ y_pad[3:end]) ./ (dx*dx)

end

function euler_step(y ::Array{Float64, 1}, dx ::Float64, dt ::Float64)::Array{Float64, 1}

    return y .+ dt .* f(y, dx)

end

function runge_kutta_step(y::Array{Float64, 1}, dy::Float64, dt::Float64)::Array{Float64, 1}

    # advance half time step
    y_star = y .+ 0.5 .* dt .* f(y, dy)

    # advance another half step using the half step values
    return y .+ dt .* f(y_star, dy)

end

function boundary(Theta::Array{Float64, 1})::Array{Float64, 1}

    Theta[1] = 0.
    Theta[end] = 0.

    return Theta

end

x::Array{Float64, 1} = range(0, 1, length=L)

Theta_euler_step = zeros(Float64, (L, steps))
Theta_euler_step[:, 1] = init_step.(x)
Theta_euler_gauss = zeros(Float64, (L, steps))
Theta_euler_gauss[:, 1] = init_gauss.(x)

Theta_runge_step = zeros(Float64, (L, steps))
Theta_runge_step[:, 1] = init_step.(x)
Theta_runge_gauss = zeros(Float64, (L, steps))
Theta_runge_gauss[:, 1] = init_gauss.(x)

difference = zeros(Float64, (L, steps))

println("Running the simulation...")
for i in 2:steps
    # update the temperature fields
    Theta_euler_step[:, i] = euler_step(Theta_euler_step[:, i-1], dy, dt)
    Theta_euler_step[:, i] = boundary(Theta_euler_step[:, i])

    Theta_euler_gauss[:, i] = euler_step(Theta_euler_gauss[:, i-1], dy, dt)
    Theta_euler_gauss[:, i] = boundary(Theta_euler_gauss[:, i])

    Theta_runge_step[:, i] = runge_kutta_step(Theta_runge_step[:, i-1], dy, dt)
    Theta_runge_step[:, i] = boundary(Theta_runge_step[:, i])

    Theta_runge_gauss[:, i] = runge_kutta_step(Theta_runge_gauss[:, i-1], dy, dt)
    Theta_runge_gauss[:, i] = boundary(Theta_runge_gauss[:, i])

    difference[:, i] = sum(abs.(Theta_euler_gauss[:, i] .- Theta_runge_gauss[:, i]))
    println("Step: $i, diff: $diff")
end

println("Creating the gif...")
anim = @animate for i in 1:1:steps
    euler = plot(x, Theta_euler_step[:, i], label="Step init", title="Heat diffusion euler")
    euler = plot!(x, Theta_euler_gauss[:, i], label="Gauss init", title="Heat diffusion euler")
    ylims!(0, 1)
    runge = plot(x, Theta_runge_step[:, i], label="Step init", title="Heat diffusion 2nd order RK")
    runge = plot!(x, Theta_runge_gauss[:, i], label="Gauss init", title="Heat diffusion 2nd order RK")
    ylims!(0, 1)
    plot(euler, runge, layout=(1, 2), size=(800, 400))
end

gif(anim, "heat_diffusion.gif", fps = 30)

println("Done!")
