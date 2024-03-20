using Plots

steps = 1000
L = 50
dx = 1.0 / L
dt = 0.0001

function init_step(x::Float64)::Float64

    if x < 0.4 || x > 0.6
        return 0.
    else
        return 1.
    end

end

function init_gauss(x::Float64, σ::Float64 = 0.1)::Float64

    return exp(-(x - 0.5)^2 / (σ^2))

end

function second_derivative(x::Array{Float64, 1}, dx::Float64)::Array{Float64, 1}

    x_ = cat([0.], x, dims=1)
    x_ = cat(x_, [0.], dims=1)

    return (x_[1:end-2] .- 2*x_[2:end-1] .+ x_[3:end]) ./ (dx*dx)

end

function euler_step(Theta ::Array{Float64, 1}, dx ::Float64, dt ::Float64)::Array{Float64, 1}

    return Theta .+ dt .* second_derivative(Theta, dx)

end

function boundary(Theta::Array{Float64, 1})::Array{Float64, 1}

    Theta[1] = 0.
    Theta[end] = 0.

    return Theta

end

x::Array{Float64, 1} = range(0, 1, length=L)

Theta_step = zeros(Float64, (L, steps))
Theta_step[:, 1] = init_step.(x)
Theta_gauss = zeros(Float64, (L, steps))
Theta_gauss[:, 1] = init_gauss.(x)

for i in 2:steps
    # update the temperature fields
    Theta_step[:, i] = euler_step(Theta_step[:, i-1], dx, dt)
    Theta_step[:, i] = boundary(Theta_step[:, i])

    Theta_gauss[:, i] = euler_step(Theta_gauss[:, i-1], dx, dt)
    Theta_gauss[:, i] = boundary(Theta_gauss[:, i])
end

p1  = plot(Theta_step[:, 1:20:end], label="", title="Step initial condition")
p2  = plot(Theta_gauss[:, 1:20:end], label="", title="Gaussian initial condition")

plot(p1, p2, layout=(2, 1), size=(800, 800))

gui()

sleep(50)