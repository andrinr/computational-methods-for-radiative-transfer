using Plots

println("Hello, world!")
x = range(0, 1, length=100)

function step(x)
    if x < 0.4
        return 0
    else if x < 0.6
        return 1
    else
        return 0
    end
end

function gauss(x, sigma = 0.1)
    return 1 + exp(-(x - 0.5)^2 / (sigma^2))

y_step = step.(x)

println(y_step)

y_gauss = gauss.(x)

display(plot(x, y_step, label="step"))
display(plot(x, y_gauss, label="gauss"))

# make sure plot is displayed
sleep(10)
