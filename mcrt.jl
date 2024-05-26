begin MCRT 

    struct State
        # position, since we are operating on a sphere only r is needed
        tau::AbstractArray{Float32, 1}
        # propagation direction, same principle as position
        mu::AbstractArray{Float32, 1}

    
end