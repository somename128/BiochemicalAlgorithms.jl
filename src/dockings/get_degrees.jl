include("helpers.jl")

function get_degrees(rotation::Matrix3{Float32})

    if (rotation[1,3] == 1 || rotation[1,3] == -1)
        β = rotation[1,3]*90
        γ = rad2deg(atan(rotation[2,1],rotation[2,2]))
        α = 0   
    else    
        β = rad2deg(asin(rotation[1,3]))
        γ = rad2deg(atan(-rotation[1,2],rotation[1,1]))
        α = rad2deg(atan(-rotation[2,3],rotation[3,3])) 
    end

    # set almost zero values to zero
    α = zero_small_one_value(α, Float32(0.5))
    β = zero_small_one_value(β, Float32(0.5))
    γ = zero_small_one_value(γ, Float32(0.5))

    return α, β, γ

end