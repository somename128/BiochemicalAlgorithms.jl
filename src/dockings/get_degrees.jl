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

    return α, β, γ

end