using Distributions
# using SpecialFunctions

# create the lookup table for normalization constant F
function create_lookup_table_F()
    # initalize 3D matrix for values
    # might not be the best decision but works for now
    F = Array{Float32}(undef,100,100,100)

    # initalize Λ values
    λ1 = Int32(-100)
    λ2 = Int32(-99)
    λ3 = Int32(-98)
    
    # initalize result variable
    result = Float32(0)

    while λ3 < 0
        # The Quaternion Bingham Distribution, 3D Object Detection, 
        # and Dynamic Manipulation
        # set z values depending on λ1,2,3
        z1 = -λ1
        z2 = λ2 - λ1
        z3 = λ3 - λ1

        # series expansion of integral for calculating 
        # lookup table for F(Λ), Λ = diag(λ1,λ2,λ3)
        for i in 1:20, j in 1:20, k in 1:20
            result += (z1^(i-1) * z2^(j-1) * z3^(k-1) * i * j * k * Float32(0.5))/
                (factorial(i-1) * factorial(j-1) * factorial(k-1) * 2 * (i+j+k)) 
        end

        F[-λ1, -λ2, -λ3] = exp(λ1) * 2 * π^2 * result

        # Monte Carlo Pose Estimation with Quaternion Kernels and 
        # the Bingham Distribution
        # so far: this stuff does not work
        #=
        for i in 1:20, j in 1:20, k in 1:20 
            result += (gamma(Float32(i-1)+Float32(0.5))*gamma(Float32(j-1)+Float32(0.5))*gamma(Float32(k-1)+Float32(0.5))/
                gamma(Float32(i-1)+Float32(j-1)+Float32(k-1)+Float32(2))) * (λ1^(Float32(i-1))*λ2^(Float32(j-1))*λ3^(Float32(k-1))/
                factorial(i-1)*factorial(j-1)*factorial(k-1))
        end

        F[-λ1, -λ2, -λ3] = 2 * sqrt(π) * result
        =#

        λ1 += one(Int32)
        λ2 += one(Int32)
        λ3 += one(Int32)
    end

    return F
end

function estimate_V(Q::Vector{QuaternionF32})
    # create X matrix out of quaternions
    X = Matrix{Float32}(undef, 4, 0)
    for i in Q
        q = [i.s, i.v1, i.v2, i.v3]
        X = hcat(X, q)
    end
    V = eigvecs(Float32(1/length(Q))*X*Transpose(X))
    return V[:,[1,2,3]]
end

function estimate_Λ(Q::Vector{QuaternionF32}, V::Matrix{Float32}, F::Array{Float32, 3})
    # create X matrix out of quaternions
    X = Matrix{Float32}(undef, 4, 0)
    for i in Q
        q = [i.s, i.v1, i.v2, i.v3]
        X = hcat(X, q)
    end
    # initalize max_L and variables for λs
    max_L = -Inf32
    max_Λ = Matrix{Int32}(undef, 3, 3)
    # set λ1,2,3
    λ1 = Int32(-100)
    λ2 = Int32(-99)
    λ3 = Int32(-98)

    while λ3 < 0
        # generate Λ
        Λ = Int32[λ1 0 0; 0 λ2 0; 0 0 λ3]
        # maximize L function
        L = tr(Λ*Transpose(V)*X*Transpose(X)*V) - length(Q)*F[-λ1, -λ2, -λ3]
        if (L > max_L)
            println("yey")
            max_L = L
            max_Λ = Λ
        end
        # update λ values
        λ1 += one(Int32)
        λ2 += one(Int32)
        λ3 += one(Int32)
    end

    return max_Λ
end

function metropolis_hastings_sampler(start_q::QuaternionF32, Λ::Matrix{Int32}, V::Matrix{Float32}, S::Matrix{Float32}, F::Array{Float32, 3}, iter::Int32)
    # transforme quaternion into vector
    v = Float32[start_q.s, start_q.v1, start_q.v2, start_q.v3]
    # generate new vector via multidimensional gauss
    for i in 1:iter
        v_next = rand(MvNormal(v, S))
        # calculate acceptance rotation
        α = bingham(v_next,Λ,V,F)/bingham(v,Λ,V,F)
        #accept or reject
        u = rand(Uniform{Float32}(0,1))
        if (u <= α)
            v = v_next
        end
    end

    # transform vector back to unit quaternion
    q = QuaternionF32(v[1],v[2],v[3],v[4])
    q *= 1/norm(q)
    
    return q
end

function bingham(x::Vector{Float32}, Λ::Matrix{Int32}, V::Matrix{Float32}, F::Array{Float32, 3})
    # extract λs
    λ1 = Λ[1,1]
    λ2 = Λ[2,2]
    λ3 = Λ[3,3]
    
    return 1/F[-λ1, -λ2, -λ3] * exp(Transpose(x) * V * Λ * Transpose(V) * x)
end

function scatter_matrix(Q::Vector{QuaternionF32})
    # create X matrix out of quaternions
    X = Matrix{Float32}(undef, 4, 0)
    for i in Q
        q = [i.s, i.v1, i.v2, i.v3]
        X = hcat(X, q)
    end
    # return scatter matrix
    return Float32(1/length(Q)) * X * Transpose(X)
end