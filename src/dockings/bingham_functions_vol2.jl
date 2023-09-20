using NLsolve
using Distributions
using LinearAlgebra
using Rotations

# Sample from the Bingham distribution in n dimensions
# Uses the rejection method from
# Kent, J.T., Ganeiber, A.M., Mardia, K.V., 2013.
# A new method to simulate the Bingham and related distributions in
# directional data analysis with applications.

function f_star_bingham(x::Vector{Float32}, A::Matrix{Float32})
    exp(-transpose(x)*A*x)
end

function f_star_acg(x::Vector{Float32}, Ω::Matrix{Float32}, n::Int32)
    (transpose(x)*Ω*x)^(-n/2.0)
end

function sample(M::Matrix{Float32}, Z::Diagonal{Float32, Vector{Float32}}, nsamples::Int32=one(Float32))
    # determine the dimensionality
    n = Int32(size(Z)[1])

    # first, determine the value of b that minimizes the bound M(b)
    # to optimize the acceptance ratio
    b₀ = Float32(nlsolve(b -> sum(1.0./(b .+ 2.0.*Z)) - 1, [1.0]).zero[1])

    # Kent uses f^*_bing = exp(-x A x^t), with A = -M Z M^T
    # the minus sign is unconventional but simplifies notation
    A = -M * Z * transpose(M)

    Ω = Matrix{Float32}(I + 2.0 .* A./b₀)
    M_star = Float32(exp(-(n-b₀)/2.0)*((n/b₀)^(n/2.0)))

    # now, perform the rejection sampling
    samples = Vector{Vector{Float32}}()
    proposal_density = MvNormal(zeros(Float32, n), Hermitian(inv(Ω)))

    current_sample = zero(Int32)
    while current_sample < nsamples
        proposal = normalize(rand(proposal_density))
          
        W = rand(Float32)
        
        if W < f_star_bingham(proposal, A) / (M_star * f_star_acg(proposal, Ω, n))
            push!(samples, proposal)

            current_sample += one(Int32)
        end
    end

    samples
end

function sample_quaternions(μ::Vector{Float32}, λ::Int32, nsamples::Int32=one(Int32))
    # construct orthonormal matrix with normalized μ as first column
    axis = normalize(μ)

    M = hcat(axis, nullspace(axis * transpose(axis)))

    Z = Diagonal{Float32}(vcat([0.0], repeat([λ], length(μ)-1)))
  
    sample(M, Z, nsamples)
end

function sample_xyz(X, Y, Z, λ, nsamples=1)
    # X,y,z, [0,2π]
    r_euler = RotXYZ(X, Y, Z)
    μ = Rotations.params(QuatRotation(r_euler))

    μ_new = sample_quaternions(μ, λ, nsamples)

    r_euler_new = [RotXYZ(QuatRotation(q)) for q in μ_new]

    convert_angles = s -> getproperty.(r_euler_new, Ref(s))

    return convert_angles(:theta1), convert_angles(:theta2), convert_angles(:theta3)
end