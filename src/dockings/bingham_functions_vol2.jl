using NLsolve
using Distributions
using LinearAlgebra
using Rotations

# Sample from the Bingham distribution in n dimensions
# Uses the rejection method from
# Kent, J.T., Ganeiber, A.M., Mardia, K.V., 2013.
# A new method to simulate the Bingham and related distributions in
# directional data analysis with applications.

function f_star_bingham(x, A)
    exp(-transpose(x)*A*x)
end

function f_star_acg(x, Ω, n)
    (transpose(x)*Ω*x)^(-n/2.0)
end

function sample(M, Z, nsamples=1)
    # determine the dimensionality
    n = size(Z)[1]

    # first, determine the value of b that minimizes the bound M(b)
    # to optimize the acceptance ratio
    b₀ = nlsolve(b -> sum(1.0./(b .+ 2.0.*Z)) - 1, [1.0]).zero[1]

    # Kent uses f^*_bing = exp(-x A x^t), with A = -M Z M^T
    # the minus sign is unconventional but simplifies notation
    A = -M * Z * transpose(M)

    Ω = I + 2.0 .* A./b₀
    M_star = exp(-(n-b₀)/2.0)*((n/b₀)^(n/2.0))

    # now, perform the rejection sampling
    samples = []

    proposal_density = MvNormal(zeros(n), Hermitian(inv(Ω)))

    current_sample = 0
    while current_sample < nsamples
        proposal = normalize(rand(proposal_density))

        W = rand()

        if W < f_star_bingham(proposal, A) / (M_star * f_star_acg(proposal, Ω, n))
            push!(samples, proposal)

            current_sample += 1
        end
    end

    samples
end

function sample_quaternions(μ, λ, nsamples=1)
    # construct orthonormal matrix with normalized μ as first column
    axis = normalize(μ)

    M = hcat(axis, nullspace(axis * transpose(axis)))

    Z = Diagonal(vcat([0.0], repeat([λ], length(μ)-1)))

    sample(M, Z, nsamples)
end
