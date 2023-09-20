function reflector(A, n)
    A - 2*n*(transpose(n)*A)/(transpose(n)*n)
end

function rotationMatrix(u, v)
    S = reflector(I, u+v)
    R = reflector(S, v)

    R
end

# Based on algorithm VM* in [Wood, A.T. (1994). Simulation of the von mises fisher distribution.]
function sampleVMFMixture(λ, m)
    b  = (-2λ + sqrt(4λ^2 + (m-1)^2))/(m-1)
    x₀ = (1-b)/(1+b)
    c  = λ*x₀ + (m-1)*log(1-x₀^2)

    converged = false
    W = undef

    while !converged
        Z = rand(Beta((m-1)/2, (m-1)/2))
        U = rand(Uniform())
        W = (1 - (1+b)*Z)/(1-(1-b)*Z)

        if λ*W + (m-1)*log(1-x₀*W) -c >= log(U)
            converged = true
        end
    end

    W
end

function sampleVMF(μ::Vector{Float32}, λ::Int32)
    m = length(μ)

    W = sampleVMFMixture(λ, m)

    ξ = normalize(randn(m-1))

    x = vcat(sqrt(1-W^2)*ξ, W)

    # at this point, x has been sampled wrt to the modal vector (0, …, 0, 1)^t
    # now, we need to rotate it to be in line with modal vector μ
    u = zeros(m)
    u[end] = 1.0

    R = rotationMatrix(u, normalize(μ))

    R * x
end

function swapcols(M, i, j)
    k = M[:, i]

    M[:, i] = M[:, j]
    M[:, j] = k
end

#=
using AbstractPlotting

function plot_samples(μ,  s)
    scene = Scene()
    cam3d!(scene)

    wireframe!(scene, Sphere(Point3f0(0), 1f0), transparency=true, color=RGBA(0.0,0.0,0.4,0.4))
    mesh!(scene, Sphere(Point3f0(0), 1f0), transparency=true, color=RGBA(0.0,0.0,0.4,0.1))

    scatter!(scene, Point3f0(axis), color=:darkred, markersize=75.0)

    scatter!(scene, map(Point3f0, s), color=:darkblue, markersize=20)

    scene
end
=#