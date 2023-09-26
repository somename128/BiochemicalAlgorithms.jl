using Rotations
using Distributions
using Base.Threads
using ProgressMeter
using BiochemicalAlgorithms
using Meshes
using DataFrames

include("generate_record.jl")
include("quaternion_functions.jl")
include("bingham_functions_vol2.jl")
include("eval_hhb.jl")
include("eval.jl")

function refine2!(results_docking::Tuple{DataFrame, Array{ComplexF32, 3}, Vector{Tuple{String, Vector3{Float32}}}, Array{Meshes.Point3f, 3}, Int32, Int32}, λ::Float32, runs::Int32, vdW::Bool)

    # store scoring table, grid of protein A, roomcoordinates of protein b,
    # centroids and gridsize
    scoring_table = results_docking[1]
    # insertcols!(scoring_table, :rmsd => Float32[0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    grid_A = results_docking[2]
    roomcoordiantes_B = results_docking[3]
    # println(typeof(roomcoordiantes_B[1][2]))
    centroids = results_docking[4]
    gridsize = results_docking[5]
    resolution = results_docking[6]

    proteinA = "src/dockings/testproteins/2hhb_alpha_chain.pdb"
    proteinB = "src/dockings/testproteins/2hhb_beta_chain.pdb"
    complexAB = "src/dockings/testproteins/2hhb.pdb"

    # extract N best quaternions (current value: five)
    Q = Vector{QuaternionF32}()
    for i in eachrow(scoring_table[1:5, :])
        q = extract_quaternion(i)
        push!(Q, q)
    end

    rotations = Vector{Vector{Float32}}()
    for i in eachindex(Q)
    # set parameters for quaternion sampling
        μ = Float32[Q[i].s, Q[i].v1, Q[i].v2, Q[i].v3]
        # sample quaternions Bingham distributed around best quaternion
        append!(rotations, sample_quaternions(μ, λ, runs))
    end
    println("got rotations...")
    # lock for threads (unsure how this really works)
    # lk = ReentrantLock()
    # for progress bar
    # p = Progress(length(rotations))
    # min rmsd
    # rmsd_min = typemax(Float32)
    # calculate scorings for sampled rotations
    println("Start loop...")
    for i in eachindex(rotations)
        # generate rotation quaternion
        R = QuaternionF32(rotations[i][1], rotations[i][2], rotations[i][3], rotations[i][4])
        # generate record with new sampled rotation
        println("generate record...")
        record = generate_record(grid_A, R, roomcoordiantes_B, centroids, gridsize, resolution, vdW)
        # calculate rmsd for record
        t = Vector3{Float32}(record.α, record.β, record.γ)
        R = (record.R[1], record.R[2], record.R[3])
        println("calculate rmsd...")
        rmsd_eval = try eval_hhb(proteinA, proteinB, complexAB, R, t)
        catch
            typemax(Float32)
        end
        record = (α=record.α, β=record.β, γ=record.γ, R=record.R, score=record.score, rmsd=rmsd_eval)
        # check if record is better than first one in
        # current results
        # lock(lk) do
        println("push record...")
        push!(scoring_table, record)
        # end
        # next!(p)
    end
    # finish progress 
    # finish!(p)
    sort!(scoring_table, [:rmsd])
    # only keep best five 
    resize!(scoring_table, 5)
    # return new scoring table
    return results_docking
end