# using Rotations
# using Distributions
# using Base.Threads
using ProgressMeter
# using BiochemicalAlgorithms
# using Meshes
# using DataFrames

include("generate_record.jl")
include("quaternion_functions.jl")
include("bingham_functions_vol2.jl")
include("eval_hhb.jl")
include("eval.jl")

function refine2!(proteinA::Molecule{Float32}, proteinB::Molecule{Float32}, complexAB::Molecule{Float32}, results_docking::Tuple{DataFrame, Array{ComplexF32, 3}, Vector{Tuple{String, Vector3{Float32}}}, Array{Meshes.Point3f, 3}, Int32, Int32}, λ::Float32, runs::Int32, vdW::Bool)
    # load protein A,B and complex AB
    # proteinA = molecules(load_pdb(protein_A))[1]
    # proteinB = molecules(load_pdb(protein_B))[1]
    # complexAB = molecules(load_pdb(protein_complex))[1]
    
    # extract N best quaternions (current value: five)
    Q = Vector{QuaternionF32}()
    for i in eachrow(results_docking[1][1:5, :])
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
    # lock for threads (unsure how this really works)
    # lk = ReentrantLock()
    # for progress bar
    p = Progress(length(rotations))
    # min rmsd
    # rmsd_min = typemax(Float32)
    # calculate scorings for sampled rotations
    for i in eachindex(rotations)
        # generate rotation quaternion
        R = QuaternionF32(rotations[i][1], rotations[i][2], rotations[i][3], rotations[i][4])
        # generate record with new sampled rotation
        record = generate_record(results_docking[2], R, results_docking[3], results_docking[4], results_docking[5], results_docking[6], vdW)
        # calculate rmsd for record
        t = Vector3{Float32}(record.α, record.β, record.γ)
        R = (record.R[1], record.R[2], record.R[3])
        # --------------
        # |CHANGE HERE!|
        # --------------
        rmsd_eval = try eval(proteinA, proteinB, complexAB, R, t)
        catch
            typemax(Float32)
        end
        rec = (α=record.α, β=record.β, γ=record.γ, R=record.R, score=record.score, rmsd=rmsd_eval)
        # check if record is better than first one in
        # current results
        # lock(lk) do
        push!(results_docking[1], rec)
        next!(p)
    end
    # finish progress 
    finish!(p)
    sort!(results_docking[1], [:rmsd])
    # only keep best five 
    resize!(results_docking[1], 5)
    # return new scoring table
    # return results_docking
end