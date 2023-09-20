using BiochemicalAlgorithms
using JLD2
using Meshes

include("eval.jl")
include("eval_hhb.jl")
include("refine2!.jl")
include("refine3!.jl")

# vdW or grid surface
vdW = true

# load result after initialization and protein a and b
sc = load_object("src/dockings/testrun_huge/2hhb_vdW_120.jld2")
proteinA = "src/dockings/testproteins/2hhb_alpha_chain.pdb"
proteinB = "src/dockings/testproteins/2hhb_beta_chain.pdb"
complexAB = "src/dockings/testproteins/2hhb.pdb"
# print initalization result
println(sc[1])
# store results
results_rmsd = Float32[]

# set rmsd to infinity
rmsd = typemax(Float32)

@time while 0 < rmsd
    # refine
    @time score_refined = refine2!(sc, Int32(100), vdW)
    println(score_refined[1][1:10, :])

    # eval
    t = Vector3{Float32}(sc[1][1, :].α, sc[1][1, :].β, sc[1][1, :].β)
    R = (sc[1][1, :].R[1], sc[1][1, :].R[2], sc[1][1, :].R[3])
    global rmsd = eval_hhb(proteinA, proteinB, complexAB, R, t)
    # store and print rmsd
    push!(results_rmsd, rmsd)
    println(rmsd)
end

save_object("src/dockings/results_rmsd.jld2", results_rmsd)
results_rmsd
