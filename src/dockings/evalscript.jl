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
insertcols!(sc[1], :rmsd => Float32[typemax(Float32), typemax(Float32), typemax(Float32), typemax(Float32), typemax(Float32), typemax(Float32), typemax(Float32), typemax(Float32), typemax(Float32), typemax(Float32)])
# print initalization result
println(sc[1])

# set rmsd to infinity
# rmsd_min = typemax(Float32)
# standard deviation
λ = Float32(100)
# counter = 1
@time while Float32(2.5) < sc[1][1, :].rmsd
    # refine
    @time global sc = refine2!(sc, λ, Int32(5), vdW)

    #=
    # eval
    for i in 1:10
        t = Vector3{Float32}(sc[1][i, :].α, sc[1][i, :].β, sc[1][i, :].γ)
        R = (sc[1][i, :].R[1], sc[1][i, :].R[2], sc[1][i, :].R[3])
        rmsd = eval_hhb(proteinA, proteinB, complexAB, R, t)
        # println(rmsd)
        
        if (rmsd_min >= rmsd)
            global rmsd_min = rmsd
        end

    end
    =#
    println(sc[1])
    println("Aktuelles Minimum: ", sc[1][1, :].rmsd)
    println(λ)

    # every second time half the standard deviation
    # global  counter += 1
    # if (counter % 2 == 0)
    global λ *= Float32(0.75)
    # end
    # else
    #     global λ += λ - 10
    # end
    
end


