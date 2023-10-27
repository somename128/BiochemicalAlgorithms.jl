# using BiochemicalAlgorithms
using JLD2
# using Meshes

# include("eval.jl")
# include("eval_hhb.jl")
include("refine2!.jl")
# include("refine3!.jl")
function evalscript(proteinA::Molecule{Float32}, proteinB::Molecule{Float32}, complexAB::Molecule{Float32}, sc::Tuple{DataFrame, Array{ComplexF32, 3}, Vector{Tuple{String, Vector3{Float32}}}, Array{Meshes.Point3f, 3}, Int32, Int32}, λ::Float32, vdW::Bool)
    
    rmsds = Vector{Float32}()
    
    # print initalization result
    println(sc[1])

    # set rmsd to infinity
    # rmsd_min = typemax(Float32)
    # standard deviation
    # counter = 1
    # while sc[1][1, :].rmsd > 0.5
    for i in 1:100
        println("------------------")
        println("#Iterations: ", i)
        # refine
        refine2!(proteinA, proteinB, complexAB, sc, λ, Int32(5), vdW)

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
        println("Current minimum: ", sc[1][1, :].rmsd)
        println("λ: ", λ)
        # push current min in vector 
        push!(rmsds, sc[1][1, :].rmsd)

        # every second time half the standard deviation
        # global  counter += 1
        if (λ < 0.3)
            λ = Float32(5) * sc[1][1, :].rmsd
        else
            λ *= Float32(0.75)
        end
    end

    println("Final result")
    println(sc[1])

    # write result as latex table format in file
    # --------------
    # |CHANGE HERE!|
    # --------------
    #=
    open("src/dockings/latex/3ts1/refine_vdW_120.txt", "a") do io
        show(io, MIME("text/latex"), sc[1])
    end
    # add new line
    # --------------
    # |CHANGE HERE!|
    # --------------
    write(open("src/dockings/latex/3ts1/refine_vdW_120.txt", "a") , "\n")
    =#
    return rmsds
end