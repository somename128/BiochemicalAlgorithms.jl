using JLD2
using DataFrames

#=
# grid 20 deg
hhb = load_object("src/dockings/testrun/2hhb_grid_20.jld2")
println(hhb)
mhb = load_object("src/dockings/testrun/2mhb_grid_20.jld2")
println(mhb)
ptc = load_object("src/dockings/testrun/2ptc_grid_20.jld2")
println(ptc)
apr = load_object("src/dockings/testrun/3apr_grid_20.jld2")
println(apr)
ts1 = load_object("src/dockings/testrun/3ts1_grid_20.jld2")
println(ts1)

# grid 120-cell
hhb = load_object("src/dockings/testrun/2hhb_grid_120.jld2")
println(hhb)
mhb = load_object("src/dockings/testrun/2mhb_grid_120.jld2")
println(mhb)
ptc = load_object("src/dockings/testrun/2ptc_grid_120.jld2")
println(ptc)
apr = load_object("src/dockings/testrun/3apr_grid_120.jld2")
println(apr)
ts1 = load_object("src/dockings/testrun/3ts1_grid_120.jld2")
println(ts1)
=#
#=
# vdW 20 deg
hhb = load_object("src/dockings/testrun/2hhb_vdW_20.jld2")
println(hhb)
mhb = load_object("src/dockings/testrun/2mhb_vdW_20.jld2")
println(mhb)
ptc = load_object("src/dockings/testrun/2ptc_vdW_20.jld2")
println(ptc)
apr = load_object("src/dockings/testrun/3apr_vdW_20.jld2")
println(apr)
ts1 = load_object("src/dockings/testrun/3ts1_vdW_20.jld2")
println(ts1)

# vdW 120-cell
hhb = load_object("src/dockings/testrun/2hhb_vdW_120.jld2")
println(hhb)
mhb = load_object("src/dockings/testrun/2mhb_vdW_120.jld2")
println(mhb)
ptc = load_object("src/dockings/testrun/2ptc_vdW_120.jld2")
println(ptc)
apr = load_object("src/dockings/testrun/3apr_vdW_120.jld2")
println(apr)
ts1 = load_object("src/dockings/testrun/3ts1_vdW_120.jld2")
println(ts1)
=#