using Meshes

function create_atomballs(atoms::Vector{Tuple{String, Vector3{Float32}}}, thickness::Float32)
    # println("Extract room coordinates...")

    # convert SVector to Meshes.Point3f in Tuple for creating balls
    atoms = [(i[1], convert(Meshes.Point3f, i[2])) for i in atoms]

    # balls with radius r and atom points as center
    # println("Build atomballs...")
    atomballs = Array{Meshes.Ball{3,Float32}}(undef, length(atoms))
    # dictionary for radii of elements
    radii = Dict("C" => 1.7, "H" => 1.0, "N" => 1.5, "O" => 1.4)
    # i is the tuple of (elementname, coordinates)
    # first set radii depending on dictionary
    # second create atomball from coordinate and radius
    # default radius 1.8 Å
    for i in eachindex(atoms)
        # for simple geometry
        # radius = Float32(1)
        # for huge cube
        # radius = Float32(0.75)
        # for testproteins default = 1.8
        radius = get(radii, atoms[i][1], Float32(1.8))
        atomballs[i] = Meshes.Ball{3,Float32}(atoms[i][2], radius+thickness)
    end

    return atomballs
end