using BiochemicalAlgorithms: _SystemNucleotideTuple, _nucleotides

@testset "Nucleotide" begin
    for T in [Float32, Float64]
        sys = System{T}()
        mol = Molecule(sys)
        mol2 = Molecule(sys)
        chain = Chain(mol)
        chain2 = Chain(mol2)

        # constructors + parent
        nuc = Nucleotide(chain, 1)
        @test nuc isa Nucleotide{T}
        @test parent(nuc) === sys
        @test parent_system(nuc) === sys

        if T == Float32
            mol_ds = Molecule()
            chain_ds = Chain(mol_ds)
            nuc_ds = Nucleotide(chain_ds, 1)
            parent(nuc_ds) === default_system()
            parent_system(nuc_ds) === default_system()

            Nucleotide(chain_ds, 1, "something", Properties(:a => "b"), Flags([:A]))
        end

        nuc2 = Nucleotide(chain2, 1, "something", Properties(:a => 1), Flags([:A, :B]))

        #=
            Make sure we test for the correct number of fields.
            Add missing tests if the following test fails!
        =#
        @test length(getfield(nuc, :_row)) == 7

        # getproperty
        @test nuc.idx isa Int
        @test nuc.number isa Int
        @test nuc.number == 1
        @test nuc.name isa String
        @test nuc.name == ""
        @test nuc.properties isa Properties
        @test nuc.properties == Properties()
        @test nuc.flags isa Flags
        @test nuc.flags == Flags()

        @test nuc._sys isa System{T}
        @test nuc._row isa DataFrameRow
        
        @test_throws ErrorException nuc.molecule_id
        @test_throws ErrorException nuc.chain_id

        @test nuc._row.molecule_id isa Int
        @test nuc._row.molecule_id == mol.idx
        @test nuc._row.chain_id isa Int
        @test nuc._row.chain_id == chain.idx

        @test nuc2.number == 1
        @test nuc2.name == "something"
        @test nuc2.properties == Properties(:a => 1)
        @test nuc2.flags == Flags([:A, :B])
        @test nuc2._row.molecule_id == mol2.idx
        @test nuc2._row.chain_id == chain2.idx

        # setproperty!
        nuc.number = 0
        @test nuc.number == 0
        nuc.name = "something else"
        @test nuc.name == "something else"
        nuc.properties = Properties(:first => "v1", :second => 99)
        @test length(nuc.properties) == 2
        @test nuc.properties[:first] == "v1"
        @test nuc.properties[:second] == 99
        nuc.flags = Flags([:C])
        @test length(nuc.flags) == 1
        @test :C in nuc.flags

        @test_throws ErrorException nuc.molecule_id = 0
        @test_throws ErrorException nuc.chain_id = 0

        # nucleotide_by_idx
        @test isnothing(nucleotide_by_idx(sys, -1))
        @test nucleotide_by_idx(sys, nuc.idx) isa Nucleotide{T}
        @test nucleotide_by_idx(sys, nuc.idx) == nuc

        # _nucleotides
        df = _nucleotides(sys)
        @test df isa AbstractDataFrame
        @test size(df) == (2, length(fieldnames(_SystemNucleotideTuple)))
        @test copy(df[1, 1:length(fieldnames(NucleotideTuple))]) isa NucleotideTuple
        @test size(_nucleotides(sys), 1) == 2
        @test size(_nucleotides(sys, molecule_id = -1), 1) == 0
        @test size(_nucleotides(sys, molecule_id = mol.idx), 1) == 1
        @test size(_nucleotides(sys, molecule_id = mol2.idx), 1) == 1
        @test size(_nucleotides(sys, molecule_id = nothing), 1) == 2
        @test size(_nucleotides(sys, chain_id = -1), 1) == 0
        @test size(_nucleotides(sys, chain_id = chain.idx), 1) == 1
        @test size(_nucleotides(sys, chain_id = chain2.idx), 1) == 1
        @test size(_nucleotides(sys, chain_id = nothing), 1) == 2
        @test size(_nucleotides(sys, molecule_id = -1, chain_id = chain.idx), 1) == 0
        @test size(_nucleotides(sys, molecule_id = mol.idx, chain_id = -1), 1) == 0
        @test size(_nucleotides(sys, molecule_id = mol.idx, chain_id = chain.idx), 1) == 1
        @test size(_nucleotides(sys, molecule_id = mol.idx, chain_id = nothing), 1) == 1
        @test size(_nucleotides(sys, molecule_id = nothing, chain_id = chain.idx), 1) == 1
        @test size(_nucleotides(sys, molecule_id = nothing, chain_id = nothing), 1) == 2

        # nucleotides_df
        df = nucleotides_df(sys)
        @test df isa AbstractDataFrame
        @test size(df) == (2, length(fieldnames(NucleotideTuple)))
        @test copy(df[1, :]) isa NucleotideTuple
        @test size(nucleotides_df(sys), 1) == 2
        @test size(nucleotides_df(sys, molecule_id = -1), 1) == 0
        @test size(nucleotides_df(sys, molecule_id = mol.idx), 1) == 1
        @test size(nucleotides_df(sys, molecule_id = mol2.idx), 1) == 1
        @test size(nucleotides_df(sys, molecule_id = nothing), 1) == 2
        @test size(nucleotides_df(sys, chain_id = -1), 1) == 0
        @test size(nucleotides_df(sys, chain_id = chain.idx), 1) == 1
        @test size(nucleotides_df(sys, chain_id = chain2.idx), 1) == 1
        @test size(nucleotides_df(sys, chain_id = nothing), 1) == 2
        @test size(nucleotides_df(sys, molecule_id = -1, chain_id = chain.idx), 1) == 0
        @test size(nucleotides_df(sys, molecule_id = mol.idx, chain_id = -1), 1) == 0
        @test size(nucleotides_df(sys, molecule_id = mol.idx, chain_id = chain.idx), 1) == 1
        @test size(nucleotides_df(sys, molecule_id = mol.idx, chain_id = nothing), 1) == 1
        @test size(nucleotides_df(sys, molecule_id = nothing, chain_id = chain.idx), 1) == 1
        @test size(nucleotides_df(sys, molecule_id = nothing, chain_id = nothing), 1) == 2

        # nucleotides
        fv = nucleotides(sys)
        @test fv isa Vector{Nucleotide{T}}
        @test length(fv) == 2
        @test length(nucleotides(sys)) == 2
        @test length(nucleotides(sys, molecule_id = -1)) == 0
        @test length(nucleotides(sys, molecule_id = mol.idx)) == 1
        @test length(nucleotides(sys, molecule_id = mol2.idx)) == 1
        @test length(nucleotides(sys, molecule_id = nothing)) == 2
        @test length(nucleotides(sys, chain_id = -1)) == 0
        @test length(nucleotides(sys, chain_id = chain.idx)) == 1
        @test length(nucleotides(sys, chain_id = chain2.idx)) == 1
        @test length(nucleotides(sys, chain_id = nothing)) == 2
        @test length(nucleotides(sys, molecule_id = -1, chain_id = chain.idx)) == 0
        @test length(nucleotides(sys, molecule_id = mol.idx, chain_id = -1)) == 0
        @test length(nucleotides(sys, molecule_id = mol.idx, chain_id = chain.idx)) == 1
        @test length(nucleotides(sys, molecule_id = mol.idx, chain_id = nothing)) == 1
        @test length(nucleotides(sys, molecule_id = nothing, chain_id = chain.idx)) == 1
        @test length(nucleotides(sys, molecule_id = nothing, chain_id = nothing)) == 2

        # eachnucleotide
        @test length(eachnucleotide(sys)) == 2
        @test length(eachnucleotide(sys, molecule_id = -1)) == 0
        @test length(eachnucleotide(sys, molecule_id = mol.idx)) == 1
        @test length(eachnucleotide(sys, molecule_id = mol2.idx)) == 1
        @test length(eachnucleotide(sys, molecule_id = nothing)) == 2
        @test length(eachnucleotide(sys, chain_id = -1)) == 0
        @test length(eachnucleotide(sys, chain_id = chain.idx)) == 1
        @test length(eachnucleotide(sys, chain_id = chain2.idx)) == 1
        @test length(eachnucleotide(sys, chain_id = nothing)) == 2
        @test length(eachnucleotide(sys, molecule_id = -1, chain_id = chain.idx)) == 0
        @test length(eachnucleotide(sys, molecule_id = mol.idx, chain_id = -1)) == 0
        @test length(eachnucleotide(sys, molecule_id = mol.idx, chain_id = chain.idx)) == 1
        @test length(eachnucleotide(sys, molecule_id = mol.idx, chain_id = nothing)) == 1
        @test length(eachnucleotide(sys, molecule_id = nothing, chain_id = chain.idx)) == 1
        @test length(eachnucleotide(sys, molecule_id = nothing, chain_id = nothing)) == 2

        # nnucleotides
        @test nnucleotides(sys) isa Int
        @test nnucleotides(sys) == 2
        @test nnucleotides(sys, molecule_id = -1) == 0
        @test nnucleotides(sys, molecule_id = mol.idx) == 1
        @test nnucleotides(sys, molecule_id = mol2.idx) == 1
        @test nnucleotides(sys, molecule_id = nothing) == 2
        @test nnucleotides(sys, chain_id = -1) == 0
        @test nnucleotides(sys, chain_id = chain.idx) == 1
        @test nnucleotides(sys, chain_id = chain2.idx) == 1
        @test nnucleotides(sys, chain_id = nothing) == 2
        @test nnucleotides(sys, molecule_id = -1, chain_id = chain.idx) == 0
        @test nnucleotides(sys, molecule_id = mol.idx, chain_id = -1) == 0
        @test nnucleotides(sys, molecule_id = mol.idx, chain_id = chain.idx) == 1
        @test nnucleotides(sys, molecule_id = mol.idx, chain_id = nothing) == 1
        @test nnucleotides(sys, molecule_id = nothing, chain_id = chain.idx) == 1
        @test nnucleotides(sys, molecule_id = nothing, chain_id = nothing) == 2

        # chain/molecule nucleotides
        mol3 = Molecule(sys)
        @test size(_nucleotides(mol3), 1) == 0
        @test _nucleotides(mol3) == _nucleotides(sys, molecule_id = mol3.idx)
        @test size(nucleotides_df(mol3), 1) == 0
        @test nucleotides_df(mol3) == nucleotides_df(sys, molecule_id = mol3.idx)
        @test size(nucleotides(mol3), 1) == 0
        @test nucleotides(mol3) == nucleotides(sys, molecule_id = mol3.idx)
        @test length(eachnucleotide(mol3)) == 0
        @test length(eachnucleotide(mol3)) == length(eachnucleotide(sys, molecule_id = mol3.idx))
        @test nnucleotides(mol3) == 0
        @test nnucleotides(mol3) == nnucleotides(sys, molecule_id = mol3.idx)

        chain3 = Chain(mol3)
        @test size(_nucleotides(chain3), 1) == 0
        @test _nucleotides(chain3) == _nucleotides(sys, chain_id = chain3.idx)
        @test size(nucleotides_df(chain3), 1) == 0
        @test nucleotides_df(chain3) == nucleotides_df(sys, chain_id = chain3.idx)
        @test size(nucleotides(chain3), 1) == 0
        @test nucleotides(chain3) == nucleotides(sys, chain_id = chain3.idx)
        @test length(eachnucleotide(chain3)) == 0
        @test length(eachnucleotide(chain3)) == length(eachnucleotide(sys, chain_id = chain3.idx))
        @test nnucleotides(chain3) == 0
        @test nnucleotides(chain3) == nnucleotides(sys, chain_id = chain3.idx)

        Nucleotide(chain3, 1)
        @test size(_nucleotides(mol3), 1) == 1
        @test size(_nucleotides(mol3)) == size(_nucleotides(sys, molecule_id = mol3.idx))
        @test size(nucleotides_df(mol3), 1) == 1
        @test size(nucleotides_df(mol3)) == size(nucleotides_df(sys, molecule_id = mol3.idx))
        @test size(nucleotides(mol3), 1) == 1
        @test nucleotides(mol3) == nucleotides(sys, molecule_id = mol3.idx)
        @test length(eachnucleotide(mol3)) == 1
        @test length(eachnucleotide(mol3)) == length(eachnucleotide(sys, molecule_id = mol3.idx))
        @test nnucleotides(mol3) == 1
        @test nnucleotides(mol3) == nnucleotides(sys, molecule_id = mol3.idx)

        @test size(_nucleotides(chain3), 1) == 1
        @test size(_nucleotides(chain3)) == size(_nucleotides(sys, chain_id = chain3.idx))
        @test size(nucleotides_df(chain3), 1) == 1
        @test size(nucleotides_df(chain3)) == size(nucleotides_df(sys, chain_id = chain3.idx))
        @test size(nucleotides(chain3), 1) == 1
        @test nucleotides(chain3) == nucleotides(sys, chain_id = chain3.idx)
        @test length(eachnucleotide(chain3)) == 1
        @test length(eachnucleotide(chain3)) == length(eachnucleotide(sys, chain_id = chain3.idx))
        @test nnucleotides(chain3) == 1
        @test nnucleotides(chain3) == nnucleotides(sys, chain_id = chain3.idx)

        # nucleotide atoms
        @test size(_atoms(nuc), 1) == 0
        @test _atoms(nuc) == _atoms(sys, nucleotide_id = nuc.idx)
        @test size(atoms_df(nuc), 1) == 0
        @test atoms_df(nuc) == atoms_df(sys, nucleotide_id = nuc.idx)
        @test length(atoms(nuc)) == 0
        @test atoms(nuc) == atoms(sys, nucleotide_id = nuc.idx)
        @test length(eachatom(nuc)) == 0
        @test length(eachatom(nuc)) == length(eachatom(sys, nucleotide_id = nuc.idx))
        @test natoms(nuc) == 0
        @test natoms(nuc) == natoms(sys, nucleotide_id = nuc.idx)

        @test push!(nuc, AtomTuple{T}(1, Elements.H)) === nuc
        @test size(_atoms(nuc), 1) == 1
        @test size(_atoms(nuc)) == size(_atoms(sys, nucleotide_id = nuc.idx))
        @test size(atoms_df(nuc), 1) == 1
        @test size(atoms_df(nuc)) == size(atoms_df(sys, nucleotide_id = nuc.idx))
        @test length(atoms(nuc)) == 1
        @test atoms(nuc) == atoms(sys, nucleotide_id = nuc.idx)
        @test length(eachatom(nuc)) == 1
        @test length(eachatom(nuc)) == length(eachatom(sys, nucleotide_id = nuc.idx))
        @test natoms(nuc) == 1
        @test natoms(nuc) == natoms(sys, nucleotide_id = nuc.idx)

        # nucleotide bonds
        @test size(_bonds(nuc), 1) == 0
        @test _bonds(nuc) == _bonds(sys, nucleotide_id = nuc.idx)
        @test size(bonds_df(nuc), 1) == 0
        @test bonds_df(nuc) == bonds_df(sys, nucleotide_id = nuc.idx)
        @test length(bonds(nuc)) == 0
        @test bonds(nuc) == bonds(sys, nucleotide_id = nuc.idx)
        @test length(eachbond(nuc)) == 0
        @test length(eachbond(nuc)) == length(eachbond(sys, nucleotide_id = nuc.idx))
        @test nbonds(nuc) == 0
        @test nbonds(nuc) == nbonds(sys, nucleotide_id = nuc.idx)

        @test push!(nuc, BondTuple(
            Atom(nuc, 1, Elements.H).idx,
            Atom(nuc, 2, Elements.C).idx,
            BondOrder.Single
        )) === nuc
        @test size(_bonds(nuc), 1) == 1
        @test size(_bonds(nuc)) == size(_bonds(sys, nucleotide_id = nuc.idx))
        @test size(bonds_df(nuc), 1) == 1
        @test size(bonds_df(nuc)) == size(bonds_df(sys, nucleotide_id = nuc.idx))
        @test length(bonds(nuc)) == 1
        @test bonds(nuc) == bonds(sys, nucleotide_id = nuc.idx)
        @test length(eachbond(nuc)) == 1
        @test length(eachbond(nuc)) == length(eachbond(sys, nucleotide_id = nuc.idx))
        @test nbonds(nuc) == 1
        @test nbonds(nuc) == nbonds(sys, nucleotide_id = nuc.idx)
    end
end
