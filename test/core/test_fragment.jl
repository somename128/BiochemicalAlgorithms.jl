using BiochemicalAlgorithms: _SystemFragmentTuple, _fragments

@testset "Fragment" begin
    for T in [Float32, Float64]
        sys = System{T}()
        mol = Molecule(sys)
        mol2 = Molecule(sys)
        chain = Chain(mol)
        chain2 = Chain(mol2)

        # constructors + parent
        frag = Fragment(chain, 1)
        @test frag isa Fragment{T}
        @test parent(frag) === sys
        @test parent_system(frag) === sys

        if T == Float32
            mol_ds = Molecule()
            chain_ds = Chain(mol_ds)
            frag_ds = Fragment(chain_ds, 1)
            parent(frag_ds) === default_system()
            parent_system(frag_ds) === default_system()

            Fragment(chain_ds, 1, "something", Properties(:a => "b"), Flags([:A]))
        end

        frag2 = Fragment(chain2, 1, "something", Properties(:a => 1), Flags([:A, :B]))

        #=
            Make sure we test for the correct number of fields.
            Add missing tests if the following test fails!
        =#
        @test length(getfield(frag, :_row)) == 7

        # getproperty
        @test frag.idx isa Int
        @test frag.number isa Int
        @test frag.number == 1
        @test frag.name isa String
        @test frag.name == ""
        @test frag.properties isa Properties
        @test frag.properties == Properties()
        @test frag.flags isa Flags
        @test frag.flags == Flags()

        @test frag._sys isa System{T}
        @test frag._row isa DataFrameRow
        
        @test_throws ErrorException frag.molecule_id
        @test_throws ErrorException frag.chain_id

        @test frag._row.molecule_id isa Int
        @test frag._row.molecule_id == mol.idx
        @test frag._row.chain_id isa Int
        @test frag._row.chain_id == chain.idx

        @test frag2.number == 1
        @test frag2.name == "something"
        @test frag2.properties == Properties(:a => 1)
        @test frag2.flags == Flags([:A, :B])
        @test frag2._row.molecule_id == mol2.idx
        @test frag2._row.chain_id == chain2.idx

        # setproperty!
        frag.number = 0
        @test frag.number == 0
        frag.name = "something else"
        @test frag.name == "something else"
        frag.properties = Properties(:first => "v1", :second => 99)
        @test length(frag.properties) == 2
        @test frag.properties[:first] == "v1"
        @test frag.properties[:second] == 99
        frag.flags = Flags([:C])
        @test length(frag.flags) == 1
        @test :C in frag.flags

        @test_throws ErrorException frag.molecule_id = 0
        @test_throws ErrorException frag.chain_id = 0

        # fragment_by_idx
        @test isnothing(fragment_by_idx(sys, -1))
        @test fragment_by_idx(sys, frag.idx) isa Fragment{T}
        @test fragment_by_idx(sys, frag.idx) == frag

        # _fragments
        df = _fragments(sys)
        @test df isa AbstractDataFrame
        @test size(df) == (2, length(fieldnames(_SystemFragmentTuple)))
        @test copy(df[1, 1:length(fieldnames(FragmentTuple))]) isa FragmentTuple
        @test size(_fragments(sys), 1) == 2
        @test size(_fragments(sys, molecule_id = -1), 1) == 0
        @test size(_fragments(sys, molecule_id = mol.idx), 1) == 1
        @test size(_fragments(sys, molecule_id = mol2.idx), 1) == 1
        @test size(_fragments(sys, molecule_id = nothing), 1) == 2
        @test size(_fragments(sys, chain_id = -1), 1) == 0
        @test size(_fragments(sys, chain_id = chain.idx), 1) == 1
        @test size(_fragments(sys, chain_id = chain2.idx), 1) == 1
        @test size(_fragments(sys, chain_id = nothing), 1) == 2
        @test size(_fragments(sys, molecule_id = -1, chain_id = chain.idx), 1) == 0
        @test size(_fragments(sys, molecule_id = mol.idx, chain_id = -1), 1) == 0
        @test size(_fragments(sys, molecule_id = mol.idx, chain_id = chain.idx), 1) == 1
        @test size(_fragments(sys, molecule_id = mol.idx, chain_id = nothing), 1) == 1
        @test size(_fragments(sys, molecule_id = nothing, chain_id = chain.idx), 1) == 1
        @test size(_fragments(sys, molecule_id = nothing, chain_id = nothing), 1) == 2

        # fragments_df
        df = fragments_df(sys)
        @test df isa AbstractDataFrame
        @test size(df) == (2, length(fieldnames(FragmentTuple)))
        @test copy(df[1, :]) isa FragmentTuple
        @test size(fragments_df(sys), 1) == 2
        @test size(fragments_df(sys, molecule_id = -1), 1) == 0
        @test size(fragments_df(sys, molecule_id = mol.idx), 1) == 1
        @test size(fragments_df(sys, molecule_id = mol2.idx), 1) == 1
        @test size(fragments_df(sys, molecule_id = nothing), 1) == 2
        @test size(fragments_df(sys, chain_id = -1), 1) == 0
        @test size(fragments_df(sys, chain_id = chain.idx), 1) == 1
        @test size(fragments_df(sys, chain_id = chain2.idx), 1) == 1
        @test size(fragments_df(sys, chain_id = nothing), 1) == 2
        @test size(fragments_df(sys, molecule_id = -1, chain_id = chain.idx), 1) == 0
        @test size(fragments_df(sys, molecule_id = mol.idx, chain_id = -1), 1) == 0
        @test size(fragments_df(sys, molecule_id = mol.idx, chain_id = chain.idx), 1) == 1
        @test size(fragments_df(sys, molecule_id = mol.idx, chain_id = nothing), 1) == 1
        @test size(fragments_df(sys, molecule_id = nothing, chain_id = chain.idx), 1) == 1
        @test size(fragments_df(sys, molecule_id = nothing, chain_id = nothing), 1) == 2

        # fragments
        fv = fragments(sys)
        @test fv isa Vector{Fragment{T}}
        @test length(fv) == 2
        @test length(fragments(sys)) == 2
        @test length(fragments(sys, molecule_id = -1)) == 0
        @test length(fragments(sys, molecule_id = mol.idx)) == 1
        @test length(fragments(sys, molecule_id = mol2.idx)) == 1
        @test length(fragments(sys, molecule_id = nothing)) == 2
        @test length(fragments(sys, chain_id = -1)) == 0
        @test length(fragments(sys, chain_id = chain.idx)) == 1
        @test length(fragments(sys, chain_id = chain2.idx)) == 1
        @test length(fragments(sys, chain_id = nothing)) == 2
        @test length(fragments(sys, molecule_id = -1, chain_id = chain.idx)) == 0
        @test length(fragments(sys, molecule_id = mol.idx, chain_id = -1)) == 0
        @test length(fragments(sys, molecule_id = mol.idx, chain_id = chain.idx)) == 1
        @test length(fragments(sys, molecule_id = mol.idx, chain_id = nothing)) == 1
        @test length(fragments(sys, molecule_id = nothing, chain_id = chain.idx)) == 1
        @test length(fragments(sys, molecule_id = nothing, chain_id = nothing)) == 2

        # eachfragment
        @test length(eachfragment(sys)) == 2
        @test length(eachfragment(sys, molecule_id = -1)) == 0
        @test length(eachfragment(sys, molecule_id = mol.idx)) == 1
        @test length(eachfragment(sys, molecule_id = mol2.idx)) == 1
        @test length(eachfragment(sys, molecule_id = nothing)) == 2
        @test length(eachfragment(sys, chain_id = -1)) == 0
        @test length(eachfragment(sys, chain_id = chain.idx)) == 1
        @test length(eachfragment(sys, chain_id = chain2.idx)) == 1
        @test length(eachfragment(sys, chain_id = nothing)) == 2
        @test length(eachfragment(sys, molecule_id = -1, chain_id = chain.idx)) == 0
        @test length(eachfragment(sys, molecule_id = mol.idx, chain_id = -1)) == 0
        @test length(eachfragment(sys, molecule_id = mol.idx, chain_id = chain.idx)) == 1
        @test length(eachfragment(sys, molecule_id = mol.idx, chain_id = nothing)) == 1
        @test length(eachfragment(sys, molecule_id = nothing, chain_id = chain.idx)) == 1
        @test length(eachfragment(sys, molecule_id = nothing, chain_id = nothing)) == 2

        # nfragments + push!
        @test nfragments(sys) isa Int
        @test nfragments(sys) == 2
        @test nfragments(sys, molecule_id = -1) == 0
        @test nfragments(sys, molecule_id = mol.idx) == 1
        @test nfragments(sys, molecule_id = mol2.idx) == 1
        @test nfragments(sys, molecule_id = nothing) == 2
        @test nfragments(sys, chain_id = -1) == 0
        @test nfragments(sys, chain_id = chain.idx) == 1
        @test nfragments(sys, chain_id = chain2.idx) == 1
        @test nfragments(sys, chain_id = nothing) == 2
        @test nfragments(sys, molecule_id = -1, chain_id = chain.idx) == 0
        @test nfragments(sys, molecule_id = mol.idx, chain_id = -1) == 0
        @test nfragments(sys, molecule_id = mol.idx, chain_id = chain.idx) == 1
        @test nfragments(sys, molecule_id = mol.idx, chain_id = nothing) == 1
        @test nfragments(sys, molecule_id = nothing, chain_id = chain.idx) == 1
        @test nfragments(sys, molecule_id = nothing, chain_id = nothing) == 2

        @test push!(chain, FragmentTuple(1)) === chain
        @test nfragments(sys) isa Int
        @test nfragments(sys) == 3
        @test nfragments(sys, molecule_id = -1) == 0
        @test nfragments(sys, molecule_id = mol.idx) == 2
        @test nfragments(sys, molecule_id = mol2.idx) == 1
        @test nfragments(sys, molecule_id = nothing) == 3
        @test nfragments(sys, chain_id = -1) == 0
        @test nfragments(sys, chain_id = chain.idx) == 2
        @test nfragments(sys, chain_id = chain2.idx) == 1
        @test nfragments(sys, chain_id = nothing) == 3
        @test nfragments(sys, molecule_id = -1, chain_id = chain.idx) == 0
        @test nfragments(sys, molecule_id = mol.idx, chain_id = -1) == 0
        @test nfragments(sys, molecule_id = mol.idx, chain_id = chain.idx) == 2
        @test nfragments(sys, molecule_id = mol.idx, chain_id = nothing) == 2
        @test nfragments(sys, molecule_id = nothing, chain_id = chain.idx) == 2
        @test nfragments(sys, molecule_id = nothing, chain_id = nothing) == 3

        # chain/molecule fragments
        mol3 = Molecule(sys)
        @test size(_fragments(mol3), 1) == 0
        @test _fragments(mol3) == _fragments(sys, molecule_id = mol3.idx)
        @test size(fragments_df(mol3), 1) == 0
        @test fragments_df(mol3) == fragments_df(sys, molecule_id = mol3.idx)
        @test size(fragments(mol3), 1) == 0
        @test fragments(mol3) == fragments(sys, molecule_id = mol3.idx)
        @test length(eachfragment(mol3)) == 0
        @test length(eachfragment(mol3)) == length(eachfragment(sys, molecule_id = mol3.idx))
        @test nfragments(mol3) == 0
        @test nfragments(mol3) == nfragments(sys, molecule_id = mol3.idx)

        chain3 = Chain(mol3)
        @test size(_fragments(chain3), 1) == 0
        @test _fragments(chain3) == _fragments(sys, chain_id = chain3.idx)
        @test size(fragments_df(chain3), 1) == 0
        @test fragments_df(chain3) == fragments_df(sys, chain_id = chain3.idx)
        @test size(fragments(chain3), 1) == 0
        @test fragments(chain3) == fragments(sys, chain_id = chain3.idx)
        @test length(eachfragment(chain3)) == 0
        @test length(eachfragment(chain3)) == length(eachfragment(sys, chain_id = chain3.idx))
        @test nfragments(chain3) == 0
        @test nfragments(chain3) == nfragments(sys, chain_id = chain3.idx)

        push!(chain3, FragmentTuple(1))
        @test size(_fragments(mol3), 1) == 1
        @test size(_fragments(mol3)) == size(_fragments(sys, molecule_id = mol3.idx))
        @test size(fragments_df(mol3), 1) == 1
        @test size(fragments_df(mol3)) == size(fragments_df(sys, molecule_id = mol3.idx))
        @test size(fragments(mol3), 1) == 1
        @test fragments(mol3) == fragments(sys, molecule_id = mol3.idx)
        @test length(eachfragment(mol3)) == 1
        @test length(eachfragment(mol3)) == length(eachfragment(sys, molecule_id = mol3.idx))
        @test nfragments(mol3) == 1
        @test nfragments(mol3) == nfragments(sys, molecule_id = mol3.idx)

        @test size(_fragments(chain3), 1) == 1
        @test size(_fragments(chain3)) == size(_fragments(sys, chain_id = chain3.idx))
        @test size(fragments_df(chain3), 1) == 1
        @test size(fragments_df(chain3)) == size(fragments_df(sys, chain_id = chain3.idx))
        @test size(fragments(chain3), 1) == 1
        @test fragments(chain3) == fragments(sys, chain_id = chain3.idx)
        @test length(eachfragment(chain3)) == 1
        @test length(eachfragment(chain3)) == length(eachfragment(sys, chain_id = chain3.idx))
        @test nfragments(chain3) == 1
        @test nfragments(chain3) == nfragments(sys, chain_id = chain3.idx)

        # fragment atoms
        @test size(_atoms(frag), 1) == 0
        @test _atoms(frag) == _atoms(sys, fragment_id = frag.idx)
        @test size(atoms_df(frag), 1) == 0
        @test atoms_df(frag) == atoms_df(sys, fragment_id = frag.idx)
        @test length(atoms(frag)) == 0
        @test atoms(frag) == atoms(sys, fragment_id = frag.idx)
        @test length(eachatom(frag)) == 0
        @test length(eachatom(frag)) == length(eachatom(sys, fragment_id = frag.idx))
        @test natoms(frag) == 0
        @test natoms(frag) == natoms(sys, fragment_id = frag.idx)

        @test push!(frag, AtomTuple{T}(1, Elements.H)) === frag
        @test size(_atoms(frag), 1) == 1
        @test size(_atoms(frag)) == size(_atoms(sys, fragment_id = frag.idx))
        @test size(atoms_df(frag), 1) == 1
        @test size(atoms_df(frag)) == size(atoms_df(sys, fragment_id = frag.idx))
        @test length(atoms(frag)) == 1
        @test atoms(frag) == atoms(sys, fragment_id = frag.idx)
        @test length(eachatom(frag)) == 1
        @test length(eachatom(frag)) == length(eachatom(sys, fragment_id = frag.idx))
        @test natoms(frag) == 1
        @test natoms(frag) == natoms(sys, fragment_id = frag.idx)

        # fragment bonds
        @test size(_bonds(frag), 1) == 0
        @test _bonds(frag) == _bonds(sys, fragment_id = frag.idx)
        @test size(bonds_df(frag), 1) == 0
        @test bonds_df(frag) == bonds_df(sys, fragment_id = frag.idx)
        @test length(bonds(frag)) == 0
        @test bonds(frag) == bonds(sys, fragment_id = frag.idx)
        @test length(eachbond(frag)) == 0
        @test length(eachbond(frag)) == length(eachbond(sys, fragment_id = frag.idx))
        @test nbonds(frag) == 0
        @test nbonds(frag) == nbonds(sys, fragment_id = frag.idx)

        @test push!(frag, BondTuple(
            Atom(frag, 1, Elements.H).idx,
            Atom(frag, 2, Elements.C).idx,
            BondOrder.Single
        )) === frag
        @test size(_bonds(frag), 1) == 1
        @test size(_bonds(frag)) == size(_bonds(sys, fragment_id = frag.idx))
        @test size(bonds_df(frag), 1) == 1
        @test size(bonds_df(frag)) == size(bonds_df(sys, fragment_id = frag.idx))
        @test length(bonds(frag)) == 1
        @test bonds(frag) == bonds(sys, fragment_id = frag.idx)
        @test length(eachbond(frag)) == 1
        @test length(eachbond(frag)) == length(eachbond(sys, fragment_id = frag.idx))
        @test nbonds(frag) == 1
        @test nbonds(frag) == nbonds(sys, fragment_id = frag.idx)
    end
end
