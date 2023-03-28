export AtomTuple, BondTuple, MoleculeTuple, ChainTuple, FragmentTuple, NucleotideTuple, ResidueTuple

"""
    const AtomTuple{T} = NamedTuple{...}

Tuple-based atom representation for `DataFrame` usage.

# Fields
 - `idx::Int`
 - `number::Int`
 - `element::ElementType`
 - `name::String`
 - `atomtype::String`
 - `r::Vector3{T}`
 - `v::Vector3{T}`
 - `F::Vector3{T}`
 - `has_velocity::Bool`
 - `has_force::Bool`
 - `properties::Properties`

# Constructors
```julia
AtomTuple(
    number::Int,
    element::ElementType;
    # keyword arguments
    idx::Int = 0,
    name::String = "",
    atomtype::String = "",
    r::Vector3{T} = Vector3{T}(0, 0, 0),
    v::Vector3{T} = Vector3{T}(0, 0, 0),
    F::Vector3{T} = Vector3{T}(0, 0, 0),
    formal_charge::Int = 0,
    charge::T = zero(T),
    radius::T = zero(T),
    has_velocity::Bool = false,
    has_force::Bool = false,
    properties::Properties = Properties()
)
```
Creates a new `AtomTuple{Float32}` with default values for all omitted fields.

```julia
AtomTuple{T}(
    number::Int,
    element::ElementType;
    # keyword arguments
    idx::Int = 0,
    name::String = "",
    atomtype::String = "",
    r::Vector3{T} = Vector3{T}(0, 0, 0),
    v::Vector3{T} = Vector3{T}(0, 0, 0),
    F::Vector3{T} = Vector3{T}(0, 0, 0),
    formal_charge::Int = 0,
    charge::T = zero(T),
    radius::T = zero(T),
    has_velocity::Bool = false,
    has_force::Bool = false,
    properties::Properties = Properties()
)
```
Creates a new `AtomTuple{T}` with default values for all omitted fields.
"""
const AtomTuple{T} = @NamedTuple begin
    idx::Int
    number::Int
    element::ElementType
    name::String
    atomtype::String
    r::Vector3{T}
    v::Vector3{T}
    F::Vector3{T}
    formal_charge::Int
    charge::T
    radius::T
    has_velocity::Bool
    has_force::Bool
    properties::Properties
end

@inline AtomTuple{T}(
    number::Int,
    element::ElementType;
    idx::Int = 0,
    name::String = "",
    atomtype::String = "",
    r::Vector3{T} = Vector3{T}(0, 0, 0),
    v::Vector3{T} = Vector3{T}(0, 0, 0),
    F::Vector3{T} = Vector3{T}(0, 0, 0),
    formal_charge::Int = 0,
    charge::T = zero(T),
    radius::T = zero(T),
    has_velocity::Bool = false,
    has_force::Bool = false,
    properties::Properties = Properties()
) where T = (
    idx = idx,
    number = number,
    element = element,
    name = name,
    atomtype = atomtype,
    r = r,
    v = v,
    F = F,
    formal_charge = formal_charge,
    charge = charge,
    radius = radius,
    has_velocity = has_velocity,
    has_force = has_force,
    properties = properties
)::AtomTuple{T}

@inline AtomTuple(
    number::Int,
    element::ElementType;
    kwargs...
) = AtomTuple{Float32}(number, element; kwargs...)

"""
    _with_idx(::AtomTuple{T}, idx::Int)
    _with_idx(::BondTuple{T}, idx::Int)
    _with_idx(::ChainTuple{T}, idx::Int)
    _with_idx(::FragmentTuple{T}, idx::Int)
    _with_idx(::MoleculeTuple{T}, idx::Int)
    _with_idx(::NucleotideTuple{T}, idx::Int)
    _with_idx(::ProteinTuple{T}, idx::Int)
    _with_idx(::ResidueTuple{T}, idx::Int)

Returns a copy of the given tuple with replaced `idx`.
"""
@inline function _with_idx(atom::AtomTuple{T}, idx::Int) where T
    ntuple(i -> i == 1 ? idx : atom[i], length(atom))
end

"""
    const BondTuple = NamedTuple{...}

Tuple-based bond representation for `DataFrame` usage.

# Fields
 - `idx::Int`
 - `a1::Int`
 - `a2::Int`
 - `order::BondOrderType`
 - `properties::Properties`

# Constructors
```julia
BondTuple(
    a1::Int,
    a2::Int,
    order::BondOrderType;
    # keyword arguments
    idx::Int = 0,
    properties::Properties = Properties()
)
```
Creates a new `BondTuple` with default values for all omitted fields.
"""
const BondTuple = @NamedTuple begin
    idx::Int
    a1::Int
    a2::Int
    order::BondOrderType
    properties::Properties
end

@inline BondTuple(
    a1::Int,
    a2::Int,
    order::BondOrderType;
    idx::Int = 0,
    properties::Properties = Properties()
) = (
    idx = idx,
    a1 = a1,
    a2 = a2,
    order = order,
    properties = properties
)::BondTuple

@inline function _with_idx(bond::BondTuple, idx::Int)
    ntuple(i -> i == 1 ? idx : bond[i], length(bond))
end

"""
    const MoleculeTuple = NamedTuple{...}

Tuple-based molecule representation for `DataFrame` usage.

# Fields
 - `idx::Int`
 - `name::String`
 - `properties::Properties`

# Constructors
```julia
MoleculeTuple(;
    # keyword arguments
    idx::Int = 0,
    name::String = "",
    properties::Properties = Properties()
)
```
Creates a new `MoleculeTuple` with default values for all omitted fields.
"""
const MoleculeTuple = @NamedTuple begin
    idx::Int
    name::String
    properties::Properties
end

@inline MoleculeTuple(;
    idx::Int = 0,
    name::String = "",
    properties::Properties = Properties()
) = (
    idx = idx,
    name = name,
    properties = properties
)::MoleculeTuple

@inline function _with_idx(mol::MoleculeTuple, idx::Int)
    ntuple(i -> i == 1 ? idx : mol[i], length(mol))
end

"""
    const ChainTuple = NamedTuple{...}

Tuple-based chain representation for `DataFrame` usage.

# Fields
 - `idx::Int`
 - `name::String`
 - `properties::Properties`

# Constructors
```julia
ChainTuple(;
    # keyword arguments
    idx::Int = 0,
    name::String = "",
    properties::Properties = Properties()
)
```
Creates a new `ChainTuple` with default values for all omitted fields.
"""
const ChainTuple = MoleculeTuple

"""
    const FragmentTuple = NamedTuple{...}

Tuple-based fragment representation for `DataFrame` usage.

# Fields
 - `idx::Int`
 - `number::Int`
 - `name::String`
 - `properties::Properties`

# Constructors
```julia
FragmentTuple(
    number::Int;
    # keyword arguments
    idx::Int = 0,
    name::String = "",
    properties::Properties = Properties()
)
```
Creates a new `FragmentTuple` with default values for all omitted fields.
"""
const FragmentTuple = @NamedTuple begin
    idx::Int
    number::Int
    name::String
    properties::Properties
end

@inline FragmentTuple(
    number::Int;
    idx::Int = 0,
    name::String = "",
    properties::Properties = Properties()
) = (
    idx = idx,
    number = number,
    name = name,
    properties = properties
)::FragmentTuple

@inline function _with_idx(frag::FragmentTuple, idx::Int)
    ntuple(i -> i == 1 ? idx : frag[i], length(frag))
end

"""
    const NucleotideTuple = NamedTuple{...}

Tuple-based nucleotide representation for `DataFrame` usage.

# Fields
 - `idx::Int`
 - `number::Int`
 - `name::String`
 - `properties::Properties`

# Constructors
```julia
NucleotideTuple(
    number::Int;
    # keyword arguments
    idx::Int = 0,
    name::String = "",
    properties::Properties = Properties()
)
```
Creates a new `NucleotideTuple` with default values for all omitted fields.
"""
const NucleotideTuple = FragmentTuple

"""
    const ResidueTuple = NamedTuple{...}

Tuple-based residue representation for `DataFrame` usage.

# Fields
 - `idx::Int`
 - `number::Int`
 - `type::AminoAcid`
 - `properties::Properties`

# Constructors
```julia
ResidueTuple(
    number::Int,
    type::AminoAcid;
    # keyword arguments
    idx::Int = 0,
    properties::Properties = Properties()
)
```
"""
const ResidueTuple = @NamedTuple begin
    idx::Int
    number::Int
    type::AminoAcid
    properties::Properties
end

@inline ResidueTuple(
    number::Int,
    type::AminoAcid;
    idx::Int = 0,
    properties::Properties = Properties()
) = (
    idx = idx,
    number = number,
    type = type,
    properties = properties
)::ResidueTuple

@inline function _with_idx(frag::ResidueTuple, idx::Int)
    ntuple(i -> i == 1 ? idx : frag[i], length(frag))
end
