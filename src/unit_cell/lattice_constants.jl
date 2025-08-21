#   Copyright 2023 Velexi Corporation
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

"""
Types and functions that support lattice constant computations
"""
# --- Exports

# Types
export LatticeConstants, LatticeConstantDeltas

# Functions
export (-)
export isapprox, convert
export standardize

# --- Types

"""
    LatticeConstants{T<:LatticeSystem}

Supertype for lattice constants for the seven lattice systems in 3D

Subtypes
========
[`TriclinicLatticeConstants`](@ref), [`MonoclinicLatticeConstants`](@ref),
[`OrthorhombicLatticeConstants`](@ref), [`TetragonalLatticeConstants`](@ref),
[`RhombohedralLatticeConstants`](@ref), [`HexagonalLatticeConstants`](@ref),
[`CubicLatticeConstants`](@ref)
"""
abstract type LatticeConstants{T<:LatticeSystem} end

"""
    LatticeConstantDeltas{T<:LatticeSystem}

Supertype for lattice constant deltas for the seven lattice systems in 3D

Subtypes
========
[`TriclinicLatticeConstantDeltas`](@ref), [`MonoclinicLatticeConstantDeltas`](@ref),
[`OrthorhombicLatticeConstantDeltas`](@ref), [`TetragonalLatticeConstantDeltas`](@ref),
[`RhombohedralLatticeConstantDeltas`](@ref), [`HexagonalLatticeConstantDeltas`](@ref),
[`CubicLatticeConstantDeltas`](@ref)
"""
abstract type LatticeConstantDeltas{T<:LatticeSystem} end

# --- Functions/Methods

# ------ Utility methods

function lattice_system(lattice_constants::LatticeConstants{T}) where {T<:LatticeSystem}
    return T()
end

function lattice_system(
    Δlattice_constant::LatticeConstantDeltas{T}
) where {T<:LatticeSystem}
    return T()
end

# ------ LatticeConstants methods

import Base.:(-)
import Base.isapprox
import Base.convert
using LinearAlgebra: LinearAlgebra

"""
    isapprox(x::LatticeConstants, y::LatticeConstants;
             atol::Real=0, rtol::Real=atol>0 ? 0 : √eps)

Inexact equality comparison between `LatticeConstants`. Two sets of lattice constants
compare equal if all lattice constant values are within the tolerance bounds. For instance,
`isapprox` returns `true` for `TetragonalLatticeConstants` if
`isapprox(x.a - y.a; atol=atol, rtol=rtol)` and `isapprox(x.b - y.b; atol=atol, rtol=rtol)`.
Returns `false` if `x` and `y` have different types.
"""
function isapprox(
    x::LatticeConstants, y::LatticeConstants; atol::Real=0, rtol::Real=atol > 0 ? 0 : √eps()
)
    # Default isapprox() implementation to allow comparison between lattice constants of
    # different lattice systems.
    return false
end

"""
    convert(::Type{<:Array}, lattice_constants::LatticeConstants)

Convert `lattice_constants` to an `Array`. For each lattice system, the order of lattice
constants in the array are in the conventional order:

* triclinic: `[a, b, c, α, β, γ]`
* monoclinic: `[a, b, c, β]`
* orthorhombic: `[a, b, c]`
* tetragonal: `[a, c]`
* rhombohedral: `[a, α]`
* hexagonal: `[a, c]`
* cubic: `[a]`
"""
function convert(type::Type{T}, lattice_constants::LatticeConstants) where {T<:Array}
    return convert(
        type,
        [
            getfield(lattice_constants, name) for
            name in fieldnames(typeof(lattice_constants))
        ],
    )
end

"""
    standardize(
        lattice_constants::LatticeConstants, centering::Centering
    ) -> (LatticeConstants, Centering)

Standardize the lattice constants and centering for the unit cell defined by
`lattice_constants` and `centering`.

Return values
=============
- standardized lattice constants and centering

Examples
========
```jldoctest
julia> standardize(OrthorhombicLatticeConstants(3, 2, 1), primitive)
(OrthorhombicLatticeConstants(1.0, 2.0, 3.0), Primitive())

julia> standardize(OrthorhombicLatticeConstants(3, 2, 1), base_centered)
(OrthorhombicLatticeConstants(2.0, 3.0, 1.0), BaseCentered())
```
"""
function standardize(lattice_constants::LatticeConstants, centering::Centering)
    # --- Check arguments

    standardize_check_args(lattice_constants, centering)

    # --- By default, no standardization is performed

    return lattice_constants, centering
end

"""
    standardize(lattice_constants::LatticeConstants) -> LatticeConstants

Standardize the lattice constants the primitive unit cell defined by `lattice_constants`.

Return values
=============
- standardized lattice constants for primitive unit cell

Examples
========
```jldoctest
julia> standardize(OrthorhombicLatticeConstants(3, 2, 1))
OrthorhombicLatticeConstants(1.0, 2.0, 3.0)
```
"""
function standardize(lattice_constants::LatticeConstants)
    # --- Check arguments

    standardize_check_args(lattice_constants, primitive)

    # Return standardized lattice constants for primitive unit cell
    standardized_lattice_constants, _ = standardize(lattice_constants, primitive)
    return standardized_lattice_constants
end

function standardize_check_args(lattice_constants::LatticeConstants, centering::Centering)
    # --- Check arguments

    if !is_bravais_lattice(lattice_system(lattice_constants), centering)
        throw(
            ArgumentError(
                "Invalid Bravais lattice: " *
                "(lattice_system=$(nameof(typeof(lattice_system(lattice_constants)))), " *
                "centering=$(nameof(typeof(centering))))",
            ),
        )
    end
end

# ------ LatticeConstantDeltas methods

import Base.isapprox
import Base.convert
using LinearAlgebra: LinearAlgebra

"""
    isapprox(Δx::LatticeConstantDeltas, Δy::LatticeConstantDeltas;
             atol::Real=0, rtol::Real=atol>0 ? 0 : √eps)

Inexact equality comparison between `LatticeConstantDeltas`. Two sets of lattice constant
deltas compare equal if all lattice constant values are within the tolerance bounds. For
instance, `isapprox` returns `true` for `TetragonalLatticeConstantDeltas` if
`isapprox(Δx.a - Δy.a; atol=atol, rtol=rtol)` and
`isapprox(Δx.b - Δy.b; atol=atol, rtol=rtol)`. Returns `false` if `Δx` and `Δy` have
different types.
"""
function isapprox(
    x::LatticeConstantDeltas,
    y::LatticeConstantDeltas;
    atol::Real=0,
    rtol::Real=atol > 0 ? 0 : √eps(),
)
    # Default isapprox() implementation to allow comparison between deltas of lattice
    # constants for different lattice systems.
    return false
end

"""
    convert(::Type{<:Array}, Δlattice_constants::LatticeConstantDeltas)

Convert `Δlattice_constants` to an `Array`. For each lattice system, the order of lattice
constant deltas in the array are in the conventional order:

* triclinic: `[Δa, Δb, Δc, Δα, Δβ, Δγ]`
* monoclinic: `[Δa, Δb, Δc, Δβ]`
* orthorhombic: `[Δa, Δb, Δc]`
* tetragonal: `[Δa, Δc]`
* rhombohedral: `[Δa, Δα]`
* hexagonal: `[Δa, Δc]`
* cubic: `[Δa]`
"""
function convert(type::Type{T}, Δlattice_constants::LatticeConstantDeltas) where {T<:Array}
    return convert(
        type,
        [
            getfield(Δlattice_constants, name) for
            name in fieldnames(typeof(Δlattice_constants))
        ],
    )
end
