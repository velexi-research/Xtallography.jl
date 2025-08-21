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
Lattice system and symmetry types
"""
# --- Exports

# ------ Types

export LatticeSystem
export Triclinic, Monoclinic, Orthorhombic, Hexagonal, Rhombohedral, Tetragonal, Cubic
export triclinic, monoclinic, orthorhombic, hexagonal, rhombohedral, tetragonal, cubic

# ------ Functions/Methods

export lattice_system

# --- Types

"""
    LatticeSystem

Supertype for the seven lattice systems in 3D

Subtypes
========
[`Triclinic`](@ref), [`Monoclinic`](@ref), [`Orthorhombic`](@ref),
[`Hexagonal`](@ref), [`Rhombohedral`](@ref), [`Tetragonal`](@ref), [`Cubic`](@ref)
"""
abstract type LatticeSystem end

"""
    Triclinic

Type representing the triclinic lattice system that is the type of [`triclinic`](@ref)

Supertype: [`LatticeSystem`](@ref)
"""
struct Triclinic <: LatticeSystem end

"""
    triclinic

The singleton instance of type [`Triclinic`](@ref)
"""
const triclinic = Triclinic()

"""
    Monoclinic

Type representing the monoclinic lattice system that is the type of [`monoclinic`](@ref)

Supertype: [`LatticeSystem`](@ref)
"""
struct Monoclinic <: LatticeSystem end

"""
    monoclinic

The singleton instance of type [`Monoclinic`](@ref)
"""
const monoclinic = Monoclinic()

"""
    Orthorhombic

Type representing the orthorhombic lattice system that is the type of [`orthorhombic`](@ref)

Supertype: [`LatticeSystem`](@ref)
"""
struct Orthorhombic <: LatticeSystem end

"""
    orthorhombic

The singleton instance of type [`Orthorhombic`](@ref)
"""
const orthorhombic = Orthorhombic()

"""
    Hexagonal

Type representing the hexagonal lattice system that is the type of [`hexagonal`](@ref)

Supertype: [`LatticeSystem`](@ref)
"""
struct Hexagonal <: LatticeSystem end

"""
    hexagonal

The singleton instance of type [`Hexagonal`](@ref)
"""
const hexagonal = Hexagonal()

"""
    Rhombohedral

Type representing the rhombohedral lattice system that is the type of [`rhombohedral`](@ref)

Supertype: [`LatticeSystem`](@ref)
"""
struct Rhombohedral <: LatticeSystem end

"""
    rhombohedral

The singleton instance of type [`Rhombohedral`](@ref)
"""
const rhombohedral = Rhombohedral()

"""
    Tetragonal

Type representing the tetragonal lattice system that is the type of [`tetragonal`](@ref)

Supertype: [`LatticeSystem`](@ref)
"""
struct Tetragonal <: LatticeSystem end

"""
    tetragonal

The singleton instance of type [`Tetragonal`](@ref)
"""
const tetragonal = Tetragonal()

"""
    Cubic

Type representing the cubic lattice system that is the type of [`cubic`](@ref)

Supertype: [`LatticeSystem`](@ref)
"""
struct Cubic <: LatticeSystem end

"""
    cubic

The singleton instance of type [`Cubic`](@ref)
"""
const cubic = Cubic()

# --- Functions/Methods

"""
    lattice_system(lattice_constants::LatticeConstants) -> LatticeSystem

    lattice_system(Δlattice_constants::LatticeConstantDeltas) -> LatticeSystem

    lattice_system(unit_cell::UnitCell) -> LatticeSystem

Return the lattice system for a set of `lattice_constants`, `Δlattice_constants` or
`unit_cell`.

Return values
=============
- lattice system

Examples
========
```jldoctest
julia> lattice_system(HexagonalLatticeConstants(2, 4))
Hexagonal()

julia> lattice_system(CubicLatticeConstants(2))
Cubic()

julia> lattice_system(UnitCell(OrthorhombicLatticeConstants(2, 3, 4), F_centering))
Orthorhombic()
```
"""
function lattice_system end
