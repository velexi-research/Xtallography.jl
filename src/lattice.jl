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
Types and functions that support lattice computations
"""
# --- Exports

# Types
export LatticeSystem
export Triclinic, Monoclinic, Orthorhombic, Hexagonal, Rhombohedral, Tetragonal, Cubic
export Centering, Primitive, BaseCentered, BodyCentered, FaceCentered

# Constants
export BRAVAIS_LATTICES
export triclinic, monoclinic, orthorhombic, hexagonal, rhombohedral, tetragonal, cubic
export primitive, base_centered, body_centered, face_centered

# Functions
export is_bravais_lattice

# --- Types

# ------ Lattice systems

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

# ------ Lattice centerings

"""
    Centering

Supertype for the four lattice centerings in 3D

Subtypes
========
[`Primitive`](@ref), [`BaseCentered`](@ref), [`BodyCentered`](@ref), [`FaceCentered`](@ref)
"""
abstract type Centering end

"""
    Primitive

Type representing no centering that is the type of [`primitive`](@ref)

Supertype: [`Centering`](@ref)
"""
struct Primitive <: Centering end

"""
    primitive

The singleton instance of type [`Primitive`](@ref)
"""
const primitive = Primitive()

"""
    BaseCentered

Type representing base centering that is the type of [`base_centered`](@ref)

!!! note

    By convention, base-centering is on the C-face of the unit cell.

Supertype: [`Centering`](@ref)
"""
struct BaseCentered <: Centering end

"""
    base_centered

The singleton instance of type [`BaseCentered`](@ref)
"""
const base_centered = BaseCentered()

"""
    BodyCentered

Type representing body centering that is the type of [`body_centered`](@ref)

Supertype: [`Centering`](@ref)
"""
struct BodyCentered <: Centering end

"""
    body_centered

The singleton instance of type [`BodyCentered`](@ref)
"""
const body_centered = BodyCentered()

"""
    FaceCentered

Type representing face centering that is the type of [`face_centered`](@ref)

Supertype: [`Centering`](@ref)
"""
struct FaceCentered <: Centering end

"""
    face_centered

The singleton instance of type [`FaceCentered`](@ref)
"""
const face_centered = FaceCentered()

# --- Constants

# Lattice Types
"""
    BRAVAIS_LATTICES

List of valid Bravais lattices.
"""
const BRAVAIS_LATTICES = (
    (lattice_system=triclinic, centering=primitive),
    (lattice_system=monoclinic, centering=primitive),
    (lattice_system=monoclinic, centering=body_centered),
    (lattice_system=monoclinic, centering=base_centered),
    (lattice_system=orthorhombic, centering=primitive),
    (lattice_system=orthorhombic, centering=body_centered),
    (lattice_system=orthorhombic, centering=face_centered),
    (lattice_system=orthorhombic, centering=base_centered),
    (lattice_system=tetragonal, centering=primitive),
    (lattice_system=tetragonal, centering=body_centered),
    (lattice_system=rhombohedral, centering=primitive),
    (lattice_system=hexagonal, centering=primitive),
    (lattice_system=cubic, centering=primitive),
    (lattice_system=cubic, centering=body_centered),
    (lattice_system=cubic, centering=face_centered),
)

# --- Functions/Methods

"""
    lattice_system(lattice_constants::LatticesConstants) -> LatticeSystem

    lattice_system(Δlattice_constants::LatticesConstantDeltas) -> LatticeSystem

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

julia> lattice_system(UnitCell(OrthorhombicLatticeConstants(2, 3, 4), face_centered))
Orthorhombic()
```
"""
function lattice_system end

"""
    is_bravais_lattice(lattice_system::LatticeSystem, centering::Centering) -> Bool

    is_bravais_lattice(unit_cell::UnitCell) -> Bool

Determine if the unit cell defined by `unit_cell` or `lattice_system` and `centering` is
a valid Bravais lattice type.

Return values
=============
- `true` if `lattice_system` and `centering` define a valid Bravais lattice type; `false`
  otherwise

Examples
========
```jldoctest
julia> is_bravais_lattice(cubic, body_centered)
true

julia> is_bravais_lattice(cubic, base_centered)
false

julia> is_bravais_lattice(UnitCell(TetragonalLatticeConstants(2, 3), primitive))
true

julia> is_bravais_lattice(UnitCell(TetragonalLatticeConstants(2, 3), face_centered))
false
```
"""
function is_bravais_lattice(lattice_system::LatticeSystem, centering::Centering)
    return (lattice_system=lattice_system, centering=centering) in BRAVAIS_LATTICES
end
