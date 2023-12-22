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

Type representing the triclinic lattice system

Supertype: [`LatticeSystem`](@ref)
"""
struct Triclinic <: LatticeSystem end
const triclinic = Triclinic()

"""
    Monoclinic

Type representing the monoclinic lattice system

Supertype: [`LatticeSystem`](@ref)
"""
struct Monoclinic <: LatticeSystem end
const monoclinic = Monoclinic()

"""
    Orthorhombic

Type representing the orthorhombic lattice system

Supertype: [`LatticeSystem`](@ref)
"""
struct Orthorhombic <: LatticeSystem end
const orthorhombic = Orthorhombic()

"""
    Hexagonal

Type representing the hexagonal lattice system

Supertype: [`LatticeSystem`](@ref)
"""
struct Hexagonal <: LatticeSystem end
const hexagonal = Hexagonal()

"""
    Rhombohedral

Type representing the rhombohedral lattice system

Supertype: [`LatticeSystem`](@ref)
"""
struct Rhombohedral <: LatticeSystem end
const rhombohedral = Rhombohedral()

"""
    Tetragonal

Type representing the tetragonal lattice system

Supertype: [`LatticeSystem`](@ref)
"""
struct Tetragonal <: LatticeSystem end
const tetragonal = Tetragonal()

"""
    Cubic

Type representing the cubic lattice system

Supertype: [`LatticeSystem`](@ref)
"""
struct Cubic <: LatticeSystem end
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

Type representing no centering

Supertype: [`Centering`](@ref)
"""
struct Primitive <: Centering end
const primitive = Primitive()

"""
    BaseCentered

Type representing base centering

!!! note

    By convention, base-centering is on the C-face of the unit cell.

Supertype: [`Centering`](@ref)
"""
struct BaseCentered <: Centering end
const base_centered = BaseCentered()

"""
    BodyCentered

Type representing body centering

Supertype: [`Centering`](@ref)
"""
struct BodyCentered <: Centering end
const body_centered = BodyCentered()

"""
    FaceCentered

Type representing face centering

Supertype: [`Centering`](@ref)
"""
struct FaceCentered <: Centering end
const face_centered = FaceCentered()

# --- Constants

# Lattice Types
const BRAVAIS_LATTICES = (
    (lattice_system=Triclinic(), centering=Primitive()),
    (lattice_system=Monoclinic(), centering=Primitive()),
    (lattice_system=Monoclinic(), centering=BodyCentered()),
    (lattice_system=Monoclinic(), centering=BaseCentered()),
    (lattice_system=Orthorhombic(), centering=Primitive()),
    (lattice_system=Orthorhombic(), centering=BodyCentered()),
    (lattice_system=Orthorhombic(), centering=FaceCentered()),
    (lattice_system=Orthorhombic(), centering=BaseCentered()),
    (lattice_system=Tetragonal(), centering=Primitive()),
    (lattice_system=Tetragonal(), centering=BodyCentered()),
    (lattice_system=Rhombohedral(), centering=Primitive()),
    (lattice_system=Hexagonal(), centering=Primitive()),
    (lattice_system=Cubic(), centering=Primitive()),
    (lattice_system=Cubic(), centering=BodyCentered()),
    (lattice_system=Cubic(), centering=FaceCentered()),
)

# --- Functions/Methods

"""
    is_bravais_lattice(lattice_system::LatticeSystem, centering::Centering) -> Bool

    is_bravais_lattice(lattice_system::Type{<:LatticeSystem}, centering::Centering) -> Bool

    is_bravais_lattice(unit_cell::UnitCell) -> Bool

Determine if the unit cell defined by `unit_cell` or `lattice_system` and `centering` is
a valid Bravais lattice type.

Return values
=============
- `true` if `lattice_system` and `centering` define a valid Bravais lattice type; `false`
  otherwise

Examples
========
TODO
"""
function is_bravais_lattice(lattice_system_::LatticeSystem, centering::Centering)
    return (lattice_system=lattice_system_, centering=centering) in BRAVAIS_LATTICES
end

function is_bravais_lattice(lattice_system_::Type{<:LatticeSystem}, centering::Centering)
    return is_bravais_lattice(lattice_system_(), centering)
end
