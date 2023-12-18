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
export Centering, Primitive, BaseCentered, BodyCentered, FaceCentered

# Functions
export is_bravais_lattice

# --- Types

"""
    LatticeSystem

Supertype for the seven lattice systems in 3D

Subtypes
========
[`Triclinic`](@ref), [`Monoclinic`](@ref), [`Orthorhombic`](@ref), [`Tetragonal`](@ref),
[`Rhombohedral`](@ref), [`Hexagonal`](@ref), [`Cubic`](@ref)
"""
abstract type LatticeSystem end

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

"""
    BaseCentered

Type representing base centering

!!! note

    By convention, base-centering is on the C-face of the unit cell.

Supertype: [`Centering`](@ref)
"""
struct BaseCentered <: Centering end

"""
    BodyCentered

Type representing body centering

Supertype: [`Centering`](@ref)
"""
struct BodyCentered <: Centering end

"""
    FaceCentered

Type representing face centering

Supertype: [`Centering`](@ref)
"""
struct FaceCentered <: Centering end

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
    if lattice_system_ == Cubic() &&
        centering in (Primitive(), BodyCentered(), FaceCentered())
        return true
    elseif lattice_system_ == Tetragonal() && centering in (Primitive(), BodyCentered())
        return true
    elseif lattice_system_ == Orthorhombic() &&
        centering in (Primitive(), BaseCentered(), BodyCentered(), FaceCentered())
        return true
    elseif lattice_system_ == Hexagonal() && centering == Primitive()
        return true
    elseif lattice_system_ == Rhombohedral() && centering == Primitive()
        return true
    elseif lattice_system_ == Monoclinic() &&
        centering in (Primitive(), BaseCentered(), BodyCentered())
        return true
    elseif lattice_system_ == Triclinic() && centering == Primitive()
        return true
    end

    return false
end

function is_bravais_lattice(lattice_system_::Type{<:LatticeSystem}, centering::Centering)
    return is_bravais_lattice(lattice_system_(), centering)
end
