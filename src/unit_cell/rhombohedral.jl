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
Rhombohedral unit cell types and functions
"""
# --- Exports

# Types
export RhombohedralUnitCell, RhombohedralUnitCellDelta

# Constants
export RHOMBOHEDRAL_MIN_ANGLE, RHOMBOHEDRAL_MAX_ANGLE

# --- Constants

const RHOMBOHEDRAL_MIN_ANGLE = 0
const RHOMBOHEDRAL_MAX_ANGLE = 2π / 3

# --- Types

# ------ RhombohedralUnitCell

"""
    RhombohedralUnitCell

Lattice constants for a rhombohedral unit cell

Fields
======
* `a`: length of the edge of the unit cell

* `α`: angle between edges of the unit cell

* `symmetry`: unit cell symmetry
"""
const RhombohedralUnitCell = UnitCell{Rhombohedral}

# Outer constructor
"""
    RhombohedralUnitCell(
        a::Real, α::Real;
        centering::Centering=primitive_centering,
        symmetry_elements::Union{Set,Vector,Nothing}=nothing
    )

Construct a RhombohedralUnitCell object from a set of lattice constants.

!!! note

    No constraints are imposed on `centering`. The unit cell does _not_ to be a valid
    Bravais lattice.

Keyword Arguments
=================
- `centering`: centering of unit cell

- `symmetry_elements`: symmetry elements of crystal
"""
function RhombohedralUnitCell(
    a::Real,
    α::Real;
    centering::Centering=primitive_centering,
    symmetry_elements::Union{Set,Vector,Nothing}=nothing,
)
    # --- Check arguments

    if a <= 0
        throw(DomainError(a, "`a` must be positive"))
    end

    if α <= RHOMBOHEDRAL_MIN_ANGLE || α >= RHOMBOHEDRAL_MAX_ANGLE
        throw(DomainError(α, "`α` must satisfy 0 < α < 2π / 3"))
    end

    # symmetry elements
    # TODO

    # --- Construct and return RhombohedralUnitCell object

    return RhombohedralUnitCell(
        (a=a, α=α); centering=centering, symmetry_elements=symmetry_elements
    )
end

# ------ RhombohedralUnitCellDelta

"""
    RhombohedralUnitCellDelta

Lattice constant deltas for a rhombohedral unit cell

Fields
======
* `Δa`: delta of the length of the edge of the unit cell

* `Δα`: delta of angle between edges of the unit cell
"""
const RhombohedralUnitCellDelta = UnitCellDelta{Rhombohedral}

# Outer constructors
"""
    RhombohedralUnitCellDelta(Δa::Real, Δα::Real)

Construct a RhombohedralUnitCellDelta object from a set of lattice constant deltas.
"""
function RhombohedralUnitCellDelta(Δa::Real, Δα::Real)
    Δlattice_constants = (Δa=Δa, Δα=Δα)
    return RhombohedralUnitCellDelta(Δlattice_constants)
end

# --- Functions/Methods

# ------ UnitCell methods

function basis(unit_cell::RhombohedralUnitCell)
    # Get lattice constants
    lattice_constants_ = lattice_constants(unit_cell)
    a = lattice_constants_.a
    α = lattice_constants_.α

    # Construct basis
    basis_a = Vector{Float64}([a, 0, 0])
    basis_b = Vector{Float64}([a * cos(α), a * sin(α), 0])
    basis_c = Vector{Float64}([
        a * cos(α),
        a * cot(α) * (1 - cos(α)),
        a / cos(0.5 * α) * sqrt(sin(1.5 * α) * sin(0.5 * α)),
    ])

    return basis_a, basis_b, basis_c
end

function volume(unit_cell::RhombohedralUnitCell)
    # Get lattice constants
    lattice_constants_ = lattice_constants(unit_cell)
    a = lattice_constants_.a
    α = lattice_constants_.α

    # Compute volume
    return 2 * a^3 * sin(0.5 * α) * sqrt(sin(1.5 * α) * sin(0.5 * α))
end

function surface_area(unit_cell::RhombohedralUnitCell)
    lattice_constants_ = lattice_constants(unit_cell)
    return 6 * lattice_constants_.a^2 * sin(lattice_constants_.α)
end

function conventional_cell(::Rhombohedral, unit_cell::UnitCell)
    # --- Check arguments

    conventional_cell_check_args(unit_cell)

    # --- Preparations

    # Get lattice constants
    lattice_constants_ = lattice_constants(unit_cell)
    a = lattice_constants_.a
    α = lattice_constants_.α

    # --- Compute IUCr conventional cell

    # Check limiting cases
    if α ≈ π / 3
        # cubic, face-centering, edge length `a` / sin(π/4)
        @debug "hR --> cF"
        return CubicUnitCell(a / SIN_PI_OVER_FOUR; centering=face_centering)

    elseif α ≈ π / 2
        # cubic, primitive-centering, edge length `a`
        @debug "hR --> cP"
        return CubicUnitCell(a; centering=primitive_centering)

    elseif α ≈ ACOS_MINUS_ONE_THIRD
        # cubic, body-centering, edge length `a` / sin(π/3)
        @debug "hR --> cI"
        return CubicUnitCell(a / SIN_PI_OVER_THREE; centering=body_centering)
    end

    # Not a limiting case, so return a copy of the original unit cell
    return UnitCell(unit_cell)
end
