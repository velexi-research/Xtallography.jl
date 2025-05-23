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
Functions that support computations specific to rhombohedral lattices
"""
# --- Imports

# Standard library
using Logging

# --- Exports

# Types
export RhombohedralLatticeConstants, RhombohedralLatticeConstantDeltas

# Functions

# Constants
export RHOMBOHEDRAL_MIN_ANGLE, RHOMBOHEDRAL_MAX_ANGLE

# --- Constants

const RHOMBOHEDRAL_MIN_ANGLE = 0
const RHOMBOHEDRAL_MAX_ANGLE = 2π / 3

# --- Types

"""
    RhombohedralLatticeConstants

Lattice constants for a rhombohedral unit cell

Fields
======
* `a`: length of the edge of the unit cell

* `α`: angle between edges of the unit cell in the plane of the faces of the unit cell

Supertype: [`LatticeConstants`](@ref)
"""
struct RhombohedralLatticeConstants <: LatticeConstants{Rhombohedral}
    # Fields
    a::Float64
    α::Float64  # radians

    """
    Construct a set of rhombohedral lattice constants.
    """
    function RhombohedralLatticeConstants(a::Real, α::Real)

        # --- Enforce constraints

        if a <= 0
            throw(DomainError(a, "`a` must be positive"))
        end

        if α <= RHOMBOHEDRAL_MIN_ANGLE || α >= RHOMBOHEDRAL_MAX_ANGLE
            throw(DomainError(α, "`α` must satisfy 0 < α < 2π / 3"))
        end

        # --- Construct and return new RhombohedralLatticeConstants

        return new(a, α)
    end
end

"""
    RhombohedralLatticeConstantDeltas

Lattice constant deltas for a rhombohedral unit cell

Fields
======
* `Δa`: delta of the length of the edge of the unit cell

* `Δα`: delta of angle between edges of the unit cell in the plane of the faces of the unit
  cell

Supertype: [`LatticeConstantDeltas`](@ref)
"""
struct RhombohedralLatticeConstantDeltas <: LatticeConstantDeltas{Rhombohedral}
    # Fields
    Δa::Float64
    Δα::Float64

    """
    Construct a set of rhombohedral lattice constant deltas.
    """
    function RhombohedralLatticeConstantDeltas(Δa::Real, Δα::Real)
        return new(Δa, Δα)
    end
end

# --- Functions/Methods

# ------ LatticeConstants functions

function isapprox(
    x::RhombohedralLatticeConstants,
    y::RhombohedralLatticeConstants;
    atol::Real=0,
    rtol::Real=atol > 0 ? 0 : √eps(),
)
    return isapprox(x.a, y.a; atol=atol, rtol=rtol) &&
           isapprox(x.α, y.α; atol=atol, rtol=rtol)
end

function -(x::RhombohedralLatticeConstants, y::RhombohedralLatticeConstants)
    return RhombohedralLatticeConstantDeltas(x.a - y.a, x.α - y.α)
end

# ------ LatticeConstantDeltas functions

function isapprox(
    x::RhombohedralLatticeConstantDeltas,
    y::RhombohedralLatticeConstantDeltas;
    atol::Real=0,
    rtol::Real=atol > 0 ? 0 : √eps(),
)
    return isapprox(x.Δa, y.Δa; atol=atol, rtol=rtol) &&
           isapprox(x.Δα, y.Δα; atol=atol, rtol=rtol)
end

# ------ Unit cell computations

function basis(lattice_constants::RhombohedralLatticeConstants)
    # Get lattice constants
    a = lattice_constants.a
    α = lattice_constants.α

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

function volume(lattice_constants::RhombohedralLatticeConstants)
    # Get lattice constants
    a = lattice_constants.a
    α = lattice_constants.α

    # Compute volume
    return 2 * a^3 * sin(0.5 * α) * sqrt(sin(1.5 * α) * sin(0.5 * α))
end

function surface_area(lattice_constants::RhombohedralLatticeConstants)
    return 6 * lattice_constants.a^2 * sin(lattice_constants.α)
end

function conventional_cell(::Rhombohedral, unit_cell::UnitCell)
    # --- Check arguments

    conventional_cell_check_args(unit_cell)

    # --- Preparations

    # Get lattice constants
    a = unit_cell.lattice_constants.a
    α = unit_cell.lattice_constants.α

    # --- Compute IUCr conventional cell

    # Check limiting cases
    if α ≈ π / 3
        # cubic, face-centered, edge length `a` / sin(π/4)
        @debug "hR --> cF"
        return UnitCell(CubicLatticeConstants(a / SIN_PI_OVER_FOUR), face_centered)

    elseif α ≈ π / 2
        # cubic, primitive, edge length `a`
        @debug "hR --> cP"
        return UnitCell(CubicLatticeConstants(a), primitive)

    elseif α ≈ ACOS_MINUS_ONE_THIRD
        # cubic, body-centered, edge length `a` / sin(π/3)
        @debug "hR --> cI"
        return UnitCell(CubicLatticeConstants(a / SIN_PI_OVER_THREE), body_centered)
    end

    # Not a limiting case, so return unit cell with original lattice constants
    return UnitCell(unit_cell.lattice_constants, unit_cell.centering)
end
