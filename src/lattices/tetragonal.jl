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
Functions that support computations specific to tetragonal lattices
"""
# --- Imports

# Standard library
using Logging

# --- Exports

# Types
export TetragonalLatticeConstants

# Functions

# --- Types

"""
    TetragonalLatticeConstants

Lattice constants for a tetragonal unit cell

Fields
======
* `a`, `c`: lengths of the edges of the unit cell

Supertype: [`LatticeConstants`](@ref)
"""
struct TetragonalLatticeConstants <: LatticeConstants
    # Fields
    a::Float64
    c::Float64

    """
    Construct a set of tetragonal lattice constants.
    """
    function TetragonalLatticeConstants(a::Real, c::Real)

        # --- Enforce constraints

        if a <= 0
            throw(ArgumentError("`a` must be positive"))
        end

        if c <= 0
            throw(ArgumentError("`c` must be positive"))
        end

        # --- Construct and return new TetragonalLatticeConstants

        return new(a, c)
    end
end

# --- Functions/Methods

# ------ LatticeConstants functions

function isapprox(
    x::TetragonalLatticeConstants,
    y::TetragonalLatticeConstants;
    atol::Real=0,
    rtol::Real=atol > 0 ? 0 : √eps(),
)
    return isapprox(x.a, y.a; atol=atol, rtol=rtol) &&
           isapprox(x.c, y.c; atol=atol, rtol=rtol)
end

function lattice_system(::TetragonalLatticeConstants)
    return tetragonal
end

# ------ Unit cell computations

function basis(lattice_constants::TetragonalLatticeConstants)
    # Construct basis
    basis_a = Vector{Float64}([lattice_constants.a, 0, 0])
    basis_b = Vector{Float64}([0, lattice_constants.a, 0])
    basis_c = Vector{Float64}([0, 0, lattice_constants.c])

    return basis_a, basis_b, basis_c
end

function volume(lattice_constants::TetragonalLatticeConstants)
    return lattice_constants.a^2 * lattice_constants.c
end

function surface_area(lattice_constants::TetragonalLatticeConstants)
    return 2 * lattice_constants.a^2 + 4 * lattice_constants.a * lattice_constants.c
end

function conventional_cell(::Tetragonal, unit_cell::UnitCell)
    # --- Check arguments

    conventional_cell_arg_checks(unit_cell)

    # --- Preparations

    # Get lattice constants and centering
    a = unit_cell.lattice_constants.a
    c = unit_cell.lattice_constants.c

    centering = unit_cell.centering

    # --- Compute IUCr conventional cell

    # Check limiting cases
    if centering === primitive
        if a ≈ c
            # cubic, primitive, edge length `a`
            @debug "tP --> cP"
            return conventional_cell(UnitCell(CubicLatticeConstants(a), primitive))
        end

    elseif centering === body_centered
        if a ≈ c
            # cubic, body-centered, edge length `a`
            @debug "tI --> cI"
            return conventional_cell(UnitCell(CubicLatticeConstants(a), body_centered))

        elseif c * SIN_PI_OVER_FOUR ≈ a
            # cubic, face-centered, edge length `c`
            @debug "tI --> cF"
            return conventional_cell(UnitCell(CubicLatticeConstants(c), face_centered))
        end
    end

    # Not a limiting case, so return unit cell with original lattice constants
    return UnitCell(unit_cell.lattice_constants, centering)
end
