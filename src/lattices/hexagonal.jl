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
Functions that support computations specific to hexagonal lattices
"""
# --- Exports

# Types
export HexagonalLatticeConstants

# --- Types

"""
    HexagonalLatticeConstants

Lattice constants for a hexagonal unit cell

Fields
======
* `a`, `c`: lengths of the edges of the unit cell

Supertype: [`LatticeConstants`](@ref)
"""
struct HexagonalLatticeConstants <: LatticeConstants
    # Fields
    a::Float64
    c::Float64

    """
    Construct a set of hexagonal lattice constants.
    """
    function HexagonalLatticeConstants(a::Real, c::Real)

        # --- Enforce constraints

        if a <= 0
            throw(ArgumentError("`a` must be positive"))
        end

        if c <= 0
            throw(ArgumentError("`c` must be positive"))
        end

        # --- Construct and return new HexagonalLatticeConstants

        return new(a, c)
    end
end

# --- Functions/Methods

# ------ LatticeConstants functions

function isapprox(
    x::HexagonalLatticeConstants,
    y::HexagonalLatticeConstants;
    atol::Real=0,
    rtol::Real=atol > 0 ? 0 : √eps(),
)
    return isapprox(x.a, y.a; atol=atol, rtol=rtol) &&
           isapprox(x.c, y.c; atol=atol, rtol=rtol)
end

function lattice_system(::HexagonalLatticeConstants)
    return Hexagonal()
end

# ------ Unit cell computations

function basis(lattice_constants::HexagonalLatticeConstants)
    # Get lattice constants
    a = lattice_constants.a
    c = lattice_constants.c

    # Construct basis
    basis_a = Vector{Float64}([a, 0, 0])
    basis_b = Vector{Float64}([-0.5 * a, SIN_PI_OVER_THREE * a, 0])
    basis_c = Vector{Float64}([0, 0, c])

    return basis_a, basis_b, basis_c
end

function volume(lattice_constants::HexagonalLatticeConstants)
    # Get lattice constants
    a = lattice_constants.a
    c = lattice_constants.c

    # Compute volume
    return a^2 * SIN_PI_OVER_THREE * c
end

function surface_area(lattice_constants::HexagonalLatticeConstants)
    # Get lattice constants
    a = lattice_constants.a
    c = lattice_constants.c

    # Compute surface area
    return 4 * a * c + 2 * a^2 * SIN_PI_OVER_THREE
end
