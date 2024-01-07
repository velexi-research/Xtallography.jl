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
Functions that support computations specific to cubic lattices
"""
# --- Exports

# Types
export CubicLatticeConstants, CubicLatticeConstantDeltas

# --- Types

"""
    CubicLatticeConstants

Lattice constants for a cubic unit cell

Fields
======
* `a`: length of the edge of the unit cell

Supertype: [`LatticeConstants`](@ref)
"""
struct CubicLatticeConstants <: LatticeConstants
    # Fields
    a::Float64

    """
    Construct a set of cubic lattice constants.
    """
    function CubicLatticeConstants(a::Real)

        # --- Enforce constraints

        if a <= 0
            throw(ArgumentError("`a` must be positive"))
        end

        # --- Construct and return new CubicLatticeConstants

        return new(a)
    end
end

"""
    CubicLatticeConstantDeltas

Lattice constant deltas for a cubic unit cell

Fields
======
* `Δa`: delta of the length of the edge of the unit cell

Supertype: [`LatticeConstantDeltas`](@ref)
"""
struct CubicLatticeConstantDeltas <: LatticeConstantDeltas
    # Fields
    Δa::Float64

    """
    Construct a set of cubic lattice constant deltas.
    """
    function CubicLatticeConstantDeltas(Δa::Real)
        return new(Δa)
    end
end

# --- Functions/Methods

# ------ LatticeConstants functions

function isapprox(
    x::CubicLatticeConstants,
    y::CubicLatticeConstants;
    atol::Real=0,
    rtol::Real=atol > 0 ? 0 : √eps(),
)
    return isapprox(x.a, y.a; atol=atol, rtol=rtol)
end

function -(x::CubicLatticeConstants, y::CubicLatticeConstants)
    return CubicLatticeConstantDeltas(x.a - y.a)
end

function lattice_system(::CubicLatticeConstants)
    return cubic
end

# ------ Unit cell computations

function basis(lattice_constants::CubicLatticeConstants)
    basis_a = Vector{Float64}([lattice_constants.a, 0, 0])
    basis_b = Vector{Float64}([0, lattice_constants.a, 0])
    basis_c = Vector{Float64}([0, 0, lattice_constants.a])

    return basis_a, basis_b, basis_c
end

function volume(lattice_constants::CubicLatticeConstants)
    return lattice_constants.a^3
end

function surface_area(lattice_constants::CubicLatticeConstants)
    return 6 * lattice_constants.a^2
end

function is_supercell(
    lattice_constants_test::CubicLatticeConstants,
    lattice_constants_ref::CubicLatticeConstants;
    tol::Real=1e-3,
    max_index::Integer=3,
)
    # --- Check arguments

    if tol <= 0
        throw(ArgumentError("`tol` must be positive"))
    end

    # --- Compare lattice constants

    # Check that a_test is an integer multiple of a_ref
    multiplier = lattice_constants_test.a / lattice_constants_ref.a
    diff_from_int = abs(multiplier - round(multiplier))
    return diff_from_int < tol && !isapprox(multiplier, 1; atol=tol)
end
