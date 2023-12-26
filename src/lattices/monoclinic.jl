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
Functions that support computations specific to monoclinic lattices
"""
# --- Imports

# Standard library
using Logging

# --- Exports

# Types
export MonoclinicLatticeConstants

# Functions
export convert_to_body_centering, convert_to_base_centering

# Constants
export MONOCLINIC_MIN_ANGLE, MONOCLINIC_MAX_ANGLE

# --- Constants

const MONOCLINIC_MIN_ANGLE = π / 2
const MONOCLINIC_MAX_ANGLE = 2π / 3

# --- Types

"""
    MonoclinicLatticeConstants

Lattice constants for a monoclinic unit cell

Fields
======
* `a`, `b`, `c`: lengths of the edges of the unit cell

* `β`: angle between edges of the unit cell in the plane of the face of the unit cell
  where the edges are not orthogonal

Supertype: [`LatticeConstants`](@ref)
"""
struct MonoclinicLatticeConstants <: LatticeConstants
    # Fields
    a::Float64  # arbitrary units
    b::Float64  # arbitrary units
    c::Float64  # arbitrary units
    β::Float64  # radians

    """
    Construct a set of monoclinic lattice constants.
    """
    function MonoclinicLatticeConstants(a::Real, b::Real, c::Real, β::Real)

        # --- Enforce constraints

        if a <= 0
            throw(ArgumentError("`a` must be positive"))
        end

        if b <= 0
            throw(ArgumentError("`b` must be positive"))
        end

        if c <= 0
            throw(ArgumentError("`c` must be positive"))
        end

        if β < 0 || β > π
            throw(ArgumentError("`β` must lie in the interval [0, π]"))
        end

        # --- Construct MonoclinicLatticeConstants object

        return new(a, b, c, β)
    end
end

# --- Functions/Methods

# ------ LatticeConstants functions

function isapprox(
    x::MonoclinicLatticeConstants,
    y::MonoclinicLatticeConstants;
    atol::Real=0,
    rtol::Real=atol > 0 ? 0 : √eps(),
)
    return isapprox(x.a, y.a; atol=atol, rtol=rtol) &&
           isapprox(x.b, y.b; atol=atol, rtol=rtol) &&
           isapprox(x.c, y.c; atol=atol, rtol=rtol) &&
           isapprox(x.β, y.β; atol=atol, rtol=rtol)
end

function lattice_system(::MonoclinicLatticeConstants)
    return monoclinic
end

function standardize(lattice_constants::MonoclinicLatticeConstants, centering::Centering)
    # --- Check arguments

    standardize_arg_checks(lattice_constants, centering)

    # --- Handle base-centering as a special case

    if centering === base_centered
        # Convert to a body-centered unit cell and standardize
        return standardize(convert_to_body_centering(lattice_constants), body_centered)
    end

    # --- Preparations

    # Extract lattice constants
    a = lattice_constants.a
    b = lattice_constants.b
    c = lattice_constants.c
    β = lattice_constants.β

    # --- Standardize lattice constants

    # Enforce IUCr conventions that a ≤ c and π/2 < β < π
    if a > c
        tmp = a
        a = c
        c = tmp
    end

    β = mod(β, π)
    if β < π / 2
        β = π - β
    end

    # --- Reduce the 2D unit cell in the plane normal to the b-axis so that the IUCr
    #     conventions for a, c, and β are satisfied.

    # Perform reduction
    if centering === primitive
        while -2 * c * cos(β) >= a

            # Attempt to compute equivalent lattice constants with smaller value of c.
            #
            # Note: calculation of c_alt assumes that a <= c and π/2 < β < π
            c_alt = sqrt(a^2 + c^2 + 2 * a * c * cos(β))

            if c_alt >= c
                break

            else
                β = π - asin_(sin(β) / c_alt * c)
                if c_alt < a
                    c = a
                    a = c_alt
                else
                    c = c_alt
                end
            end
        end
    elseif centering === body_centered
        while -c * cos(β) >= a

            # Attempt to compute equivalent lattice constants with smaller value of c.
            #
            # Note: calculation of c_alt assumes that a <= c and π/2 < β < π
            c_alt = sqrt(4 * a^2 + c^2 + 4 * a * c * cos(β))
            if c_alt >= c
                break

            else
                β = π - asin_(sin(β) / c_alt * c)
                if c_alt < a
                    c = a
                    a = c_alt
                else
                    c = c_alt
                end
            end
        end
    end

    return MonoclinicLatticeConstants(a, b, c, β), centering
end

# ------ Unit cell computations

import LinearAlgebra: det

function basis(lattice_constants::MonoclinicLatticeConstants)
    # Get lattice constants
    a = lattice_constants.a
    b = lattice_constants.b
    c = lattice_constants.c
    β = lattice_constants.β

    # Construct basis
    basis_a = Vector{Float64}([a, 0, 0])
    basis_b = Vector{Float64}([0, b, 0])
    basis_c = Vector{Float64}([c * cos(β), 0, c * sin(β)])

    return basis_a, basis_b, basis_c
end

function volume(lattice_constants::MonoclinicLatticeConstants)
    # Get lattice constants
    a = lattice_constants.a
    b = lattice_constants.b
    c = lattice_constants.c
    β = lattice_constants.β

    # Compute volume
    return a * b * c * sin(β)
end

function surface_area(lattice_constants::MonoclinicLatticeConstants)
    # Get lattice constants
    a = lattice_constants.a
    b = lattice_constants.b
    c = lattice_constants.c
    β = lattice_constants.β

    # Compute surface area
    return 2 * (a * b + b * c + c * a * sin(lattice_constants.β))
end

function conventional_cell(::Monoclinic, unit_cell::UnitCell)
    # --- Check arguments

    conventional_cell_arg_checks(unit_cell)

    # --- Preparations

    # Get standardized lattice constants and centering
    lattice_constants, centering = standardize(
        unit_cell.lattice_constants, unit_cell.centering
    )

    # Get lattice constants and centering
    a = lattice_constants.a
    b = lattice_constants.b
    c = lattice_constants.c
    β = lattice_constants.β

    # Initialize return values to the standardized monoclinic lattice constants
    iucr_lattice_constants = lattice_constants
    iucr_centering = centering

    # --- Compute IUCr conventional cell

    # Check limiting cases
    if centering === primitive
        if a ≈ -2 * c * cos(β)
            # Orthorhombic, C-centered
            @debug "mP --> oC (a = -2 * c * cos(β))"
            return conventional_cell(
                UnitCell(OrthorhombicLatticeConstants(a, 2 * c * sin(β), b), base_centered)
            )

        elseif a ≈ c
            # Orthorhombic, C-centered
            @debug "mP --> oC (a = c)"
            return conventional_cell(
                UnitCell(
                    OrthorhombicLatticeConstants(2 * c * cos(β / 2), 2 * c * sin(β / 2), b),
                    base_centered,
                ),
            )

        elseif β ≈ π / 2
            # Orthorhombic, primitive
            @debug "mP --> oP"
            return conventional_cell(
                UnitCell(OrthorhombicLatticeConstants(a, b, c), primitive)
            )
        end

    elseif centering === body_centered
        if β ≈ π / 2
            # Orthorhombic, body-centered
            @debug "mI --> oI"
            return conventional_cell(
                UnitCell(OrthorhombicLatticeConstants(a, b, c), body_centered)
            )

        elseif a ≈ c
            # Orthorhombic, face-centered
            @debug "mI --> oF"
            return conventional_cell(
                UnitCell(
                    OrthorhombicLatticeConstants(2 * a * cos(β / 2), 2 * a * sin(β / 2), b),
                    face_centered,
                ),
            )

        elseif a ≈ -c * cos(β)
            # Orthorhombic, C-centered
            @debug "mI --> oC"
            return conventional_cell(
                UnitCell(OrthorhombicLatticeConstants(b, c * sin(β), a), base_centered)
            )

        elseif a^2 + b^2 ≈ c^2 && a^2 + a * c * cos(β) ≈ b^2
            # Rhomohedral, primitive: α < π/3
            @debug "mI --> hR (α < π/3)"
            α = acos_(1 - 0.5 * b^2 / a^2)
            return conventional_cell(
                UnitCell(RhombohedralLatticeConstants(a, α), primitive)
            )

        elseif a^2 + b^2 ≈ c^2 && b^2 + a * c * cos(β) ≈ a^2
            # Rhomohedral, primitive: π/3 < α < π/2
            @debug "mI --> hR (π/3 < α < π/2)"
            α = acos_(1 - 0.5 * b^2 / a^2)
            return conventional_cell(
                UnitCell(RhombohedralLatticeConstants(a, α), primitive)
            )

        elseif c^2 + 3 * b^2 ≈ 9 * a^2 && c ≈ -3 * a * cos(β)
            # Rhomohedral, primitive: π/2 < α < acos(-1/3)
            @debug "mI --> hR (π/2 < α < acos(-1/3))"
            α = acos_((c^2 / a^2 - 3) / 6)
            return conventional_cell(
                UnitCell(RhombohedralLatticeConstants(a, α), primitive)
            )

        elseif a^2 + 3 * b^2 ≈ 9 * c^2 && a ≈ -3 * c * cos(β)
            # Rhomohedral, primitive: acos(-1/3) < α
            @debug "mI --> hR (acos(-1/3) < α)"
            α = acos_((a^2 / c^2 - 3) / 6)
            return conventional_cell(
                UnitCell(RhombohedralLatticeConstants(c, α), primitive)
            )
        end
    end

    # Not a limiting case, so return unit cell with standardized lattice constants
    return UnitCell(lattice_constants, centering)
end

"""
    convert_to_body_centering(
        lattice_constants::MonoclinicLatticeConstants
    ) -> MonoclinicLatticeConstants

Convert a base-centered monoclinic unit cell to body-centered monoclinic unit cell.

Return values
=============
- lattice constants for equivalent body-centered unit cell

Examples
========
```jldoctest
julia> lattice_constants = MonoclinicLatticeConstants(1.0, 2.0, 3.0, 3π / 5);

julia> body_centered_lattice_constants = convert_to_body_centering(lattice_constants);

julia> body_centered_lattice_constants.a ≈ 2.8541019662496847
true

julia> body_centered_lattice_constants.b ≈ 2
true

julia> body_centered_lattice_constants.c ≈ 3
true

julia> body_centered_lattice_constants.β ≈ 2.8018712454717734
true
```
"""
function convert_to_body_centering(lattice_constants::MonoclinicLatticeConstants)

    # Get lattice constants
    a = lattice_constants.a
    b = lattice_constants.b
    c = lattice_constants.c
    β = lattice_constants.β

    # Compute lattice constants for base-centered unit cell
    a_base = sqrt(a^2 + c^2 - 2 * a * c * abs(cos(β)))
    β_base = π - asin_(sin(β) / a_base * a)
    c_base = c

    if a_base < c_base
        return MonoclinicLatticeConstants(a_base, b, c_base, β_base)
    else
        return MonoclinicLatticeConstants(c_base, b, a_base, β_base)
    end
end

"""
    convert_to_base_centering(
        lattice_constants::MonoclinicLatticeConstants
    ) -> MonoclinicLatticeConstants

Convert a body-centered monoclinic unit cell to base-centered monoclinic unit cell.

Return values
=============
- lattice constants for equivalent base-centered unit cell

Examples
========
```jldoctest
julia> lattice_constants = MonoclinicLatticeConstants(1.0, 2.0, 3.0, 3π / 5);

julia> base_centered_lattice_constants = convert_to_base_centering(lattice_constants);

julia> base_centered_lattice_constants.a ≈ 2.8541019662496847
true

julia> base_centered_lattice_constants.b ≈ 2
true

julia> base_centered_lattice_constants.c ≈ 1
true

julia> base_centered_lattice_constants.β ≈ 1.5963584695539381
true
```
"""
function convert_to_base_centering(lattice_constants::MonoclinicLatticeConstants)

    # Get lattice constants
    a = lattice_constants.a
    b = lattice_constants.b
    c = lattice_constants.c
    β = lattice_constants.β

    # Compute lattice constants for base-centered unit cell
    a_base = sqrt(a^2 + c^2 - 2 * a * c * abs(cos(β)))
    c_base = a
    β_base = π - asin_(sin(β) / a_base * c)

    return MonoclinicLatticeConstants(a_base, b, c_base, β_base)
end

function is_supercell(
    lattice_constants_test::MonoclinicLatticeConstants,
    lattice_constants_ref::MonoclinicLatticeConstants;
    tol::Real=1e-3,
    max_index::Integer=3,
)
    # --- Check arguments

    if tol <= 0
        throw(ArgumentError("`tol` must be positive"))
    end

    if max_index <= 0
        throw(ArgumentError("`max_index` must be positive"))
    end

    # --- Preparations

    # Extract lattice constants
    a_test = lattice_constants_test.a
    b_test = lattice_constants_test.b
    c_test = lattice_constants_test.c
    β_test = lattice_constants_test.β

    a_ref = lattice_constants_ref.a
    b_ref = lattice_constants_ref.b
    c_ref = lattice_constants_ref.c
    β_ref = lattice_constants_ref.β

    # Construct metric tensors for 2D unit cells (a, c) and (a_ref, c_ref)
    U = [
        a_test c_test*cos(β_test)
        0 c_test*sin(β_test)
    ]
    G = U' * U

    U_ref = [
        a_ref c_ref*cos(β_ref)
        0 c_ref*sin(β_ref)
    ]
    G_ref = U_ref' * U_ref

    # --- Compare lattice constants

    # Check that b_test is an integer multiple of b_ref
    b_supercell_found = false
    multiplier = b_test / b_ref
    diff_from_int = abs(multiplier - round(multiplier))
    if diff_from_int < tol && !isapprox(multiplier, 1; atol=tol)
        b_supercell_found = true
    end

    # Check if (a_test, c_test, β_test) is a supercell of (a_ref, c_ref, β_ref)
    ac_supercell_found = false
    for h_a in (-max_index):max_index
        for l_a in (-max_index):max_index
            for h_c in (-max_index):max_index
                for l_c in (-max_index):max_index
                    # Construct transformation matrix from basis (a_ref, c_ref)
                    P = [
                        h_a h_c
                        l_a l_c
                    ]

                    # Skip (0, 0) and indices for unit cells with the same volume as
                    # (a_ref, c_ref, β_ref)
                    if abs(det(P)) <= 1
                        continue
                    end

                    # Compare metric tensors of (a, c) and test supercell basis
                    G_test = P' * G_ref * P
                    if G_test[2, 2] < G[1, 1]
                        tmp = G_test[1, 1]
                        G_test[1, 1] = G_test[2, 2]
                        G_test[2, 2] = tmp
                    end

                    G_diff = abs.(G_test - G)
                    if all(isapprox.(G_diff, 0; atol=tol))
                        ac_supercell_found = true
                    end
                end
            end
        end
    end

    if b_supercell_found || ac_supercell_found
        return true
    end

    return false
end
