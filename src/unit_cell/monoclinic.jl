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
Monoclinic unit cell types and functions
"""
# --- Exports

# Types
export MonoclinicUnitCell, MonoclinicUnitCellDelta

# Functions
export convert_to_body_centering, convert_to_base_centering

# Constants
export MONOCLINIC_MIN_ANGLE, MONOCLINIC_MAX_ANGLE

# --- Constants

const MONOCLINIC_MIN_ANGLE = 0
const MONOCLINIC_MAX_ANGLE = π

# --- Types

# ------ MonoclinicUnitCell

"""
    MonoclinicUnitCell

Lattice constants for a monoclinic unit cell

Fields
======
* `a`, `b`, `c`: lengths of the edges of the unit cell

* `β`: angle between edges of the unit cell in the plane of the face of the unit cell
  where the edges are not orthogonal

* `symmetry`: unit cell symmetry
"""
const MonoclinicUnitCell = UnitCell{Monoclinic}

# Outer constructor
"""
    MonoclinicUnitCell(
        a::Real, b::Real, c::Real, β::Real;
        centering::Centering=primitive_centering,
        symmetry_elements::Union{Set,Vector,Nothing}=nothing
    )

Construct a MonoclinicUnitCell object from a set of lattice constants.

!!! note

    No constraints are imposed on `centering`. The unit cell does _not_ to be a valid
    Bravais lattice.

Keyword Arguments
=================
- `centering`: centering of unit cell

- `symmetry_elements`: symmetry elements of crystal
"""
function MonoclinicUnitCell(
    a::Real,
    b::Real,
    c::Real,
    β::Real;
    centering::Centering=primitive_centering,
    symmetry_elements::Union{Set,Vector,Nothing}=nothing,
)
    # --- Check arguments

    if a <= 0
        throw(DomainError(a, "`a` must be positive"))
    end

    if b <= 0
        throw(DomainError(b, "`b` must be positive"))
    end

    if c <= 0
        throw(DomainError(c, "`c` must be positive"))
    end

    if β <= MONOCLINIC_MIN_ANGLE || β >= MONOCLINIC_MAX_ANGLE
        throw(DomainError(β, "`β` must satisfy 0 < β < π"))
    end

    # symmetry elements
    # TODO

    # --- Construct and return MonoclinicUnitCell object

    return MonoclinicUnitCell(
        (a=a, b=b, c=c, β=β); centering=centering, symmetry_elements=symmetry_elements
    )
end

# ------ MonoclinicUnitCellDelta

"""
    MonoclinicUnitCellDelta

Lattice constant deltas for a monoclinic unit cell

Fields
======
* `Δa`, `Δb`, `Δc`: deltas of the lengths of the edges of the unit cell

* `Δβ`: delta of the angle between edges of the unit cell in the plane of the face of the
  unit cell where the edges are not orthogonal
"""
const MonoclinicUnitCellDelta = UnitCellDelta{Monoclinic}

# Outer constructors
"""
    MonoclinicUnitCellDelta(Δa::Real, Δb::Real, Δc::Real, Δβ::Real)

Construct a MonoclinicUnitCellDelta object from a set of lattice constant deltas.
"""
function MonoclinicUnitCellDelta(Δa::Real, Δb::Real, Δc::Real, Δβ::Real)
    Δlattice_constants = (Δa=Δa, Δb=Δb, Δc=Δc, Δβ=Δβ)
    return MonoclinicUnitCellDelta(Δlattice_constants)
end

# --- Functions/Methods

# ------ UnitCell methods

function standardize(unit_cell::MonoclinicUnitCell)
    # --- Check arguments

    standardize_check_args(unit_cell)

    # --- Preparations

    # Extract lattice constants
    lattice_constants_ = lattice_constants(unit_cell)
    a = lattice_constants_.a
    b = lattice_constants_.b
    c = lattice_constants_.c
    β = lattice_constants_.β

    # Extract centering
    centering_ = centering(unit_cell)

    # --- Handle base-centering as a special case

    if centering_ === base_centering
        # Convert to a body-centered unit cell and standardize
        return standardize(convert_to_body_centering(unit_cell))
    end

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
    if centering_ === primitive_centering
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
    elseif centering_ === body_centering
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

    return MonoclinicUnitCell(a, b, c, β; centering=centering_)
end

function basis(unit_cell::MonoclinicUnitCell)
    # Get lattice constants
    lattice_constants_ = lattice_constants(unit_cell)
    a = lattice_constants_.a
    b = lattice_constants_.b
    c = lattice_constants_.c
    β = lattice_constants_.β

    # Construct basis
    basis_a = Vector{Float64}([a, 0, 0])
    basis_b = Vector{Float64}([0, b, 0])
    basis_c = Vector{Float64}([c * cos(β), 0, c * sin(β)])

    return basis_a, basis_b, basis_c
end

function volume(unit_cell::MonoclinicUnitCell)
    # Get lattice constants
    lattice_constants_ = lattice_constants(unit_cell)
    a = lattice_constants_.a
    b = lattice_constants_.b
    c = lattice_constants_.c
    β = lattice_constants_.β

    # Compute volume
    return a * b * c * sin(β)
end

function surface_area(unit_cell::MonoclinicUnitCell)
    # Get lattice constants
    lattice_constants_ = lattice_constants(unit_cell)
    a = lattice_constants_.a
    b = lattice_constants_.b
    c = lattice_constants_.c
    β = lattice_constants_.β

    # Compute surface area
    return 2 * (a * b + b * c + c * a * sin(β))
end

function conventional_cell(::Monoclinic, unit_cell::UnitCell)
    # --- Check arguments

    conventional_cell_check_args(unit_cell)

    # --- Preparations

    # Standardize unit cell
    unit_cell = standardize(unit_cell)

    # Get lattice constants and centering
    lattice_constants_ = lattice_constants(unit_cell)
    a = lattice_constants_.a
    b = lattice_constants_.b
    c = lattice_constants_.c
    β = lattice_constants_.β

    centering_ = centering(unit_cell)

    # --- Compute IUCr conventional cell

    # Check limiting cases
    if centering_ === primitive_centering
        if a ≈ -2 * c * cos(β)
            # Orthorhombic, C-centering
            @debug "mP --> oC (a = -2 * c * cos(β))"
            return conventional_cell(
                OrthorhombicUnitCell(a, 2 * c * sin(β), b; centering=base_centering)
            )

        elseif a ≈ c
            # Orthorhombic, C-centering
            @debug "mP --> oC (a = c)"
            return conventional_cell(
                OrthorhombicUnitCell(
                    2 * c * cos(β / 2), 2 * c * sin(β / 2), b; centering=base_centering
                ),
            )

        elseif β ≈ π / 2
            # Orthorhombic, primitive-centering
            @debug "mP --> oP"
            return conventional_cell(
                OrthorhombicUnitCell(a, b, c; centering=primitive_centering)
            )
        end

    elseif centering_ === body_centering
        if β ≈ π / 2
            # Orthorhombic, body-centering
            @debug "mI --> oI"
            return conventional_cell(
                OrthorhombicUnitCell(a, b, c; centering=body_centering)
            )

        elseif a ≈ c
            # Orthorhombic, face-centering
            @debug "mI --> oF"
            return conventional_cell(
                OrthorhombicUnitCell(
                    2 * a * cos(β / 2), 2 * a * sin(β / 2), b; centering=face_centering
                ),
            )

        elseif a ≈ -c * cos(β)
            # Orthorhombic, C-centering
            @debug "mI --> oC"
            return conventional_cell(
                OrthorhombicUnitCell(b, c * sin(β), a; centering=base_centering)
            )

        elseif a^2 + b^2 ≈ c^2 && a^2 + a * c * cos(β) ≈ b^2
            # Rhomohedral, primitive-centering: α < π/3
            @debug "mI --> hR (α < π/3)"
            α = acos_(1 - 0.5 * b^2 / a^2)
            return conventional_cell(
                RhombohedralUnitCell(a, α; centering=primitive_centering)
            )

        elseif a^2 + b^2 ≈ c^2 && b^2 + a * c * cos(β) ≈ a^2
            # Rhomohedral, primitive-centering: π/3 < α < π/2
            @debug "mI --> hR (π/3 < α < π/2)"
            α = acos_(1 - 0.5 * b^2 / a^2)
            return conventional_cell(
                RhombohedralUnitCell(a, α; centering=primitive_centering)
            )

        elseif c^2 + 3 * b^2 ≈ 9 * a^2 && c ≈ -3 * a * cos(β)
            # Rhomohedral, primitive-centering: π/2 < α < acos(-1/3)
            @debug "mI --> hR (π/2 < α < acos(-1/3))"
            α = acos_((c^2 / a^2 - 3) / 6)
            return conventional_cell(
                RhombohedralUnitCell(a, α; centering=primitive_centering)
            )

        elseif a^2 + 3 * b^2 ≈ 9 * c^2 && a ≈ -3 * c * cos(β)
            # Rhomohedral, primitive-centering: acos(-1/3) < α
            @debug "mI --> hR (acos(-1/3) < α)"
            α = acos_((a^2 / c^2 - 3) / 6)
            return conventional_cell(
                RhombohedralUnitCell(c, α; centering=primitive_centering)
            )
        end
    end

    # Not a limiting case, so return standardized unit cell
    return unit_cell
end

"""
    convert_to_body_centering(unit_cell::MonoclinicUnitCell) -> MonoclinicUnitCell

Convert a base-centered monoclinic unit cell to body-centered monoclinic unit cell.

Return values
=============
unit cell of equivalent body-centered unit cell

Examples
========
```jldoctest
julia> unit_cell = MonoclinicUnitCell(1.0, 2.0, 3.0, 3π / 5);

julia> body_centered_unit_cell = convert_to_body_centering(unit_cell);

julia> lattice_constants_ = lattice_constants(body_centered_unit_cell);

julia> lattice_constants_.a ≈ 2.8541019662496847
true

julia> lattice_constants_.b ≈ 2
true

julia> lattice_constants_.c ≈ 3
true

julia> lattice_constants_.β ≈ 2.8018712454717734
true
```
"""
function convert_to_body_centering(unit_cell::MonoclinicUnitCell)

    # Get lattice constants
    lattice_constants_ = lattice_constants(unit_cell)
    a = lattice_constants_.a
    b = lattice_constants_.b
    c = lattice_constants_.c
    β = lattice_constants_.β

    # Compute lattice constants for base-centered unit cell
    a_base = sqrt(a^2 + c^2 - 2 * a * c * abs(cos(β)))
    β_base = π - asin_(sin(β) / a_base * a)
    c_base = c

    if a_base < c_base
        return MonoclinicUnitCell(a_base, b, c_base, β_base; centering=body_centering)
    else
        return MonoclinicUnitCell(c_base, b, a_base, β_base; centering=body_centering)
    end
end

"""
    convert_to_base_centering(unit_cell::MonoclinicUnitCell) -> MonoclinicUnitCell

Convert a body-centered monoclinic unit cell to base-centered monoclinic unit cell.

Return values
=============
lattice constants for equivalent base-centered unit cell

Examples
========
```jldoctest
julia> unit_cell = MonoclinicUnitCell(1.0, 2.0, 3.0, 3π / 5);

julia> base_centered_unit_cell = convert_to_base_centering(unit_cell);

julia> lattice_constants_ = lattice_constants(base_centered_unit_cell);

julia> lattice_constants_.a ≈ 2.8541019662496847
true

julia> lattice_constants_.b ≈ 2
true

julia> lattice_constants_.c ≈ 1
true

julia> lattice_constants_.β ≈ 1.5963584695539381
true
```
"""
function convert_to_base_centering(unit_cell::MonoclinicUnitCell)

    # Get lattice constants
    lattice_constants_ = lattice_constants(unit_cell)
    a = lattice_constants_.a
    b = lattice_constants_.b
    c = lattice_constants_.c
    β = lattice_constants_.β

    # Compute lattice constants for base-centered unit cell
    a_base = sqrt(a^2 + c^2 - 2 * a * c * abs(cos(β)))
    if a < c
        c_base = a
        β_base = π - asin_(sin(β) / a_base * c)
    else
        c_base = c
        β_base = π - asin_(sin(β) / a_base * a)
    end

    return MonoclinicUnitCell(a_base, b, c_base, β_base; centering=base_centering)
end

using LinearAlgebra: det

function is_supercell(
    unit_cell_test::MonoclinicUnitCell,
    unit_cell_ref::MonoclinicUnitCell;
    tol::Real=1e-3,
    max_index::Integer=3,
)
    # --- Check arguments

    if tol <= 0
        throw(DomainError(tol, "`tol` must be positive"))
    end

    if max_index <= 0
        throw(DomainError(max_index, "`max_index` must be positive"))
    end

    # --- Preparations

    # Extract lattice constants
    lattice_constants_test = lattice_constants(unit_cell_test)
    a_test = lattice_constants_test.a
    b_test = lattice_constants_test.b
    c_test = lattice_constants_test.c
    β_test = lattice_constants_test.β

    lattice_constants_ref = lattice_constants(unit_cell_ref)
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
