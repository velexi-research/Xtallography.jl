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
Functions that support computations specific to triclinic lattices
"""
# --- Imports

# Standard library
using Logging

# External packages
using AngleBetweenVectors: angle
using Combinatorics: permutations

# --- Exports

# Types
export Triclinic
export TriclinicLatticeConstants

# Functions
export satisfies_triclinic_angle_constraints, is_triclinic_type_I_cell
export convert_to_mP, convert_to_mI, convert_to_mS

# Constants
export TRICLINIC_MIN_ANGLE, TRICLINIC_MAX_ANGLE

# --- Constants

const TRICLINIC_MIN_ANGLE = π / 3
const TRICLINIC_MAX_ANGLE = 2π / 3

# --- Types

"""
    Triclinic

Type representing the triclinic lattice system

Supertype: [`LatticeSystem`](@ref)
"""
struct Triclinic <: LatticeSystem end

"""
    TriclinicLatticeConstants

Lattice constants for a triclinic unit cell

Fields
======
* `a`, `b`, `c`: lengths of the edges of the unit cell

* `α`, `β`, `γ`: angles between edges of the unit cell in the planes of the
  faces of the unit cell

!!! note

    The constraints that valid triclinic unit cell angles must satisfy are _not_ enforced
    by the constructor. It is acceptable to construct `TriclinicLatticeConstants` with
    invalid values for `α`, `β`, and `γ`.

Supertype: [`LatticeConstants`](@ref)
"""
struct TriclinicLatticeConstants <: LatticeConstants
    # Fields
    a::Float64
    b::Float64
    c::Float64
    α::Float64  # radians
    β::Float64  # radians
    γ::Float64  # radians

    """
    Construct a set of triclinic lattice constants.
    """
    function TriclinicLatticeConstants(a::Real, b::Real, c::Real, α::Real, β::Real, γ::Real)

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

        if α < 0 || α > π
            throw(ArgumentError("`α` must lie in the interval [0, π]"))
        end

        if β < 0 || β > π
            throw(ArgumentError("`β` must lie in the interval [0, π]"))
        end

        if γ < 0 || γ > π
            throw(ArgumentError("`γ` must lie in the interval [0, π]"))
        end

        # --- Construct and return new TriclinicLatticeConstants

        return new(a, b, c, α, β, γ)
    end
end

# --- Functions/Methods

# ------ LatticeConstants functions

function isapprox(
    x::TriclinicLatticeConstants,
    y::TriclinicLatticeConstants;
    atol::Real=0,
    rtol::Real=atol > 0 ? 0 : √eps(),
)
    return isapprox(x.a, y.a; atol=atol, rtol=rtol) &&
           isapprox(x.b, y.b; atol=atol, rtol=rtol) &&
           isapprox(x.c, y.c; atol=atol, rtol=rtol) &&
           isapprox(x.α, y.α; atol=atol, rtol=rtol) &&
           isapprox(x.β, y.β; atol=atol, rtol=rtol) &&
           isapprox(x.γ, y.γ; atol=atol, rtol=rtol)
end

function lattice_system(::TriclinicLatticeConstants)
    return Triclinic
end

function standardize(lattice_constants::TriclinicLatticeConstants, centering::Centering)
    # --- Check arguments

    standardize_arg_checks(lattice_constants, centering)

    # --- Preparations

    # Extract lattice constants
    a = lattice_constants.a
    b = lattice_constants.b
    c = lattice_constants.c
    α = lattice_constants.α
    β = lattice_constants.β
    γ = lattice_constants.γ

    # --- Standardize lattice constants

    # Shift origin of basis to "homogeneous corner"
    if is_triclinic_type_I_cell(lattice_constants)
        # Type I

        if α > π / 2
            α = π - α
        end

        if β > π / 2
            β = π - β
        end

        if γ > π / 2
            γ = π - γ
        end

    else
        # Type II

        if α ≤ π / 2
            α = π - α
        end

        if β ≤ π / 2
            β = π - β
        end

        if γ ≤ π / 2
            γ = π - γ
        end
    end

    # Reorder edge lengths so that a < b < c
    #
    # Sorting algorithm: bubble sort
    if a > b
        # Swap edge lengths
        tmp = a
        a = b
        b = tmp

        # Swap angles
        tmp = α
        α = β
        β = tmp
    end

    if b > c
        # Swap edge lengths
        tmp = b
        b = c
        c = tmp

        # Swap angles
        tmp = β
        β = γ
        γ = tmp
    end

    if a > b
        # Swap edge lengths
        tmp = a
        a = b
        b = tmp

        # Swap angles
        tmp = α
        α = β
        β = tmp
    end

    # Reorder angles so that they are increasing order if the edge lengths are the same
    #
    # Sorting algorithm: bubble sort
    if a ≈ b ≈ c
        if α > β
            tmp = α
            α = β
            β = tmp
        end

        if β > γ
            tmp = β
            β = γ
            γ = tmp
        end

        if α > β
            tmp = α
            α = β
            β = tmp
        end

    elseif a ≈ b
        if α > β
            tmp = α
            α = β
            β = tmp
        end

    elseif b ≈ c
        if β > γ
            tmp = β
            β = γ
            γ = tmp
        end
    end

    return TriclinicLatticeConstants(a, b, c, α, β, γ), PRIMITIVE
end

"""
    satisfies_triclinic_angle_constraints(α::Real, β::Real, γ::Real) -> Bool

Determine whether `α`, `β`, and `γ` satisfy the angle constraints for triclinic lattices.

Return values
=============
- `true` if (`α`, `β`, `γ`) form a valid triple of angles for a triclinic unit cell;
  `false` otherwise

Examples
========
```jldoctest
julia> satisfies_triclinic_angle_constraints(π/4, π/5, π/6)
true
julia> satisfies_triclinic_angle_constraints(3π/4, 4π/5, 5π/6)
false
```
"""
function satisfies_triclinic_angle_constraints(α::Real, β::Real, γ::Real)
    return (0 < α + β + γ < 2π) &&
           (0 < α + β - γ < 2π) &&
           (0 < α - β + γ < 2π) &&
           (0 < -α + β + γ < 2π)
end

# ------ Unit cell computations

using LinearAlgebra: dot

function basis(lattice_constants::TriclinicLatticeConstants)
    # Get lattice constants
    a = lattice_constants.a
    b = lattice_constants.b
    c = lattice_constants.c
    α = lattice_constants.α
    β = lattice_constants.β
    γ = lattice_constants.γ

    # Construct basis
    basis_a = Vector{Float64}([a, 0, 0])
    basis_b = Vector{Float64}([b * cos(γ), b * sin(γ), 0])
    basis_c = Vector{Float64}([
        c * cos(β),
        c / sin(γ) * (cos(α) - cos(β) * cos(γ)),
        2 * c / sin(γ) * sqrt(
            sin(0.5 * (α + β + γ)) *
            sin(0.5 * (α + β - γ)) *
            sin(0.5 * (α - β + γ)) *
            sin(0.5 * (-α + β + γ)),
        ),
    ])

    return basis_a, basis_b, basis_c
end

function volume(lattice_constants::TriclinicLatticeConstants)
    # Get lattice constants
    a = lattice_constants.a
    b = lattice_constants.b
    c = lattice_constants.c
    α = lattice_constants.α
    β = lattice_constants.β
    γ = lattice_constants.γ

    # Compute volume
    return (2 * a * b * c) * sqrt(
        sin(0.5 * (α + β + γ)) *
        sin(0.5 * (α + β - γ)) *
        sin(0.5 * (α - β + γ)) *
        sin(0.5 * (-α + β + γ)),
    )
end

function surface_area(lattice_constants::TriclinicLatticeConstants)
    # Get lattice constants
    a = lattice_constants.a
    b = lattice_constants.b
    c = lattice_constants.c
    α = lattice_constants.α
    β = lattice_constants.β
    γ = lattice_constants.γ

    # Compute surface area
    return 2 * (a * b * sin(γ) + b * c * sin(α) + c * a * sin(β))
end

function iucr_conventional_cell(::Triclinic, unit_cell::UnitCell)
    # --- Check arguments

    iucr_conventional_cell_arg_checks(unit_cell)

    # --- Preparations

    # Get standardized lattice constants and centering
    lattice_constants, centering = standardize(
        unit_cell.lattice_constants, unit_cell.centering
    )

    # --- Compute IUCr conventional cell

    # Check limiting case: monoclinic, primitive
    try
        @debug "aP --> mP"
        return iucr_conventional_cell(UnitCell(convert_to_mP(lattice_constants), PRIMITIVE))
    catch error
        if !(error isa ErrorException) || (
            error.msg !=
            "The triclinic unit cell defined by `lattice_constants` is not equivalent " *
            "to a primitive monoclinic unit cell."
        )
            rethrow(error)
        end
    end

    # Check limiting case: monoclinic, body-centered
    try
        @debug "aP --> mI"
        return iucr_conventional_cell(UnitCell(convert_to_mI(lattice_constants), BODY))
    catch error
        if !(error isa ErrorException) || (
            error.msg !=
            "The triclinic unit cell defined by `lattice_constants` is not equivalent " *
            "to a body-centered monoclinic unit cell."
        )
            rethrow(error)
        end
    end

    # Check limiting case: monoclinic, base-centered
    try
        @debug "aP --> mS"
        body_centered_lattice_constants, centering = standardize(
            convert_to_mS(lattice_constants), BASE
        )

        return iucr_conventional_cell(UnitCell(body_centered_lattice_constants, BODY))
    catch error
        if !(error isa ErrorException) || (
            error.msg !=
            "The triclinic unit cell defined by `lattice_constants` is not equivalent " *
            "to a base-centered monoclinic unit cell."
        )
            rethrow(error)
        end
    end

    # Not a limiting case, so return unit cell with standardized lattice constants
    return UnitCell(lattice_constants, PRIMITIVE)
end

"""
    convert_to_mP(
        lattice_constants::TriclinicLatticeConstants
    ) -> MonoclinicLatticeConstants

Attempt to convert the triclinic unit cell defined by `lattice_constants` to an equivalent
primitive monoclinic unit cell.

Return values
=============
- lattice constants for the equivalent primitive monoclinic unit cell if one exists

Exceptions
==========
Throws an `ErrorException` if the triclinic unit cell defined by `lattice_constants` is not
equivalent to a primitive monoclinic unit cell.
"""
function convert_to_mP(lattice_constants::TriclinicLatticeConstants)
    # --- Preparations

    # Get lattice constants
    a = lattice_constants.a
    b = lattice_constants.b
    c = lattice_constants.c
    α = lattice_constants.α
    β = lattice_constants.β
    γ = lattice_constants.γ

    # --- Attempt to convert the triclinic unit cell to a primitive monoclinic unit cell

    if β ≈ π / 2 && γ ≈ π / 2
        # `a` is the unique symmetry direction of the monoclinic unit cell
        return standardize(MonoclinicLatticeConstants(b, a, c, α))

    elseif α ≈ π / 2 && γ ≈ π / 2
        # `b` is the unique symmetry direction of the monoclinic unit cell
        return standardize(MonoclinicLatticeConstants(a, b, c, β))

    elseif α ≈ π / 2 && β ≈ π / 2
        # `c` is the unique symmetry direction of the monoclinic unit cell
        return standardize(MonoclinicLatticeConstants(a, c, b, γ))
    end

    throw(
        ErrorException(
            "The triclinic unit cell defined by `lattice_constants` is not equivalent " *
            "to a primitive monoclinic unit cell.",
        ),
    )
end

"""
    convert_to_mI(
        lattice_constants::TriclinicLatticeConstants
    ) -> MonoclinicLatticeConstants

Attempt to convert the triclinic unit cell defined by `lattice_constants` to an equivalent
body-centered monoclinic unit cell.

Return values
=============
- lattice constants for the equivalent body-centered monoclinic unit cell if one exists;
  `nothing` otherwise

Exceptions
==========
Throws an `ErrorException` if the triclinic unit cell defined by `lattice_constants` is not
equivalent to a body-centered monoclinic unit cell.
"""
function convert_to_mI(lattice_constants::TriclinicLatticeConstants)
    # Notes
    # =====
    # This method adopts the following variable conventions.
    #
    # - Unless otherwise noted, lattice constants and basis vectors refer to the triclinic
    #   (not monoclinic) unit cell.
    #
    # - Lattice constants and basis vectors for the monoclinic unit cell are indicated by
    #   the "m_" prefix.

    # --- Preparations

    # Get basis for triclinic unit cell
    basis_a, basis_b, basis_c = basis(lattice_constants)

    # Get lattice constants for triclinic unit cell
    a = lattice_constants.a
    b = lattice_constants.b
    c = lattice_constants.c
    α = lattice_constants.α
    β = lattice_constants.β
    γ = lattice_constants.γ

    # Compute dot products
    a_dot_a = a^2
    b_dot_b = b^2
    c_dot_c = c^2

    a_dot_b = dot(basis_a, basis_b)
    b_dot_c = dot(basis_b, basis_c)
    c_dot_a = dot(basis_c, basis_a)

    # Initialize monoclinic lattice constants
    m_a = nothing
    m_b = nothing
    m_c = nothing
    m_β = nothing

    # Initialize monoclinic basis vectors
    m_basis_a = nothing
    m_basis_b = nothing
    m_basis_c = nothing

    # --- Attempt to convert the triclinic unit cell to a body-centered monoclinic unit cell

    # ------ Case: triclinic unit cell basis contains 2 monoclinic unit cell basis vectors
    #        and 1 body-centered lattice vector

    # Case: triclinic basis contains unique monoclinic symmetry direction
    try
        # Compute monoclinic basis vectors
        m_basis_a, m_basis_b, m_basis_c = convert_to_mI_case_1(lattice_constants)
        return convert_to_mI_basis_to_lattice_constants(m_basis_a, m_basis_b, m_basis_c)

    catch error
        if !(error isa ErrorException) || (
            error.msg !=
            "The triclinic basis vectors defined by `lattice_constants` do not " *
            "include the unique monoclinic symmetry direction, one other monoclinic " *
            "basis vector, and one body-centered lattice vector."
        )
            rethrow(error)
        end
    end

    # Case: triclinic basis does not contain unique monoclinic symmetry direction
    try
        # Compute monoclinic basis vectors
        m_basis_a, m_basis_b, m_basis_c = convert_to_mI_case_2(lattice_constants)
        return convert_to_mI_basis_to_lattice_constants(m_basis_a, m_basis_b, m_basis_c)

    catch error
        if !(error isa ErrorException) || (
            error.msg !=
            "The triclinic basis vectors defined by `lattice_constants` do not " *
            "include two monoclinic basis vectors in the B-face and one " *
            "body-centered lattice vector."
        )
            rethrow(error)
        end
    end

    # ------ Case: triclinic unit cell basis contains 1 monoclinic unit cell basis vector
    #        2 body-centered lattice vectors

    # Case: triclinic basis contains the unique monoclinic symmetry direction
    try
        # Compute monoclinic basis vectors
        m_basis_a, m_basis_b, m_basis_c = convert_to_mI_case_3(lattice_constants)
        return convert_to_mI_basis_to_lattice_constants(m_basis_a, m_basis_b, m_basis_c)

    catch error
        if !(error isa ErrorException) || (
            error.msg !=
            "The triclinic basis vectors defined by `lattice_constants` do not " *
            "include the unique monoclinic symmetry direction and two body-centered " *
            "lattice vectors."
        )
            rethrow(error)
        end
    end

    # Case: triclinic basis contains the unique monoclinic symmetry direction
    try
        # Compute monoclinic basis vectors
        m_basis_a, m_basis_b, m_basis_c = convert_to_mI_case_4(lattice_constants)
        return convert_to_mI_basis_to_lattice_constants(m_basis_a, m_basis_b, m_basis_c)

    catch error
        if !(error isa ErrorException) || (
            error.msg !=
            "The triclinic basis vectors defined by `lattice_constants` do not " *
            "include one monoclinic basis vector in the B-face and two body-centered " *
            "lattice vectors."
        )
            rethrow(error)
        end
    end

    # ------ Case: triclinic unit cell basis contains 3 body-centered lattice vectors

    try
        # Compute monoclinic basis vectors
        m_basis_a, m_basis_b, m_basis_c = convert_to_mI_case_5(lattice_constants)
        return convert_to_mI_basis_to_lattice_constants(m_basis_a, m_basis_b, m_basis_c)

    catch error
        if !(error isa ErrorException) || (
            error.msg !=
            "The triclinic basis vectors defined by `lattice_constants` do not " *
            "include three body-centered lattice vectors."
        )
            rethrow(error)
        end
    end

    # --- Throw exception if unit cell is not equivalent to a body-centered monoclinic unit

    if isnothing(m_basis_a) || isnothing(m_basis_b) || isnothing(m_basis_c)
        throw(
            ErrorException(
                "The triclinic unit cell defined by `lattice_constants` is not " *
                "equivalent to a body-centered monoclinic unit cell.",
            ),
        )
    end
end

function convert_to_mI_basis_to_lattice_constants(
    m_basis_a::Vector{<:Real}, m_basis_b::Vector{<:Real}, m_basis_c::Vector{<:Real}
)
    # Compute monoclinic lattice constants
    m_a = norm(m_basis_a)
    m_b = norm(m_basis_b)
    m_c = norm(m_basis_c)
    m_β = angle(m_basis_a, m_basis_c)

    # Check IUCr conventions
    if m_β < π / 2
        throw(ErrorException("m_β > 0"))
    end

    m_lattice_constants, _ = standardize(
        MonoclinicLatticeConstants(m_a, m_b, m_c, m_β), BODY
    )
    return m_lattice_constants
end

function convert_to_mI_case_1(lattice_constants::TriclinicLatticeConstants)
    # Notes
    # =====
    # - Case
    #   - triclinic unit cell basis contains 2 monoclinic unit cell basis vectors and
    #     1 body-centered lattice vector
    #
    #   - triclinic basis contains the unique monoclinic symmetry direction m_basis_b
    #
    # - This method adopts the same variable conventions as convert_to_mI().

    # --- Preparations

    # Get basis for triclinic unit cell
    basis_a, basis_b, basis_c = basis(lattice_constants)

    # Get lattice constants for triclinic unit cell
    a = lattice_constants.a
    b = lattice_constants.b
    c = lattice_constants.c
    α = lattice_constants.α
    β = lattice_constants.β
    γ = lattice_constants.γ

    # Compute dot products
    a_dot_a = a^2
    b_dot_b = b^2
    c_dot_c = c^2

    a_dot_b = dot(basis_a, basis_b)
    b_dot_c = dot(basis_b, basis_c)
    c_dot_a = dot(basis_c, basis_a)

    # Initialize monoclinic basis vectors
    m_basis_a = nothing
    m_basis_b = nothing
    m_basis_c = nothing

    # --- Attempt to convert the triclinic unit cell to a body-centered monoclinic unit cell

    if abs(a_dot_b) < a * b * COS_APPROX_ZERO
        # Case: `a` and `b` are basis vectors of the monoclinic unit cell

        if 2 * abs(c_dot_a) ≈ a_dot_a
            # Case: `a` is the unique monoclinic symmetry direction, `b` is basis vector
            #       on B-face of the monoclinic unit cell

            # Compute monoclinic basis
            m_basis_a = basis_b
            m_basis_b = basis_a

            local basis_c_minus_m_basis_b =
                basis_c - dot(basis_c, m_basis_b) / a_dot_a * m_basis_b

            local m_basis_c
            if 2 * abs(b_dot_c) > b_dot_b
                if b_dot_c > 0
                    m_basis_c = m_basis_a - 2 * basis_c_minus_m_basis_b
                else
                    m_basis_c = m_basis_a + 2 * basis_c_minus_m_basis_b
                end
            else
                m_basis_c = -m_basis_a + 2 * basis_c_minus_m_basis_b
            end

        elseif 2 * abs(b_dot_c) ≈ b_dot_b
            # Case: `b` is the unique monoclinic symmetry direction, `a` is basis vector
            #       on B-face of the monoclinic unit cell

            # Compute monoclinic basis
            m_basis_a = basis_a
            m_basis_b = basis_b

            local basis_c_minus_m_basis_b =
                basis_c - dot(basis_c, m_basis_b) / b_dot_b * m_basis_b

            local m_basis_c
            if 2 * abs(c_dot_a) > a_dot_a
                if c_dot_a > 0
                    m_basis_c = m_basis_a - 2 * basis_c_minus_m_basis_b
                else
                    m_basis_c = m_basis_a + 2 * basis_c_minus_m_basis_b
                end
            else
                m_basis_c = -m_basis_a + 2 * basis_c_minus_m_basis_b
            end
        end

    elseif abs(b_dot_c) < b * c * COS_APPROX_ZERO
        # Case: `b` and `c` are basis vectors of the monoclinic unit cell

        if 2 * abs(a_dot_b) ≈ b_dot_b
            # Case: `b` is the unique monoclinic symmetry direction, `c` is basis vector
            #       on B-face of the monoclinic unit cell

            # Compute monoclinic basis
            m_basis_a = basis_c
            m_basis_b = basis_b

            local basis_a_minus_m_basis_b =
                basis_a - dot(basis_a, m_basis_b) / b_dot_b * m_basis_b

            local m_basis_c
            if 2 * abs(c_dot_a) > c_dot_c
                if c_dot_a > 0
                    m_basis_c = m_basis_a - 2 * basis_a_minus_m_basis_b
                else
                    m_basis_c = m_basis_a + 2 * basis_a_minus_m_basis_b
                end
            else
                m_basis_c = -m_basis_a + 2 * basis_a_minus_m_basis_b
            end

        elseif 2 * abs(c_dot_a) ≈ c_dot_c
            # Case: `c` is the unique monoclinic symmetry direction, `b` is basis vector
            #       on B-face of the monoclinic unit cell

            # Compute monoclinic basis
            m_basis_a = basis_b
            m_basis_b = basis_c

            local basis_a_minus_m_basis_b =
                basis_a - dot(basis_a, m_basis_b) / c_dot_c * m_basis_b

            local m_basis_c
            if 2 * abs(a_dot_b) > b_dot_b
                if a_dot_b > 0
                    m_basis_c = m_basis_a - 2 * basis_a_minus_m_basis_b
                else
                    m_basis_c = m_basis_a + 2 * basis_a_minus_m_basis_b
                end
            else
                m_basis_c = -m_basis_a + 2 * basis_a_minus_m_basis_b
            end
        end

    elseif abs(c_dot_a) < c * a * COS_APPROX_ZERO
        # Case: `c` and `a` are basis vectors of the monoclinic unit cell

        if 2 * abs(b_dot_c) ≈ c_dot_c
            # Case: `c` is the unique monoclinic symmetry direction, `a` is basis vector
            #       on B-face of the monoclinic unit cell

            # Compute monoclinic basis
            m_basis_a = basis_a
            m_basis_b = basis_c

            local basis_b_minus_m_basis_b =
                basis_b - dot(basis_b, m_basis_b) / c_dot_c * m_basis_b

            local m_basis_c
            if 2 * abs(a_dot_b) > a_dot_a
                if a_dot_b > 0
                    m_basis_c = m_basis_a - 2 * basis_b_minus_m_basis_b
                else
                    m_basis_c = m_basis_a + 2 * basis_b_minus_m_basis_b
                end
            else
                m_basis_c = -m_basis_a + 2 * basis_b_minus_m_basis_b
            end

        elseif 2 * abs(a_dot_b) ≈ a_dot_a
            # Case: `a` is the unique monoclinic symmetry direction, `c` is basis vector
            #       on B-face of the monoclinic unit cell

            # Compute monoclinic basis
            m_basis_a = basis_c
            m_basis_b = basis_a

            local basis_b_minus_m_basis_b =
                basis_b - dot(basis_b, m_basis_b) / a_dot_a * m_basis_b

            local m_basis_c
            if 2 * abs(b_dot_c) > c_dot_c
                if b_dot_c > 0
                    m_basis_c = m_basis_a - 2 * basis_b_minus_m_basis_b
                else
                    m_basis_c = m_basis_a + 2 * basis_b_minus_m_basis_b
                end
            else
                m_basis_c = -m_basis_a + 2 * basis_b_minus_m_basis_b
            end
        end
    end

    # --- Construct return value

    if isnothing(m_basis_a) || isnothing(m_basis_b) || isnothing(m_basis_c)
        throw(
            ErrorException(
                "The triclinic basis vectors defined by `lattice_constants` do not " *
                "include the unique monoclinic symmetry direction, one other monoclinic " *
                "basis vector, and one body-centered lattice vector.",
            ),
        )
    end

    return m_basis_a, m_basis_b, m_basis_c
end

function convert_to_mI_case_2(lattice_constants::TriclinicLatticeConstants)

    # --- Preparations

    # Initialize monoclinic basis vectors
    m_basis_a = nothing
    m_basis_b = nothing
    m_basis_c = nothing

    # --- Attempt to convert the triclinic unit cell to a body-centered monoclinic unit cell

    try
        # Compute monoclinic basis vectors
        m_basis_a, m_basis_b, m_basis_c = convert_to_mI_case_2a(lattice_constants)

        # Compute monoclinic lattice constants
        m_a = norm(m_basis_a)
        m_b = norm(m_basis_b)
        m_c = norm(m_basis_c)
        m_β = angle(m_basis_a, m_basis_c)

    catch error
        if !(error isa ErrorException) || (
            error.msg !=
            "The triclinic basis vectors defined by `lattice_constants` do not " *
            "satisfy conditions for case 2a."
        )
            rethrow(error)
        end
    end

    try
        # Compute monoclinic basis vectors
        m_basis_a, m_basis_b, m_basis_c = convert_to_mI_case_2b(lattice_constants)

        # Compute monoclinic lattice constants
        m_a = norm(m_basis_a)
        m_b = norm(m_basis_b)
        m_c = norm(m_basis_c)
        m_β = angle(m_basis_a, m_basis_c)

    catch error
        if !(error isa ErrorException) || (
            error.msg !=
            "The triclinic basis vectors defined by `lattice_constants` do not " *
            "satisfy conditions for case 2b."
        )
            rethrow(error)
        end
    end

    # --- Construct return value

    if isnothing(m_basis_a) || isnothing(m_basis_b) || isnothing(m_basis_c)
        throw(
            ErrorException(
                "The triclinic basis vectors defined by `lattice_constants` do not " *
                "include two monoclinic basis vectors in the B-face and one " *
                "body-centered lattice vector.",
            ),
        )
    end

    return m_basis_a, m_basis_b, m_basis_c
end

function convert_to_mI_case_2a(lattice_constants::TriclinicLatticeConstants)
    # Notes
    # =====
    # - Case
    #   - triclinic unit cell basis contains 2 monoclinic unit cell basis vectors and
    #     1 body-centered lattice vector
    #
    #   - triclinic basis does not contain the unique monoclinic symmetry direction
    #     m_basis_b
    #
    # - This method adopts the same variable conventions as convert_to_mI().

    # --- Preparations

    # Get basis for triclinic unit cell
    basis_a, basis_b, basis_c = basis(lattice_constants)

    # Get lattice constants for triclinic unit cell
    a = lattice_constants.a
    b = lattice_constants.b
    c = lattice_constants.c
    α = lattice_constants.α
    β = lattice_constants.β
    γ = lattice_constants.γ

    # Compute dot products
    a_dot_a = a^2
    b_dot_b = b^2
    c_dot_c = c^2

    a_dot_b = dot(basis_a, basis_b)
    b_dot_c = dot(basis_b, basis_c)
    c_dot_a = dot(basis_c, basis_a)

    # Initialize monoclinic basis vectors
    m_basis_a = nothing
    m_basis_b = nothing
    m_basis_c = nothing

    # --- Attempt to convert the triclinic unit cell to a body-centered monoclinic unit cell

    if (
        2 * abs(a_dot_b) ≈ abs(b_dot_b + b_dot_c) &&
        2 * abs(c_dot_a) ≈ abs(c_dot_c + b_dot_c) &&
        abs(b_dot_c) > b * c * COS_APPROX_ZERO
    )
        # --- Case: `a` is vector to body-centered lattice point in monoclinic unit cell,
        #           `b` and `c` are `m_a` and `m_c`. `m_a` and `m_c` have the same sign
        #           in `a`.

        # Compute monoclinic basis vectors
        m_basis_a = basis_b

        if 2 * a_dot_b ≈ b^2 + b_dot_c
            m_basis_b = 2 * basis_a - basis_b - basis_c
        else
            m_basis_b = 2 * basis_a + basis_b + basis_c
        end

        if b_dot_c < 0
            m_basis_c = basis_c
        else
            m_basis_c = -basis_c
        end

    elseif (
        2 * abs(b_dot_c) ≈ abs(c_dot_c + c_dot_a) &&
        2 * abs(a_dot_b) ≈ abs(a_dot_a + c_dot_a) &&
        abs(c_dot_a) > c * a * COS_APPROX_ZERO
    )
        # --- Case: `b` is vector to body-centered lattice point in monoclinic unit cell,
        #           `a` and `c` are `m_a` and `m_c`. `m_a` and `m_c` have the same sign
        #           in `b`.

        # Compute monoclinic basis vectors
        m_basis_a = basis_c

        if 2 * b_dot_c ≈ c^2 + c_dot_a
            m_basis_b = 2 * basis_b - basis_c - basis_a
        elseif 2 * b_dot_c ≈ -c^2 - c_dot_a
            m_basis_b = 2 * basis_b + basis_c + basis_a
        end

        if c_dot_a < 0
            m_basis_c = basis_a
        else
            m_basis_c = -basis_a
        end

    elseif (
        2 * abs(c_dot_a) ≈ abs(a_dot_a + a_dot_b) &&
        2 * abs(b_dot_c) ≈ abs(b_dot_b + a_dot_b) &&
        abs(a_dot_b) > a * b * COS_APPROX_ZERO
    )
        # --- Case: `c` is vector to body-centered lattice point in monoclinic unit cell,
        #           `a` and `b` are `m_a` and `m_c`. `m_a` and `m_c` have the same sign
        #           in `c`.

        # Compute monoclinic basis vectors
        m_basis_a = basis_a

        if 2 * c_dot_a ≈ a^2 + a_dot_b
            m_basis_b = 2 * basis_c - basis_a - basis_b
        else
            m_basis_b = 2 * basis_c + basis_a + basis_b
        end

        if a_dot_b < 0
            m_basis_c = basis_b
        else
            m_basis_c = -basis_b
        end
    end

    # --- Construct return value

    if isnothing(m_basis_a) || isnothing(m_basis_b) || isnothing(m_basis_c)
        throw(
            ErrorException(
                "The triclinic basis vectors defined by `lattice_constants` do not " *
                "satisfy conditions for case 2a.",
            ),
        )
    end

    return m_basis_a, m_basis_b, m_basis_c
end

function convert_to_mI_case_2b(lattice_constants::TriclinicLatticeConstants)
    # Notes
    # =====
    # - Case
    #   - triclinic unit cell basis contains 2 monoclinic unit cell basis vectors and
    #     1 body-centered lattice vector
    #
    #   - triclinic basis does not contain the unique monoclinic symmetry direction
    #     m_basis_b
    #
    # - This method adopts the same variable conventions as convert_to_mI().

    # --- Preparations

    # Get basis for triclinic unit cell
    basis_a, basis_b, basis_c = basis(lattice_constants)

    # Get lattice constants for triclinic unit cell
    a = lattice_constants.a
    b = lattice_constants.b
    c = lattice_constants.c
    α = lattice_constants.α
    β = lattice_constants.β
    γ = lattice_constants.γ

    # Compute dot products
    a_dot_a = a^2
    b_dot_b = b^2
    c_dot_c = c^2

    a_dot_b = dot(basis_a, basis_b)
    b_dot_c = dot(basis_b, basis_c)
    c_dot_a = dot(basis_c, basis_a)

    # Initialize monoclinic basis vectors
    m_basis_a = nothing
    m_basis_b = nothing
    m_basis_c = nothing

    # --- Attempt to convert the triclinic unit cell to a body-centered monoclinic unit cell

    if (
        2 * abs(a_dot_b) ≈ abs(b_dot_b - b_dot_c) &&
        2 * abs(c_dot_a) ≈ abs(c_dot_c - b_dot_c) &&
        abs(b_dot_c) > b * c * COS_APPROX_ZERO
    )
        # --- Case: `a` is vector to body-centered lattice point in monoclinic unit cell,
        #           `b` and `c` are `m_a` and `m_c`. `m_a` and `m_c` have opposite signs
        #           in `a`.

        # Compute monoclinic basis vectors
        m_basis_a = basis_b

        if 2 * a_dot_b ≈ b^2 - b_dot_c
            m_basis_b = 2 * basis_a - basis_b + basis_c
        else
            m_basis_b = 2 * basis_a + basis_b - basis_c
        end

        if b_dot_c < 0
            m_basis_c = basis_c
        else
            m_basis_c = -basis_c
        end

    elseif (
        2 * abs(b_dot_c) ≈ abs(c_dot_c - c_dot_a) &&
        2 * abs(a_dot_b) ≈ abs(a_dot_a - c_dot_a) &&
        abs(c_dot_a) > c * a * COS_APPROX_ZERO
    )
        # --- Case: `b` is vector to body-centered lattice point in monoclinic unit cell,
        #           `a` and `c` are `m_a` and `m_c`. `m_a` and `m_c` have opposite signs
        #           in `b`.

        # Compute monoclinic basis vectors
        m_basis_a = basis_c

        if 2 * b_dot_c ≈ c^2 - c_dot_a
            m_basis_b = 2 * basis_b - basis_c + basis_a
        else
            m_basis_b = 2 * basis_b + basis_c - basis_a
        end

        if c_dot_a < 0
            m_basis_c = basis_a
        else
            m_basis_c = -basis_a
        end

    elseif (
        2 * abs(c_dot_a) ≈ abs(a_dot_a - a_dot_b) &&
        2 * abs(b_dot_c) ≈ abs(b_dot_b - a_dot_b) &&
        abs(a_dot_b) > a * b * COS_APPROX_ZERO
    )
        # --- Case: `c` is vector to body-centered lattice point in monoclinic unit cell,
        #           `a` and `b` are `m_a` and `m_c`. `m_a` and `m_c` have opposite signs
        #           in `c`.

        # Compute monoclinic basis vectors
        m_basis_a = basis_a

        if 2 * c_dot_a ≈ a^2 - a_dot_b
            m_basis_b = 2 * basis_c - basis_a + basis_b
        else
            m_basis_b = 2 * basis_c + basis_a - basis_b
        end

        if a_dot_b < 0
            m_basis_c = basis_b
        else
            m_basis_c = -basis_b
        end
    end

    # --- Construct return value

    if isnothing(m_basis_a) || isnothing(m_basis_b) || isnothing(m_basis_c)
        throw(
            ErrorException(
                "The triclinic basis vectors defined by `lattice_constants` do not " *
                "satisfy conditions for case 2b.",
            ),
        )
    end

    return m_basis_a, m_basis_b, m_basis_c
end

function convert_to_mI_case_3(lattice_constants::TriclinicLatticeConstants)
    # Notes
    # =====
    # - Case
    #   - triclinic unit cell basis contains 1 monoclinic unit cell basis vector and
    #     2 body-centered lattice vectors
    #
    #   - triclinic basis contains the unique monoclinic symmetry direction m_basis_b
    #
    # - This method adopts the same variable conventions as convert_to_mI().

    # --- Preparations

    # Get basis for triclinic unit cell
    basis_a, basis_b, basis_c = basis(lattice_constants)

    # Get lattice constants for triclinic unit cell
    a = lattice_constants.a
    b = lattice_constants.b
    c = lattice_constants.c
    α = lattice_constants.α
    β = lattice_constants.β
    γ = lattice_constants.γ

    # Compute dot products
    a_dot_a = a^2
    b_dot_b = b^2
    c_dot_c = c^2

    a_dot_b = dot(basis_a, basis_b)
    b_dot_c = dot(basis_b, basis_c)
    c_dot_a = dot(basis_c, basis_a)

    # Initialize monoclinic basis vectors
    m_basis_a = nothing
    m_basis_b = nothing
    m_basis_c = nothing

    # --- Attempt to convert the triclinic unit cell to a body-centered monoclinic unit cell

    if abs(a_dot_b) ≈ abs(c_dot_a) ≈ 0.5 * a_dot_a
        # Case: `a` is the unique symmetry direction of the monoclinic unit cell

        # Compute monoclinic basis vectors
        m_basis_b = basis_a

        m_basis_a = basis_b + basis_c
        m_basis_a -= dot(m_basis_a, m_basis_b) / a_dot_a * m_basis_b

        if b < c
            m_basis_c = basis_b - basis_c
        else
            m_basis_c = -basis_b + basis_c
        end
        m_basis_c -= dot(m_basis_c, m_basis_b) / a_dot_a * m_basis_b

    elseif abs(b_dot_c) ≈ abs(a_dot_b) ≈ 0.5 * b_dot_b
        # Case: `b` is the unique symmetry direction of the monoclinic unit cell

        # Compute monoclinic basis vectors
        m_basis_b = basis_b

        m_basis_a = basis_a + basis_c
        m_basis_a -= dot(m_basis_a, m_basis_b) / b_dot_b * m_basis_b

        if a < c
            m_basis_c = basis_a - basis_c
        else
            m_basis_c = -basis_a + basis_c
        end
        m_basis_c -= dot(m_basis_c, m_basis_b) / b_dot_b * m_basis_b

    elseif abs(c_dot_a) ≈ abs(b_dot_c) ≈ 0.5 * c_dot_c
        # Case: `c` is the unique symmetry direction of the monoclinic unit cell

        # Compute monoclinic basis vectors
        m_basis_b = basis_c

        m_basis_a = basis_a + basis_b
        m_basis_a -= dot(m_basis_a, m_basis_b) / c_dot_c * m_basis_b

        if a < b
            m_basis_c = basis_a - basis_b
        else
            m_basis_c = -basis_a + basis_b
        end
        m_basis_c -= dot(m_basis_c, m_basis_b) / c_dot_c * m_basis_b
    end

    # --- Construct return value

    if isnothing(m_basis_a) || isnothing(m_basis_b) || isnothing(m_basis_c)
        throw(
            ErrorException(
                "The triclinic basis vectors defined by `lattice_constants` do not " *
                "include the unique monoclinic symmetry direction and two body-centered " *
                "lattice vectors.",
            ),
        )
    end

    return m_basis_a, m_basis_b, m_basis_c
end

function convert_to_mI_case_4(lattice_constants::TriclinicLatticeConstants)
    # --- Preparations

    # Initialize monoclinic basis vectors
    m_basis_a = nothing
    m_basis_b = nothing
    m_basis_c = nothing

    # --- Attempt to convert the triclinic unit cell to a body-centered monoclinic unit cell

    try
        # Compute monoclinic basis vectors
        m_basis_a, m_basis_b, m_basis_c = convert_to_mI_case_4a(lattice_constants)

        # Compute monoclinic lattice constants
        m_a = norm(m_basis_a)
        m_b = norm(m_basis_b)
        m_c = norm(m_basis_c)
        m_β = angle(m_basis_a, m_basis_c)

    catch error
        if !(error isa ErrorException) || (
            error.msg !=
            "The triclinic basis vectors defined by `lattice_constants` do not " *
            "satisfy conditions for case 4a."
        )
            rethrow(error)
        end
    end

    try
        # Compute monoclinic basis vectors
        m_basis_a, m_basis_b, m_basis_c = convert_to_mI_case_4b(lattice_constants)

        # Compute monoclinic lattice constants
        m_a = norm(m_basis_a)
        m_b = norm(m_basis_b)
        m_c = norm(m_basis_c)
        m_β = angle(m_basis_a, m_basis_c)

    catch error
        if !(error isa ErrorException) || (
            error.msg !=
            "The triclinic basis vectors defined by `lattice_constants` do not " *
            "satisfy conditions for case 4b."
        )
            rethrow(error)
        end
    end

    try
        # Compute monoclinic basis vectors
        m_basis_a, m_basis_b, m_basis_c = convert_to_mI_case_4c(lattice_constants)

        # Compute monoclinic lattice constants
        m_a = norm(m_basis_a)
        m_b = norm(m_basis_b)
        m_c = norm(m_basis_c)
        m_β = angle(m_basis_a, m_basis_c)

    catch error
        if !(error isa ErrorException) || (
            error.msg !=
            "The triclinic basis vectors defined by `lattice_constants` do not " *
            "satisfy conditions for case 4c."
        )
            rethrow(error)
        end
    end

    # --- Construct return value

    if isnothing(m_basis_a) || isnothing(m_basis_b) || isnothing(m_basis_c)
        throw(
            ErrorException(
                "The triclinic basis vectors defined by `lattice_constants` do not " *
                "include one monoclinic basis vector in the B-face and two body-centered " *
                "lattice vectors.",
            ),
        )
    end

    return m_basis_a, m_basis_b, m_basis_c
end

function convert_to_mI_case_4a(lattice_constants::TriclinicLatticeConstants)
    # Notes
    # =====
    # - Case
    #   - triclinic unit cell basis contains 1 monoclinic unit cell basis vector and
    #     2 body-centered lattice vectors
    #
    #   - triclinic basis does not contain the unique monoclinic symmetry direction
    #     m_basis_b
    #
    #   - norm(basis_b) = norm(basis_c) and
    #     |dot(basis_a, basis_b)| = |dot(basis_a, basis_c)| != 0.5 a^2
    #
    # - This method adopts the same variable conventions as convert_to_mI().

    # --- Preparations

    # Get basis for triclinic unit cell
    basis_a, basis_b, basis_c = basis(lattice_constants)

    # Get lattice constants for triclinic unit cell
    a = lattice_constants.a
    b = lattice_constants.b
    c = lattice_constants.c
    α = lattice_constants.α
    β = lattice_constants.β
    γ = lattice_constants.γ

    # Compute dot products
    a_dot_a = a^2
    b_dot_b = b^2
    c_dot_c = c^2

    a_dot_b = dot(basis_a, basis_b)
    b_dot_c = dot(basis_b, basis_c)
    c_dot_a = dot(basis_c, basis_a)

    # Initialize monoclinic basis vectors
    m_basis_a = nothing
    m_basis_b = nothing
    m_basis_c = nothing

    # --- Attempt to convert the triclinic unit cell to a body-centered monoclinic unit cell

    if (b_dot_b ≈ c_dot_c) && (abs(a_dot_b) ≈ abs(c_dot_a) ≉ 0.5 * a_dot_a)
        # Case: `a` is monoclinic basis vector

        # Compute m_basis_a
        m_basis_a = basis_a

        # Compute m_basis_b and m_basis_c
        if a_dot_b ≈ c_dot_a
            # Case: coefficients of `m_basis_a` and `m_basis_c` are the same in `basis_b`
            #       and `basis_c`, so the sign of `m_basis_b` must be opposite in `basis_b`
            #       and `basis_c`

            # Compute m_basis_b
            m_basis_b = basis_b - basis_c

            # Compute m_basis_c
            if 2 * a_dot_b < a_dot_a
                m_basis_c = basis_b + basis_c - basis_a
            else
                m_basis_c = -(basis_b + basis_c - basis_a)
            end
        else
            # Case: coefficients of `m_basis_a` and `m_basis_c` are opposite in `basis_b`
            #       and `basis_c`, so the sign of `m_basis_b` must be the same in `basis_b`
            #       and `basis_c`

            # Compute m_basis_b
            m_basis_b = basis_b + basis_c

            # Compute m_basis_c
            if 2 * a_dot_b < a_dot_a
                m_basis_c = basis_b - basis_c - basis_a
            else
                m_basis_c = -(basis_b - basis_c - basis_a)
            end
        end

    elseif abs(b_dot_c) ≈ abs(a_dot_b) ≉ 0.5 * b_dot_b
        # Case: `b` is monoclinic basis vector

        # Compute m_basis_a
        m_basis_a = basis_b

        # Compute m_basis_b and m_basis_c
        if b_dot_c ≈ a_dot_b
            # Case: coefficients of `m_basis_a` and `m_basis_c` are the same in `basis_c`
            #       and `basis_a`, so the sign of `m_basis_b` must be opposite in `basis_c`
            #       and `basis_a`

            # Compute m_basis_b
            m_basis_b = basis_c - basis_a

            # Compute m_basis_c
            if 2 * b_dot_c < b_dot_b
                m_basis_c = basis_c + basis_a - basis_b
            else
                m_basis_c = -(basis_c + basis_a - basis_b)
            end
        else
            # Case: coefficients of `m_basis_a` and `m_basis_c` are opposite in `basis_c`
            #       and `basis_a`, so the sign of `m_basis_b` must be the same in `basis_c`
            #       and `basis_a`

            # Compute m_basis_b
            m_basis_b = basis_c + basis_a

            # Compute m_basis_c
            if 2 * b_dot_c < b_dot_b
                m_basis_c = basis_c - basis_a - basis_b
            else
                m_basis_c = -(basis_c - basis_a - basis_b)
            end
        end

    elseif abs(c_dot_a) ≈ abs(b_dot_c) ≉ 0.5 * c_dot_c
        # Case: `c` is monoclinic basis vector

        # Compute m_basis_a
        m_basis_a = basis_c

        # Compute m_basis_b and m_basis_c
        if c_dot_a ≈ b_dot_c
            # Case: coefficients of `m_basis_a` and `m_basis_c` are the same in `basis_a`
            #       and `basis_b`, so the sign of `m_basis_b` must be opposite in `basis_a`
            #       and `basis_b`

            # Compute m_basis_b
            m_basis_b = basis_a - basis_b

            # Compute m_basis_c
            if 2 * c_dot_a < c_dot_c
                m_basis_c = basis_a + basis_b - basis_c
            else
                m_basis_c = -(basis_a + basis_b - basis_c)
            end
        else
            # Case: coefficients of `m_basis_a` and `m_basis_c` are opposite in `basis_a`
            #       and `basis_b`, so the sign of `m_basis_b` must be the same in `basis_a`
            #       and `basis_b`

            # Compute m_basis_b
            m_basis_b = basis_a + basis_b

            # Compute m_basis_c
            if 2 * c_dot_a < c_dot_c
                m_basis_c = basis_a - basis_b - basis_c
            else
                m_basis_c = -(basis_a - basis_b - basis_c)
            end
        end
    end

    # --- Construct return value

    if isnothing(m_basis_a) || isnothing(m_basis_b) || isnothing(m_basis_c)
        throw(
            ErrorException(
                "The triclinic basis vectors defined by `lattice_constants` do not " *
                "satisfy conditions for case 4a.",
            ),
        )
    end

    return m_basis_a, m_basis_b, m_basis_c
end

function convert_to_mI_case_4b(lattice_constants::TriclinicLatticeConstants)
    # Notes
    # =====
    # - Case
    #   - triclinic unit cell basis contains 1 monoclinic unit cell basis vector and
    #     2 body-centered lattice vectors
    #
    #   - triclinic basis does not contain the unique monoclinic symmetry direction
    #     m_basis_b
    #
    #   - |dot(basis_a, basis_b) + dot(basis_a, basis_c)| = a^2
    #     and |dot(basis_a, basis_b)| != |dot(basis_a, basis_c)|
    #
    # - This method adopts the same variable conventions as convert_to_mI().

    # --- Preparations

    # Get basis for triclinic unit cell
    basis_a, basis_b, basis_c = basis(lattice_constants)

    # Get lattice constants for triclinic unit cell
    a = lattice_constants.a
    b = lattice_constants.b
    c = lattice_constants.c
    α = lattice_constants.α
    β = lattice_constants.β
    γ = lattice_constants.γ

    # Compute dot products
    a_dot_a = a^2
    b_dot_b = b^2
    c_dot_c = c^2

    a_dot_b = dot(basis_a, basis_b)
    b_dot_c = dot(basis_b, basis_c)
    c_dot_a = dot(basis_c, basis_a)

    # Initialize monoclinic basis vectors
    m_basis_a = nothing
    m_basis_b = nothing
    m_basis_c = nothing

    # --- Attempt to convert the triclinic unit cell to a body-centered monoclinic unit cell

    if (abs(a_dot_b + c_dot_a) ≈ a_dot_a) && (abs(a_dot_b) ≉ abs(c_dot_a))
        # --- Case: `a` is monoclinic basis vector, coefficients of `m_a` have same
        #     sign in the 2 body-centered lattice vectors, coefficients of `m_c` have
        #     opposite signs in the 2 body-centered lattice vectors

        m_basis_a = basis_a

        if a_dot_b + c_dot_a < 0
            m_basis_b = basis_b + basis_c + basis_a
        else
            m_basis_b = basis_b + basis_c - basis_a
        end

        if a_dot_b < c_dot_a
            m_basis_c = basis_b - basis_c
        else
            m_basis_c = -basis_b + basis_c
        end

    elseif (abs(b_dot_c + a_dot_b) ≈ b_dot_b) && (abs(b_dot_c) ≉ abs(a_dot_b))
        # --- Case: `b` is monoclinic basis vector, coefficients of `m_b` have same
        #     sign in the 2 body-centered lattice vectors, coefficients of `m_a` have
        #     opposite signs in the 2 body-centered lattice vectors

        m_basis_a = basis_b

        if b_dot_c + a_dot_b < 0
            m_basis_b = basis_c + basis_a + basis_b
        else
            m_basis_b = basis_c + basis_a - basis_b
        end

        if b_dot_c < a_dot_b
            m_basis_c = basis_c - basis_a
        else
            m_basis_c = -basis_c + basis_a
        end

    elseif (abs(c_dot_a + b_dot_c) ≈ c_dot_c) && (abs(c_dot_a) ≉ abs(b_dot_c))
        # --- Case: `c` is monoclinic basis vector, coefficients of `m_c` have same
        #     sign in the 2 body-centered lattice vectors, coefficients of `m_b` have
        #     opposite signs in the 2 body-centered lattice vectors

        m_basis_a = basis_c

        if c_dot_a + b_dot_c < 0
            m_basis_b = basis_a + basis_b + basis_c
        else
            m_basis_b = basis_a + basis_b - basis_c
        end

        if c_dot_a < b_dot_c
            m_basis_c = basis_a - basis_b
        else
            m_basis_c = -basis_a + basis_b
        end
    end

    # --- Construct return value

    if isnothing(m_basis_a) || isnothing(m_basis_b) || isnothing(m_basis_c)
        throw(
            ErrorException(
                "The triclinic basis vectors defined by `lattice_constants` do not " *
                "satisfy conditions for case 4b.",
            ),
        )
    end

    return m_basis_a, m_basis_b, m_basis_c
end

function convert_to_mI_case_4c(lattice_constants::TriclinicLatticeConstants)
    # Notes
    # =====
    # - Case
    #   - triclinic unit cell basis contains 1 monoclinic unit cell basis vector and
    #     2 body-centered lattice vectors
    #
    #   - triclinic basis does not contain the unique monoclinic symmetry direction
    #     m_basis_b
    #
    #   - |dot(basis_a, basis_b) - dot(basis_a, basis_c)| = a^2
    #     and |dot(basis_a, basis_b)| != |dot(basis_a, basis_c)|
    #
    #
    # - This method adopts the same variable conventions as convert_to_mI().

    # --- Preparations

    # Get basis for triclinic unit cell
    basis_a, basis_b, basis_c = basis(lattice_constants)

    # Get lattice constants for triclinic unit cell
    a = lattice_constants.a
    b = lattice_constants.b
    c = lattice_constants.c
    α = lattice_constants.α
    β = lattice_constants.β
    γ = lattice_constants.γ

    # Compute dot products
    a_dot_a = a^2
    b_dot_b = b^2
    c_dot_c = c^2

    a_dot_b = dot(basis_a, basis_b)
    b_dot_c = dot(basis_b, basis_c)
    c_dot_a = dot(basis_c, basis_a)

    # Initialize monoclinic basis vectors
    m_basis_a = nothing
    m_basis_b = nothing
    m_basis_c = nothing

    # --- Attempt to convert the triclinic unit cell to a body-centered monoclinic unit cell

    if (abs(a_dot_b - c_dot_a) ≈ a_dot_a) && (abs(a_dot_b) ≉ abs(c_dot_a))
        # --- Case: `a` is monoclinic basis vector, coefficients of `m_a` have opposite
        #     signs in the 2 body-centered lattice vectors, coefficients of `m_c` have
        #     the same sign in the 2 body-centered lattice vectors

        m_basis_a = basis_a

        if a_dot_b - c_dot_a < 0
            m_basis_b = basis_b - basis_c + basis_a
        else
            m_basis_b = -basis_b + basis_c + basis_a
        end

        if a_dot_b < -c_dot_a
            m_basis_c = basis_b + basis_c
        else
            m_basis_c = -(basis_b + basis_c)
        end

    elseif (abs(b_dot_c - a_dot_b) ≈ b_dot_b) && (abs(b_dot_c) ≉ abs(a_dot_b))
        # --- Case: `b` is monoclinic basis vector, coefficients of `m_b` have opposite
        #     signs in the 2 body-centered lattice vectors, coefficients of `m_a` have
        #     the same sign in the 2 body-centered lattice vectors

        m_basis_a = basis_b

        if b_dot_c - a_dot_b < 0
            m_basis_b = basis_c - basis_a + basis_b
        else
            m_basis_b = -basis_c + basis_a + basis_b
        end

        if b_dot_c < -a_dot_b
            m_basis_c = basis_c + basis_a
        else
            m_basis_c = -(basis_c + basis_a)
        end

    elseif (abs(c_dot_a - b_dot_c) ≈ c_dot_c) && (abs(c_dot_a) ≉ abs(b_dot_c))
        # --- Case: `c` is monoclinic basis vector, coefficients of `m_c` have opposite
        #     signs in the 2 body-centered lattice vectors, coefficients of `m_b` have
        #     the same sign in the 2 body-centered lattice vectors

        m_basis_a = basis_c

        if c_dot_a - b_dot_c < 0
            m_basis_b = basis_a - basis_b + basis_c
        else
            m_basis_b = -basis_a + basis_b + basis_c
        end

        if c_dot_a < -b_dot_c
            m_basis_c = basis_a + basis_b
        else
            m_basis_c = -(basis_a + basis_b)
        end
    end

    # --- Construct return value

    if isnothing(m_basis_a) || isnothing(m_basis_b) || isnothing(m_basis_c)
        throw(
            ErrorException(
                "The triclinic basis vectors defined by `lattice_constants` do not " *
                "satisfy conditions for case 4c.",
            ),
        )
    end

    return m_basis_a, m_basis_b, m_basis_c
end

function convert_to_mI_case_5(lattice_constants::TriclinicLatticeConstants)
    # Notes
    # =====
    # - Case
    #   - triclinic unit cell basis contains 3 body-centered lattice vectors
    #
    # - This method adopts the same variable conventions as convert_to_mI().

    # --- Preparations

    # Get basis for triclinic unit cell
    basis_a, basis_b, basis_c = basis(lattice_constants)

    # Get lattice constants for triclinic unit cell
    a = lattice_constants.a
    b = lattice_constants.b
    c = lattice_constants.c
    α = lattice_constants.α
    β = lattice_constants.β
    γ = lattice_constants.γ

    # Compute dot products
    a_dot_a = a^2
    b_dot_b = b^2
    c_dot_c = c^2

    a_dot_b = dot(basis_a, basis_b)
    b_dot_c = dot(basis_b, basis_c)
    c_dot_a = dot(basis_c, basis_a)

    # Initialize monoclinic basis vectors
    m_basis_a = nothing
    m_basis_b = nothing
    m_basis_c = nothing

    m_basis_d_1 = nothing
    m_basis_d_2 = nothing

    # --- Attempt to convert the triclinic unit cell to a body-centered monoclinic unit cell

    if a_dot_a ≈ b_dot_b && abs(c_dot_a) ≉ abs(b_dot_c)
        # --- Case: `a` and `b` are the basis vector with the same length

        if 2 * abs(c_dot_a + b_dot_c) ≈ a_dot_a + b_dot_b + 2 * a_dot_b
            m_basis_b = basis_a + basis_b
            m_basis_d_1 = basis_a - basis_b

        elseif 2 * abs(c_dot_a - b_dot_c) ≈ a_dot_a + b_dot_b - 2 * a_dot_b
            m_basis_b = basis_a - basis_b
            m_basis_d_1 = basis_a + basis_b
        end

        m_basis_d_2 =
            2 * (basis_c - dot(basis_c, m_basis_b) / dot(m_basis_b, m_basis_b) * m_basis_b)

    elseif b_dot_b ≈ c_dot_c && abs(a_dot_b) ≉ abs(c_dot_a)
        # --- Case: `b` and `c` are the basis vector with the same length

        if 2 * abs(a_dot_b + c_dot_a) ≈ b_dot_b + c_dot_c + 2 * b_dot_c
            m_basis_b = basis_b + basis_c
            m_basis_d_1 = basis_b - basis_c

        elseif 2 * abs(a_dot_b - c_dot_a) ≈ b_dot_b + c_dot_c - 2 * b_dot_c
            m_basis_b = basis_b - basis_c
            m_basis_d_1 = basis_b + basis_c
        end

        m_basis_d_2 =
            2 * (basis_a - dot(basis_a, m_basis_b) / dot(m_basis_b, m_basis_b) * m_basis_b)

    elseif c_dot_c ≈ a_dot_a && abs(b_dot_c) ≉ abs(a_dot_b)
        # --- Case: `c` and `a` are the basis vector with the same length

        if 2 * abs(b_dot_c + a_dot_b) ≈ c_dot_c + a_dot_a + 2 * c_dot_a
            m_basis_b = basis_c + basis_a
            m_basis_d_1 = basis_c - basis_a

        elseif 2 * abs(b_dot_c - a_dot_b) ≈ c_dot_c + a_dot_a - 2 * c_dot_a
            m_basis_b = basis_c - basis_a
            m_basis_d_1 = basis_c + basis_a
        end

        m_basis_d_2 =
            2 * (basis_b - dot(basis_b, m_basis_b) / dot(m_basis_b, m_basis_b) * m_basis_b)
    end

    if !isnothing(m_basis_d_1) && !isnothing(m_basis_d_2)
        m_basis_a = 0.5 * (m_basis_d_1 + m_basis_d_2)

        if dot(m_basis_d_1, m_basis_d_1) < dot(m_basis_d_2, m_basis_d_2)
            m_basis_c = 0.5 * (m_basis_d_1 - m_basis_d_2)
        else
            m_basis_c = 0.5 * (-m_basis_d_1 + m_basis_d_2)
        end
    end

    # --- Construct return value

    if isnothing(m_basis_a) || isnothing(m_basis_b) || isnothing(m_basis_c)
        throw(
            ErrorException(
                "The triclinic basis vectors defined by `lattice_constants` do not " *
                "include three body-centered lattice vectors.",
            ),
        )
    end

    return m_basis_a, m_basis_b, m_basis_c
end

"""
    convert_to_mS(
        lattice_constants::TriclinicLatticeConstants
    ) -> MonoclinicLatticeConstants

Attempt to convert the triclinic unit cell defined by `lattice_constants` to an equivalent
base-centered monoclinic unit cell.

Return values
=============
- lattice constants for the equivalent base-centered monoclinic unit cell if one exists;
  `nothing` otherwise

  !!! warn

      Returned lattice constants are _not_ guaranteed to be standardized.

Exceptions
==========
Throws an `ErrorException` if the triclinic unit cell defined by `lattice_constants` is not
equivalent to a base-centered monoclinic unit cell.
"""
function convert_to_mS(lattice_constants::TriclinicLatticeConstants)
    # Notes
    # =====
    # This method adopts the following variable conventions.
    #
    # - Unless otherwise noted, lattice constants and basis vectors refer to the triclinic
    #   (not monoclinic) unit cell.
    #
    # - Lattice constants and basis vectors for the monoclinic unit cell are indicated by
    #   the "m_" prefix.

    # --- Preparations

    # Get basis for triclinic unit cell
    basis_a, basis_b, basis_c = basis(lattice_constants)

    # Get lattice constants for triclinic unit cell
    a = lattice_constants.a
    b = lattice_constants.b
    c = lattice_constants.c
    α = lattice_constants.α
    β = lattice_constants.β
    γ = lattice_constants.γ

    # Compute dot products
    a_dot_a = a^2
    b_dot_b = b^2
    c_dot_c = c^2

    a_dot_b = dot(basis_a, basis_b)
    b_dot_c = dot(basis_b, basis_c)
    c_dot_a = dot(basis_c, basis_a)

    # Initialize monoclinic lattice constants
    m_a = nothing
    m_b = nothing
    m_c = nothing
    m_β = nothing

    # Initialize monoclinic basis vectors
    m_basis_a = nothing
    m_basis_b = nothing
    m_basis_c = nothing

    # --- Attempt to convert the triclinic unit cell to a base-centered monoclinic unit cell

    # ------ Case: triclinic unit cell basis contains 2 monoclinic unit cell basis vectors
    #        and 1 base-centered lattice vector

    # Case: triclinic basis contains m_basis_a and m_basis_c
    try
        # Compute monoclinic basis vectors
        m_basis_a, m_basis_b, m_basis_c = convert_to_mS_case_1(lattice_constants)
        return convert_to_mS_basis_to_lattice_constants(m_basis_a, m_basis_b, m_basis_c)

    catch error
        if !(error isa ErrorException) || (
            error.msg !=
            "The triclinic basis vectors defined by `lattice_constants` do not " *
            "include the mononclinic basis vectors that are not the unique " *
            "monoclinic symmetry direction and one base-centered lattice vector."
        )
            rethrow(error)
        end
    end

    # ------ Case: triclinic unit cell basis contains 1 monoclinic unit cell basis vector
    #        and 2 base-centered lattice vectors

    # Case: triclinic basis contains the m_basis_a
    try
        # Compute monoclinic basis vectors
        m_basis_a, m_basis_b, m_basis_c = convert_to_mS_case_2(lattice_constants)
        return convert_to_mS_basis_to_lattice_constants(m_basis_a, m_basis_b, m_basis_c)

    catch error
        if !(error isa ErrorException) || (
            error.msg !=
            "The triclinic basis vectors defined by `lattice_constants` do not " *
            "include the mononclinic basis vector that is the intersection of the " *
            "B-face and C-face and two base-centered lattice vectors."
        )
            rethrow(error)
        end
    end

    # --- Throw exception if unit cell is not equivalent to a base-centered monoclinic unit

    if isnothing(m_basis_a) || isnothing(m_basis_b) || isnothing(m_basis_c)
        throw(
            ErrorException(
                "The triclinic unit cell defined by `lattice_constants` is not " *
                "equivalent to a base-centered monoclinic unit cell.",
            ),
        )
    end
end

function convert_to_mS_basis_to_lattice_constants(
    m_basis_a::Vector{<:Real}, m_basis_b::Vector{<:Real}, m_basis_c::Vector{<:Real}
)
    # Compute monoclinic lattice constants
    m_a = norm(m_basis_a)
    m_b = norm(m_basis_b)
    m_c = norm(m_basis_c)
    m_β = angle(m_basis_a, m_basis_c)

    # Check IUCr conventions
    if m_β < π / 2
        throw(ErrorException("m_β > 0"))
    end

    return MonoclinicLatticeConstants(m_a, m_b, m_c, m_β)
end

#=
function convert_to_mS_OLD(lattice_constants::TriclinicLatticeConstants)
    # Notes
    # =====
    # This method adopts the following variable conventions.
    #
    # - Unless otherwise noted, lattice constants and basis vectors refer to the triclinic
    #   (not monoclinic) unit cell.
    #
    # - Lattice constants and basis vectors for the monoclinic unit cell are indicated by
    #   the "m_" prefix.

    # --- Preparations

    # Get basis vectors of triclinic unit cell
    basis_a, basis_b, basis_c = basis(lattice_constants)

    # Get lattice constants for triclinic unit cell
    a = norm(basis_a)
    b = norm(basis_b)
    c = norm(basis_c)
    α = dot(basis_b, basis_c) / b / c
    β = dot(basis_c, basis_a) / c / a
    γ = dot(basis_a, basis_b) / a / b

    # Compute dot products
    a_dot_a = a^2
    b_dot_b = b^2
    c_dot_c = c^2

    a_dot_b = dot(basis_a, basis_b)
    b_dot_c = dot(basis_b, basis_c)
    c_dot_a = dot(basis_c, basis_a)

    # Initialize monoclinic lattice constants
    m_a = nothing
    m_b = nothing
    m_c = nothing
    m_β = nothing

    # --- Attempt to convert the triclinic unit cell to a base-centered monoclinic unit cell

    # ------ Case: triclinic unit cell basis contains 2 monoclinic unit cell basis vectors
    #        and 1 base-centered vector

    # --------- Case: triclinic basis contains unique monoclinic symmetry direction

    if abs(a_dot_b) < a * b * COS_APPROX_ZERO
        # Case: `a` and `b` are basis vectors of the monoclinic unit cell

        if abs(b_dot_c) ≈ 0.5 * b_dot_b
            # Case: `b` is the unique monoclinic symmetry direction, `a` is basis vector
            #       on B-face of the monoclinic unit cell

            # Compute monoclinic lattice constants
            m_a = a
            m_b = b
            m_c = sqrt(4 * c_dot_c - b_dot_b)
            m_β = acos(2 * c_dot_a / m_a / m_c)

        elseif abs(c_dot_a) ≈ 0.5 * a_dot_a
            # Case: `a` is the unique monoclinic symmetry direction, `b` is basis vector
            #       on B-face of the monoclinic unit cell

            # Compute monoclinic lattice constants
            m_a = b
            m_b = a
            m_c = sqrt(4 * c_dot_c - a_dot_a)
            m_β = acos(2 * b_dot_c / m_a / m_c)
        end

    elseif abs(b_dot_c) < b * c * COS_APPROX_ZERO
        # Case: `b` and `c` are basis vectors of the monoclinic unit cell

        if abs(c_dot_a) ≈ 0.5 * c_dot_c
            # Case: `c` is the unique monoclinic symmetry direction, `b` is basis vector
            #       on B-face of the monoclinic unit cell

            # Compute monoclinic lattice constants
            m_a = b
            m_b = c
            m_c = sqrt(4 * a_dot_a - c_dot_c)
            m_β = acos(2 * a_dot_b / m_a / m_c)

        elseif abs(a_dot_b) ≈ 0.5 * b_dot_b
            # Case: `b` is the unique monoclinic symmetry direction, `c` is basis vector
            #       on B-face of the monoclinic unit cell

            # Compute monoclinic lattice constants
            m_a = c
            m_b = b
            m_c = sqrt(4 * a_dot_a - b_dot_b)
            m_β = acos(2 * c_dot_a / m_a / m_c)
        end

    elseif abs(c_dot_a) < c * a * COS_APPROX_ZERO
        # Case: `c` and `a` are basis vectors of the monoclinic unit cell

        if abs(a_dot_b) ≈ 0.5 * a_dot_a
            # Case: `a` is the unique monoclinic symmetry direction, `c` is basis vector
            #       on B-face of the monoclinic unit cell

            # Compute monoclinic lattice constants
            m_a = c
            m_b = a
            m_c = sqrt(4 * b_dot_b - a_dot_a)
            m_β = acos(2 * b_dot_c / m_a / m_c)

        elseif abs(b_dot_c) ≈ 0.5 * c_dot_c
            # Case: `c` is the unique monoclinic symmetry direction, `a` is basis vector
            #       on B-face of the monoclinic unit cell

            # Compute monoclinic lattice constants
            m_a = a
            m_b = c
            m_c = sqrt(4 * b_dot_b - c_dot_c)
            m_β = acos(2 * a_dot_b / m_a / m_c)
        end
    end

    # --------- Case: triclinic basis does not contain unique monoclinic symmetry direction

    if (abs(a_dot_b) ≈ 2 * abs(c_dot_a) || abs(a_dot_b) ≈ 2 * abs(b_dot_c))
        # Case: `a` and `b` are basis vectors of the monoclinic unit cell

        if abs(a_dot_b) ≈ 2 * abs(b_dot_c)
            # Case: `a` is basis vector on C-face of the monoclinic unit cell

            # Compute monoclinic lattice constants
            m_a = a
            m_b = sqrt(4 * c_dot_c - a_dot_a)
            m_c = b
            m_β = acos(a_dot_b / m_a / m_c)

        else  # abs(a_dot_b) ≈ 2 * abs(c_dot_a)
            # Case: `b` is basis vector on C-face of the monoclinic unit cell

            # Compute monoclinic lattice constants
            m_a = b
            m_b = sqrt(4 * c_dot_c - b_dot_b)
            m_c = a
            m_β = acos(a_dot_b / m_a / m_c)
        end

    elseif (abs(b_dot_c) ≈ 2 * abs(a_dot_b) || abs(b_dot_c) ≈ 2 * abs(c_dot_a))
        # Case: `b` and `c` are basis vectors of the monoclinic unit cell

        if abs(b_dot_c) ≈ 2 * abs(c_dot_a)
            # Case: `b` is basis vector on C-face of the monoclinic unit cell

            # Compute monoclinic lattice constants
            m_a = b
            m_b = sqrt(4 * a_dot_a - b_dot_b)
            m_c = c
            m_β = acos(b_dot_c / m_a / m_c)

        else  # abs(b_dot_c) ≈ 2 * abs(a_dot_b)
            # Case: `c` is basis vector on C-face of the monoclinic unit cell

            # Compute monoclinic lattice constants
            m_a = c
            m_b = sqrt(4 * a_dot_a - c_dot_c)
            m_c = b
            m_β = acos(b_dot_c / m_a / m_c)
        end

    elseif (abs(c_dot_a) ≈ 2 * abs(b_dot_c) || abs(c_dot_a) ≈ 2 * abs(a_dot_b))
        # Case: `c` and `a` are basis vectors of the monoclinic unit cell

        if abs(c_dot_a) ≈ 2 * abs(a_dot_b)
            # Case: `c` is basis vector on C-face of the monoclinic unit cell

            # Compute monoclinic lattice constants
            m_a = c
            m_b = sqrt(4 * b_dot_b - c_dot_c)
            m_c = a
            m_β = acos(c_dot_a / m_a / m_c)

        else  # abs(c_dot_a) ≈ 2 * abs(b_dot_c)
            # Case: `a` is basis vector on C-face of the monoclinic unit cell

            # Compute monoclinic lattice constants
            m_a = a
            m_b = sqrt(4 * b_dot_b - a_dot_a)
            m_c = c
            m_β = acos(c_dot_a / m_a / m_c)
        end
    end

    # ------ Case: triclinic unit cell basis contains 1 monoclinic unit cell basis vector
    #        2 base-centered vectors

    # --------- Case: triclinic basis contains the unique monoclinic symmetry direction

    # Impossible because basis is not linearly independent

    # --------- Case: triclinic basis does not contain unique monoclinic symmetry direction

    if (abs(a_dot_b) ≈ abs(c_dot_a) && b_dot_b ≈ c_dot_c)
        # Case: `a` is C-basis vector of the monoclinic unit cell

        # --- Compute monoclinic basis vectors

        m_basis_a = basis_b + basis_c
        m_basis_b = basis_b - basis_c
        m_basis_c = basis_a

        # Swap m_basis_a and m_basis_b if they are reversed
        if abs(dot(m_basis_b, m_basis_c)) >= a * norm(m_basis_b) * COS_APPROX_ZERO
            m_basis_tmp = m_basis_a
            m_basis_a = m_basis_b
            m_basis_b = m_basis_tmp
        end

        # --- Compute monoclinic lattice constants

        m_a = norm(m_basis_a)
        m_b = norm(m_basis_b)
        m_c = a
        m_β = angle(m_basis_a, m_basis_c)

    elseif (abs(b_dot_c) ≈ abs(a_dot_b) && c_dot_c ≈ a_dot_a)
        # Case: `b` is C-basis vector of the monoclinic unit cell

        # --- Compute monoclinic basis vectors

        m_basis_a = basis_c + basis_a
        m_basis_b = basis_c - basis_a
        m_basis_c = basis_b

        # Swap m_basis_a and m_basis_b if they are reversed
        if abs(dot(m_basis_b, m_basis_c)) >= b * norm(m_basis_b) * COS_APPROX_ZERO
            m_basis_tmp = m_basis_a
            m_basis_a = m_basis_b
            m_basis_b = m_basis_tmp
        end

        # --- Compute monoclinic lattice constants

        m_a = norm(m_basis_a)
        m_b = norm(m_basis_b)
        m_c = b
        m_β = angle(m_basis_a, m_basis_c)

    elseif (abs(c_dot_a) ≈ abs(b_dot_c) && a_dot_a ≈ b_dot_b)
        # Case: `c` is C-basis vector of the monoclinic unit cell

        # --- Compute monoclinic basis vectors

        m_basis_a = basis_a + basis_b
        m_basis_b = basis_a - basis_b
        m_basis_c = basis_c

        # Swap m_basis_a and m_basis_b if they are reversed
        if abs(dot(m_basis_b, m_basis_c)) >= c * norm(m_basis_b) * COS_APPROX_ZERO
            m_basis_tmp = m_basis_a
            m_basis_a = m_basis_b
            m_basis_b = m_basis_tmp
        end

        # --- Compute monoclinic lattice constants

        m_a = norm(m_basis_a)
        m_b = norm(m_basis_b)
        m_c = c
        m_β = angle(m_basis_a, m_basis_c)
    end

    # ------ Case: triclinic unit cell basis contains 3 base-centered vectors

    if (a_dot_a ≈ b_dot_b && abs(c_dot_a) ≉ abs(b_dot_c))
        # Case: `a` and `b` are base-centered basis vectors that combine only two
        #       monoclinic basis vectors

        # --- Compute monoclinic basis vectors

        m_basis_a = basis_a + basis_b
        m_basis_b = basis_a - basis_b

        # Swap m_basis_a and m_basis_a if they are reversed
        #
        # m_basis_b should satisfy the following equation
        #
        #   abs(dot(m_basis_b, basis_c)) = 0.5 * dot(m_basis_b, m_basis_b)
        if abs(dot(m_basis_b, basis_c)) ≉ 0.5 * dot(m_basis_b, m_basis_b)
            m_basis_tmp = m_basis_a
            m_basis_a = m_basis_b
            m_basis_b = m_basis_tmp
        end

        # Compute m_basis_c by (1) removing m_basis_b from basis_c and (2) adding
        # 0.5 * m_basis_a
        m_basis_c =
            basis_c - dot(basis_c, m_basis_b) / dot(m_basis_b, m_basis_b) * m_basis_b +
            0.5 * m_basis_a

        # --- Compute monoclinic lattice constants

        m_a = norm(m_basis_a)
        m_b = norm(m_basis_b)
        m_c = norm(m_basis_c)
        m_β = angle(m_basis_a, m_basis_c)

    elseif (b_dot_b ≈ c_dot_c && abs(a_dot_b) ≉ abs(c_dot_a))
        # Case: `b` and `c` are base-centered basis vectors that combine only two
        #       monoclinic basis vectors

        # --- Compute monoclinic basis vectors

        m_basis_a = basis_b + basis_c
        m_basis_b = basis_b - basis_c

        # Swap m_basis_a and m_basis_a if they are reversed
        #
        # m_basis_b should satisfy the following equation
        #
        #   abs(dot(m_basis_b, basis_a)) = 0.5 * dot(m_basis_b, m_basis_b)
        if abs(dot(m_basis_b, basis_a)) ≉ 0.5 * dot(m_basis_b, m_basis_b)
            m_basis_tmp = m_basis_a
            m_basis_a = m_basis_b
            m_basis_b = m_basis_tmp
        end

        # Compute m_basis_c by (1) removing m_basis_b from basis_a and (2) adding
        # 0.5 * m_basis_a
        m_basis_c =
            basis_a - dot(basis_a, m_basis_b) / dot(m_basis_b, m_basis_b) * m_basis_b +
            0.5 * m_basis_a

        # --- Compute monoclinic lattice constants

        m_a = norm(m_basis_a)
        m_b = norm(m_basis_b)
        m_c = norm(m_basis_c)
        m_β = angle(m_basis_a, m_basis_c)

    elseif (c_dot_c ≈ a_dot_a && abs(b_dot_c) ≉ abs(a_dot_b))
        # Case: `c` and `a` are base-centered basis vectors that combine only two
        #       monoclinic basis vectors

        # --- Compute monoclinic basis vectors

        m_basis_a = basis_c + basis_a
        m_basis_b = basis_c - basis_a

        m_basis_a = basis_c - basis_a
        m_basis_b = basis_c + basis_a

        # Swap m_basis_a and m_basis_a if they are reversed
        #
        # m_basis_b should satisfy the following equation
        #
        #   abs(dot(m_basis_b, basis_b)) = 0.5 * dot(m_basis_b, m_basis_b)
        if abs(dot(m_basis_b, basis_b)) ≉ 0.5 * dot(m_basis_b, m_basis_b)
            m_basis_tmp = m_basis_a
            m_basis_a = m_basis_b
            m_basis_b = m_basis_tmp
        end

        # Compute m_basis_c by (1) removing m_basis_b from basis_b and (2) adding
        # 0.5 * m_basis_a
        m_basis_c =
            basis_b - dot(basis_b, m_basis_b) / dot(m_basis_b, m_basis_b) * m_basis_b +
            0.5 * m_basis_a

        # --- Compute monoclinic lattice constants

        m_a = norm(m_basis_a)
        m_b = norm(m_basis_b)
        m_c = norm(m_basis_c)
        m_β = angle(m_basis_a, m_basis_c)
    end

    # --- Construct return value

    if isnothing(m_a) || isnothing(m_b) || isnothing(m_c) || isnothing(m_β)
        throw(
            ErrorException(
                "The triclinic unit cell defined by `lattice_constants` is not " *
                "equivalent to a base-centered monoclinic unit cell.",
            ),
        )
    end

    mI_lattice_constants, _ = standardize(
        MonoclinicLatticeConstants(m_a, m_b, m_c, m_β), BASE
    )
    mS_lattice_constants = convert_to_base_centering(mI_lattice_constants)

    return mS_lattice_constants
end
=#

function convert_to_mS_case_1(lattice_constants::TriclinicLatticeConstants)
    # --- Preparations

    # Initialize monoclinic basis vectors
    m_basis_a = nothing
    m_basis_b = nothing
    m_basis_c = nothing

    # --- Attempt to convert the triclinic unit cell to a base-centered monoclinic unit cell

    try
        # Compute monoclinic basis vectors
        m_basis_a, m_basis_b, m_basis_c = convert_to_mS_case_1a(lattice_constants)

        # Compute monoclinic lattice constants
        m_a = norm(m_basis_a)
        m_b = norm(m_basis_b)
        m_c = norm(m_basis_c)
        m_β = angle(m_basis_a, m_basis_c)

    catch error
        if !(error isa ErrorException) || (
            error.msg !=
            "The triclinic basis vectors defined by `lattice_constants` do not " *
            "satisfy conditions for case 1a."
        )
            rethrow(error)
        end
    end

    try
        # Compute monoclinic basis vectors
        m_basis_a, m_basis_b, m_basis_c = convert_to_mS_case_1b(lattice_constants)

        # Compute monoclinic lattice constants
        m_a = norm(m_basis_a)
        m_b = norm(m_basis_b)
        m_c = norm(m_basis_c)
        m_β = angle(m_basis_a, m_basis_c)

    catch error
        if !(error isa ErrorException) || (
            error.msg !=
            "The triclinic basis vectors defined by `lattice_constants` do not " *
            "satisfy conditions for case 1b."
        )
            rethrow(error)
        end
    end

    # --- Construct return value

    if isnothing(m_basis_a) || isnothing(m_basis_b) || isnothing(m_basis_c)
        throw(
            ErrorException(
                "The triclinic basis vectors defined by `lattice_constants` do not " *
                "include the mononclinic basis vectors that are not the unique " *
                "monoclinic symmetry direction and one base-centered lattice vector.",
            ),
        )
    end

    return m_basis_a, m_basis_b, m_basis_c
end

function convert_to_mS_case_1a(lattice_constants::TriclinicLatticeConstants)
    # Notes
    # =====
    # - Case
    #   - triclinic unit cell basis contains m_basis_a, m_basis_c, and 1 base-centered
    #     lattice vector
    #
    #   - base-centered vector does not include a contribution from m_basis_c
    #
    #   - triclinic basis contains the m_basis_a and m_basis_c
    #
    # - This method adopts the same variable conventions as convert_to_mS().

    # --- Attempt to convert the triclinic unit cell to a base-centered monoclinic unit cell

    for (basis_a, basis_b, basis_c) in permutations(basis(lattice_constants))
        # Initialize monoclinic basis vectors
        m_basis_a = nothing
        m_basis_b = nothing
        m_basis_c = nothing

        a_dot_a = dot(basis_a, basis_a)
        c_dot_c = dot(basis_c, basis_c)

        a_dot_b = dot(basis_a, basis_b)
        b_dot_c = dot(basis_b, basis_c)
        c_dot_a = dot(basis_c, basis_a)

        if 2 * abs(a_dot_b) ≈ a_dot_a && 2 * abs(b_dot_c) ≈ abs(c_dot_a)
            # --- Case: `a` and `c` are the monoclinic basis vectors in the triclinic basis
            #     and `a` is m_basis_a

            m_basis_a = basis_a

            if 2 * a_dot_b ≈ a_dot_a
                m_basis_b = 2 * basis_b - basis_a
            else
                m_basis_b = 2 * basis_b + basis_a
            end

            if c_dot_a < 0
                m_basis_c = basis_c
            else
                m_basis_c = -basis_c
            end
        end

        if !isnothing(m_basis_a) && !isnothing(m_basis_b) && !isnothing(m_basis_c)
            return m_basis_a, m_basis_b, m_basis_c
        end
    end

    # --- Failed to convert the triclinic unit cell to a base-centered monoclinic unit cell

    throw(
        ErrorException(
            "The triclinic basis vectors defined by `lattice_constants` do not " *
            "satisfy conditions for case 1a.",
        ),
    )
end

function convert_to_mS_case_1b(lattice_constants::TriclinicLatticeConstants)
    # Notes
    # =====
    # - Case
    #   - triclinic unit cell basis contains m_basis_a, m_basis_c, and 1 base-centered
    #     lattice vector
    #
    #   - base-centered vector does includes a contribution from m_basis_c
    #
    #   - triclinic basis contains the m_basis_a and m_basis_c
    #
    # - This method adopts the same variable conventions as convert_to_mS().

    # --- Attempt to convert the triclinic unit cell to a base-centered monoclinic unit cell

    for (basis_a, basis_b, basis_c) in permutations(basis(lattice_constants))
        # Initialize monoclinic basis vectors
        m_basis_a = nothing
        m_basis_b = nothing
        m_basis_c = nothing

        a_dot_a = dot(basis_a, basis_a)
        c_dot_c = dot(basis_c, basis_c)

        a_dot_b = dot(basis_a, basis_b)
        b_dot_c = dot(basis_b, basis_c)
        c_dot_a = dot(basis_c, basis_a)

        if abs(a_dot_b - 0.5 * a_dot_a) ≈ abs(c_dot_a)
            # --- Case: the sign of m_basis_a in `b` is positive

            # Compute m_basis_a
            m_basis_a = basis_a

            # Compute m_basis_b
            if 2 * abs(b_dot_c - c_dot_c) ≈ abs(c_dot_a)
                # Case: signs of the coefficients of m_basis_c in `b` and `c` are the same
                m_basis_b = 2 * basis_b - 2 * basis_c - basis_a

            elseif 2 * abs(b_dot_c + c_dot_c) ≈ abs(c_dot_a)
                # Case: signs of the coefficients of m_basis_c in `b` and `c` are opposite
                m_basis_b = 2 * basis_b + 2 * basis_c - basis_a
            end

            # Compute m_basis_c
            if c_dot_a < 0
                m_basis_c = basis_c
            else
                m_basis_c = -basis_c
            end

        elseif abs(a_dot_b + 0.5 * a_dot_a) ≈ abs(c_dot_a)
            # --- Case: the sign of m_basis_a in `b` is negative

            # Compute m_basis_a
            m_basis_a = basis_a

            # Compute m_basis_b
            if 2 * abs(b_dot_c - c_dot_c) ≈ abs(c_dot_a)
                # Case: signs of the coefficients of m_basis_c in `b` and `c` are the same
                m_basis_b = 2 * basis_b - 2 * basis_c + basis_a

            elseif 2 * abs(b_dot_c + c_dot_c) ≈ abs(c_dot_a)
                # Case: signs of the coefficients of m_basis_c in `b` and `c` are opposite
                m_basis_b = 2 * basis_b + 2 * basis_c + basis_a
            end

            # Compute m_basis_c
            if c_dot_a < 0
                m_basis_c = basis_c
            else
                m_basis_c = -basis_c
            end
        end

        if !isnothing(m_basis_a) && !isnothing(m_basis_b) && !isnothing(m_basis_c)
            return m_basis_a, m_basis_b, m_basis_c
        end
    end

    # --- Failed to convert the triclinic unit cell to a base-centered monoclinic unit cell

    throw(
        ErrorException(
            "The triclinic basis vectors defined by `lattice_constants` do not " *
            "satisfy conditions for case 1b.",
        ),
    )
end

function convert_to_mS_case_2(lattice_constants::TriclinicLatticeConstants)
    # --- Preparations

    # Initialize monoclinic basis vectors
    m_basis_a = nothing
    m_basis_b = nothing
    m_basis_c = nothing

    # --- Attempt to convert the triclinic unit cell to a base-centered monoclinic unit cell

    try
        # Compute monoclinic basis vectors
        m_basis_a, m_basis_b, m_basis_c = convert_to_mS_case_2a(lattice_constants)

        # Compute monoclinic lattice constants
        m_a = norm(m_basis_a)
        m_b = norm(m_basis_b)
        m_c = norm(m_basis_c)
        m_β = angle(m_basis_a, m_basis_c)

    catch error
        if !(error isa ErrorException) || (
            error.msg !=
            "The triclinic basis vectors defined by `lattice_constants` do not " *
            "satisfy conditions for case 2a."
        )
            rethrow(error)
        end
    end

    try
        # Compute monoclinic basis vectors
        m_basis_a, m_basis_b, m_basis_c = convert_to_mS_case_2b(lattice_constants)

        # Compute monoclinic lattice constants
        m_a = norm(m_basis_a)
        m_b = norm(m_basis_b)
        m_c = norm(m_basis_c)
        m_β = angle(m_basis_a, m_basis_c)

    catch error
        if !(error isa ErrorException) || (
            error.msg !=
            "The triclinic basis vectors defined by `lattice_constants` do not " *
            "satisfy conditions for case 2b."
        )
            rethrow(error)
        end
    end

    # --- Construct return value

    if isnothing(m_basis_a) || isnothing(m_basis_b) || isnothing(m_basis_c)
        throw(
            ErrorException(
                "The triclinic basis vectors defined by `lattice_constants` do not " *
                "include the mononclinic basis vector that is the intersection of the " *
                "B-face and C-face and two base-centered lattice vectors.",
            ),
        )
    end

    return m_basis_a, m_basis_b, m_basis_c
end

function convert_to_mS_case_2a(lattice_constants::TriclinicLatticeConstants)
    # Notes
    # =====
    # - Case
    #   - triclinic unit cell basis contains m_basis_a and 2 base-centered lattice vectors
    #
    #   - relative signs of coefficients of m_basis_a and m_basis_c in basis_b and basis_c
    #     are the same
    #
    # - This method adopts the same variable conventions as convert_to_mS().

    # --- Attempt to convert the triclinic unit cell to a base-centered monoclinic unit cell

    for (basis_a, basis_b, basis_c) in permutations(basis(lattice_constants))
        # Initialize monoclinic basis vectors
        m_basis_a = nothing
        m_basis_b = nothing
        m_basis_c = nothing

        a_dot_a = dot(basis_a, basis_a)
        b_dot_b = dot(basis_b, basis_b)

        a_dot_b = dot(basis_a, basis_b)
        b_dot_c = dot(basis_b, basis_c)
        c_dot_a = dot(basis_c, basis_a)

        if 2 * abs(a_dot_b) ≈ a_dot_a
            # Compute m_basis_a
            m_basis_a = basis_a

            # Compute m_basis_b
            if 2 * a_dot_b ≈ a_dot_a
                m_basis_b = 2 * basis_b - basis_a
            elseif 2 * a_dot_b ≈ -a_dot_a
                m_basis_b = 2 * basis_b + basis_a
            end

            # Compute m_basis_c
            if 2 * abs(b_dot_b + b_dot_c + a_dot_b) ≈ abs(a_dot_b - c_dot_a)
                # Case:
                #
                # (1) relative signs of m_basis_a and m_basis_c are opposite in
                #     basis_b and basis_c
                #
                # (2) signs of m_basis_a are the same in basis_b and basis_c
                #
                # (3) sign of m_basis_a in basis_b and basis_c is negative

                if a_dot_b + c_dot_a < -a_dot_a
                    m_basis_c = basis_b + basis_c + basis_a
                else
                    m_basis_c = -basis_b - basis_c - basis_a
                end

            elseif 2 * abs(b_dot_b + b_dot_c - a_dot_b) ≈ abs(a_dot_b - c_dot_a)
                # Case:
                #
                # (1) relative signs of m_basis_a and m_basis_c are opposite in
                #     basis_b and basis_c
                #
                # (2) signs of m_basis_a are the same in basis_b and basis_c
                #
                # (3) sign of m_basis_a in basis_b and basis_c is positive

                if a_dot_b + c_dot_a < a_dot_a
                    # Case: signs of m_basis_a are the same in basis_b and basis_c
                    m_basis_c = basis_b + basis_c - basis_a
                else
                    m_basis_c = -basis_b - basis_c + basis_a
                end

            elseif 2 * abs(b_dot_b - b_dot_c + a_dot_b) ≈ abs(a_dot_b + c_dot_a)
                # Case:
                #
                # (1) relative signs of m_basis_a and m_basis_c are opposite in
                #     basis_b and basis_c
                #
                # (2) signs of m_basis_a are opposite in basis_b and basis_c
                #
                # (3) sign of m_basis_a in basis_b and basis_c is negative

                if a_dot_b - c_dot_a < -a_dot_a
                    m_basis_c = basis_b - basis_c + basis_a
                else
                    m_basis_c = -basis_b + basis_c - basis_a
                end

            elseif 2 * abs(b_dot_b - b_dot_c - a_dot_b) ≈ abs(a_dot_b + c_dot_a)
                # Case:
                #
                # (1) relative signs of m_basis_a and m_basis_c are opposite in
                #     basis_b and basis_c
                #
                # (2) signs of m_basis_a are opposite in basis_b and basis_c
                #
                # (3) sign of m_basis_a in basis_b and basis_c is positive

                if a_dot_b - c_dot_a < a_dot_a
                    # Case: signs of m_basis_a are the same in basis_b and basis_c
                    m_basis_c = basis_b - basis_c - basis_a
                else
                    m_basis_c = -basis_b + basis_c + basis_a
                end
            end
        end

        if !isnothing(m_basis_a) && !isnothing(m_basis_b) && !isnothing(m_basis_c)
            return m_basis_a, m_basis_b, m_basis_c
        end
    end

    # --- Failed to convert the triclinic unit cell to a base-centered monoclinic unit cell

    throw(
        ErrorException(
            "The triclinic basis vectors defined by `lattice_constants` do not " *
            "satisfy conditions for case 2a.",
        ),
    )
end

function convert_to_mS_case_2b(lattice_constants::TriclinicLatticeConstants)
    # Notes
    # =====
    # - Case
    #   - triclinic unit cell basis contains m_basis_a and 2 base-centered lattice vectors
    #
    #   - relative signs of coefficients of m_basis_a and m_basis_c in basis_b and basis_c
    #     are opposite
    #
    # - This method adopts the same variable conventions as convert_to_mS().

    # --- Attempt to convert the triclinic unit cell to a base-centered monoclinic unit cell

    for (basis_a, basis_b, basis_c) in permutations(basis(lattice_constants))
        # Initialize monoclinic basis vectors
        m_basis_a = nothing
        m_basis_b = nothing
        m_basis_c = nothing

        a_dot_a = dot(basis_a, basis_a)
        b_dot_b = dot(basis_b, basis_b)

        a_dot_b = dot(basis_a, basis_b)
        b_dot_c = dot(basis_b, basis_c)
        c_dot_a = dot(basis_c, basis_a)

        if 2 * abs(a_dot_b) ≈ a_dot_a
            # Compute m_basis_a
            m_basis_a = basis_a

            # Compute m_basis_b
            if 2 * a_dot_b ≈ a_dot_a
                m_basis_b = 2 * basis_b - basis_a
            elseif 2 * a_dot_b ≈ -a_dot_a
                m_basis_b = 2 * basis_b + basis_a
            end

            # Compute m_basis_c
            if 2 * abs(b_dot_c - b_dot_b) ≈ abs(c_dot_a - a_dot_b)
                # Case:
                #
                # (1) relative signs of m_basis_a and m_basis_c are the same in
                #     basis_b and basis_c
                #
                # (2) signs of m_basis_a are the same in basis_b and basis_c

                if a_dot_b < c_dot_a
                    m_basis_c = basis_b - basis_c
                else
                    m_basis_c = -basis_b + basis_c
                end

            elseif 2 * abs(b_dot_c + b_dot_b) ≈ abs(c_dot_a + a_dot_b)
                # Case:
                #
                # (1) relative signs of m_basis_a and m_basis_c are the same in
                #     basis_b and basis_c
                #
                # (2) signs of m_basis_a are opposite in basis_b and basis_c

                if a_dot_b < -c_dot_a
                    m_basis_c = basis_b + basis_c
                else
                    m_basis_c = -basis_b - basis_c
                end
            end
        end

        if !isnothing(m_basis_a) && !isnothing(m_basis_b) && !isnothing(m_basis_c)
            return m_basis_a, m_basis_b, m_basis_c
        end
    end

    # --- Failed to convert the triclinic unit cell to a base-centered monoclinic unit cell

    throw(
        ErrorException(
            "The triclinic basis vectors defined by `lattice_constants` do not " *
            "satisfy conditions for case 2b.",
        ),
    )
end

"""
    is_triclinic_type_I_cell(lattice_constants::TriclinicLatticeConstants) -> Bool

Determine whether the unit cell defined by `lattice_constants` is a Type I or Type II cell.

Return values
=============
- `true` if `lattice_constants` defines a Type I cell; `false` if `lattice_constants`
  defines a Type II cell.
"""
function is_triclinic_type_I_cell(lattice_constants::TriclinicLatticeConstants)
    α = lattice_constants.α
    β = lattice_constants.β
    γ = lattice_constants.γ

    return cos(α) * cos(β) * cos(γ) >= 0 || abs(cos(α) * cos(β) * cos(γ)) < COS_APPROX_ZERO
end
