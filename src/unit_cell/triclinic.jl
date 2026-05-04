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
Triclinic unit cell types and functions
"""
# --- Imports

# External packages
using AngleBetweenVectors: angle
using Combinatorics: permutations

# --- Exports

# Types
export TriclinicUnitCell, TriclinicUnitCellDelta

# Functions
export satisfies_triclinic_angle_constraints, is_triclinic_type_I_cell
export convert_to_mP, convert_to_mI, convert_to_mC

# Constants
export TRICLINIC_MIN_ANGLE, TRICLINIC_MAX_ANGLE

# --- Constants

const TRICLINIC_MIN_ANGLE = 0
const TRICLINIC_MAX_ANGLE = π

# --- Types

# ------ TriclinicUnitCell

"""
    TriclinicUnitCell

Lattice constants and symmetry for a triclinic unit cell

Fields
======
* `a`, `b`, `c`: lengths of the edges of the unit cell

* `α`, `β`, `γ`: angles between edges of the unit cell
  unit cell

* `symmetry`: unit cell symmetry
"""
const TriclinicUnitCell = UnitCell{Triclinic}

# Outer constructor
"""
    TriclinicUnitCell(
        a::Real, b::Real, c::Real, α::Real, β::Real, γ::Real;
        centering::Centering=primitive_centering,
        symmetry_elements::Union{Set,Vector,Nothing}=nothing,
        check_angle_constraints::Bool=true
    )

Construct a TriclinicUnitCell object from a set of lattice constants.

!!! note

    No constraints are imposed on `centering`. The unit cell does _not_ to be a valid
    Bravais lattice.

Keyword Arguments
=================
- `centering`: centering of unit cell

- `symmetry_elements`: symmetry elements of crystal

- `check_angle_constraints`: if `true`, throw an error if `α`, `β`, and `γ` do not satisfy
  the angle constraints for a valid triclinic unit cell. Otherwise, ignore angle
  constraints.
"""
function TriclinicUnitCell(
    a::Real,
    b::Real,
    c::Real,
    α::Real,
    β::Real,
    γ::Real;
    centering::Centering=primitive_centering,
    symmetry_elements::Union{Set,Vector,Nothing}=nothing,
    check_angle_constraints::Bool=true,
)

    # --- Check arguments

    # lattice constants
    if a <= 0
        throw(DomainError(a, "`a` must be positive"))
    end

    if b <= 0
        throw(DomainError(b, "`b` must be positive"))
    end

    if c <= 0
        throw(DomainError(c, "`c` must be positive"))
    end

    if α <= TRICLINIC_MIN_ANGLE || α >= TRICLINIC_MAX_ANGLE
        throw(DomainError(α, "`α` must satisfy 0 < α < π"))
    end

    if β <= TRICLINIC_MIN_ANGLE || β >= TRICLINIC_MAX_ANGLE
        throw(DomainError(β, "`β` must satisfy 0 < β < π"))
    end

    if γ <= TRICLINIC_MIN_ANGLE || γ >= TRICLINIC_MAX_ANGLE
        throw(DomainError(γ, "`γ` must satisfy 0 < γ < π"))
    end

    # symmetry elements
    # TODO

    # --- Construct and return TriclinicUnitCell object

    return TriclinicUnitCell(
        (a=a, b=b, c=c, α=α, β=β, γ=γ);
        centering=centering,
        symmetry_elements=symmetry_elements,
        check_angle_constraints=check_angle_constraints,
    )
end

"""
    UnitCell{Triclinic}(
        lattice_constants::NamedTuple;
        centering::Centering=primitive_centering,
        symmetry_elements::Union{Set,Vector,Nothing}=nothing,
        check_angle_constraints::Bool=true
    ) -> UnitCell{Triclinic}

Construct a UnitCell from a set of lattice constants. The lattice system `T` must be
consistent with the fields present in `lattice_constants`.

Keyword Arguments
=================
- `centering`: centering of unit cell

- `symmetry_elements`: symmetry elements of crystal

- `check_angle_constraints`: if `true`, throw an error if `α`, `β`, and `γ` do not satisfy
  the angle constraints for a valid triclinic unit cell. Otherwise, ignore angle
  constraints.
"""
@inline function UnitCell{Triclinic}(
    lattice_constants::NamedTuple;
    centering::Centering=primitive_centering,
    symmetry_elements::Union{Set,Vector,Nothing}=nothing,
    check_angle_constraints::Bool=true,
)
    # --- Check arguments

    if (
        check_angle_constraints &&
        !satisfies_triclinic_angle_constraints(
            lattice_constants.α, lattice_constants.β, lattice_constants.γ
        )
    )
        throw(
            DomainError(
                (lattice_constants.α, lattice_constants.β, lattice_constants.γ),
                "`α`, `β`, and `γ` do not satisfy the angle constraints for " *
                "a triclnic unit cell",
            ),
        )
    end

    # --- Construct and return TriclinicUnitCell object

    return UnitCell{Triclinic}(
        lattice_constants, UnitCellSymmetry(centering; symmetry_elements=symmetry_elements)
    )
end

# ------ TriclinicUnitCellDelta

"""
    TriclinicUnitCellDelta

Lattice constant deltas for a triclinic unit cell

Fields
======
* `Δa`, `Δb`, `Δc`: deltas of the lengths of the edges of the unit cell

* `Δα`, `Δβ`, `Δγ`: deltas of the angles between edges of the unit cell
"""
const TriclinicUnitCellDelta = UnitCellDelta{Triclinic}

# Outer constructors
"""
    TriclinicUnitCellDelta(Δa::Real, Δb::Real, Δc::Real, Δα::Real, Δβ::Real, Δγ::Real)

Construct a TriclinicUnitCellDelta object from a set of lattice constant deltas.
"""
function TriclinicUnitCellDelta(Δa::Real, Δb::Real, Δc::Real, Δα::Real, Δβ::Real, Δγ::Real)
    Δlattice_constants = (Δa=Δa, Δb=Δb, Δc=Δc, Δα=Δα, Δβ=Δβ, Δγ=Δγ)
    return TriclinicUnitCellDelta(Δlattice_constants)
end

# --- Functions/Methods

# ------ UnitCell methods

using LinearAlgebra: dot

function standardize(unit_cell::TriclinicUnitCell)
    # --- Check arguments

    standardize_check_args(unit_cell)

    # --- Preparations

    # Extract lattice constants
    lattice_constants_ = lattice_constants(unit_cell)
    a = lattice_constants_.a
    b = lattice_constants_.b
    c = lattice_constants_.c
    α = lattice_constants_.α
    β = lattice_constants_.β
    γ = lattice_constants_.γ

    # --- Standardize lattice constants

    # Shift origin of basis to "homogeneous corner"
    if is_triclinic_type_I_cell(unit_cell)
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

    return TriclinicUnitCell(a, b, c, α, β, γ; centering=primitive_centering)
end

function basis(unit_cell::TriclinicUnitCell)
    # Get lattice constants
    lattice_constants_ = lattice_constants(unit_cell)
    a = lattice_constants_.a
    b = lattice_constants_.b
    c = lattice_constants_.c
    α = lattice_constants_.α
    β = lattice_constants_.β
    γ = lattice_constants_.γ

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

function volume(unit_cell::TriclinicUnitCell)
    # Get lattice constants
    lattice_constants_ = lattice_constants(unit_cell)
    a = lattice_constants_.a
    b = lattice_constants_.b
    c = lattice_constants_.c
    α = lattice_constants_.α
    β = lattice_constants_.β
    γ = lattice_constants_.γ

    # Compute volume
    return (2 * a * b * c) * sqrt(
        sin(0.5 * (α + β + γ)) *
        sin(0.5 * (α + β - γ)) *
        sin(0.5 * (α - β + γ)) *
        sin(0.5 * (-α + β + γ)),
    )
end

function surface_area(unit_cell::TriclinicUnitCell)
    # Get lattice constants
    lattice_constants_ = lattice_constants(unit_cell)
    a = lattice_constants_.a
    b = lattice_constants_.b
    c = lattice_constants_.c
    α = lattice_constants_.α
    β = lattice_constants_.β
    γ = lattice_constants_.γ

    # Compute surface area
    return 2 * (a * b * sin(γ) + b * c * sin(α) + c * a * sin(β))
end

function conventional_cell(::Triclinic, unit_cell::UnitCell)
    # --- Check arguments

    conventional_cell_check_args(unit_cell)

    # --- Preparations

    # Standardize unit cell
    unit_cell = standardize(unit_cell)

    # --- Compute IUCr conventional cell

    # Check limiting case: monoclinic, primitive-centering
    try
        monoclinic_unit_cell = convert_to_mP(unit_cell)

        @debug "aP --> mP"
        return conventional_cell(monoclinic_unit_cell)

    catch error
        if !(error isa ErrorException) || (
            error.msg !=
            "The triclinic unit cell defined by `unit_cell` is not equivalent " *
            "to a primitive monoclinic unit cell."
        )
            rethrow(error)
        end
    end

    # Check limiting case: monoclinic, body-centering
    try
        monoclinic_unit_cell = convert_to_mI(unit_cell)

        @debug "aP --> mI"
        return conventional_cell(monoclinic_unit_cell)

    catch error
        if !(error isa ErrorException) || (
            error.msg !=
            "The triclinic unit cell defined by `unit_cell` is not equivalent " *
            "to a body-centered monoclinic unit cell."
        )
            rethrow(error)
        end
    end

    # Check limiting case: monoclinic, base-centering
    try
        monoclinic_unit_cell = convert_to_mC(unit_cell)

        @debug "aP --> mC"
        body_centered_unit_cell = standardize(monoclinic_unit_cell)
        return conventional_cell(body_centered_unit_cell)

    catch error
        if !(error isa ErrorException) || (
            error.msg !=
            "The triclinic unit cell defined by `unit_cell` is not equivalent " *
            "to a base-centered monoclinic unit cell."
        )
            rethrow(error)
        end
    end

    # Not a limiting case, so return unit cell with standardized lattice constants
    return UnitCell(unit_cell)
end

"""
    convert_to_mP(unit_cell::TriclinicUnitCell) -> MonoclinicUnitCell

Attempt to convert the triclinic unit cell defined by `unit_cell` to an equivalent
primitive monoclinic unit cell.

Return values
=============
- unit cell for the equivalent primitive monoclinic unit cell if one exists

Exceptions
==========
Throws an `ErrorException` if the triclinic unit cell defined by `unit_cell` is not
equivalent to a primitive monoclinic unit cell.
"""
function convert_to_mP(unit_cell::TriclinicUnitCell)
    # --- Preparations

    # Get lattice constants
    lattice_constants_ = lattice_constants(unit_cell)
    a = lattice_constants_.a
    b = lattice_constants_.b
    c = lattice_constants_.c
    α = lattice_constants_.α
    β = lattice_constants_.β
    γ = lattice_constants_.γ

    # --- Attempt to convert the triclinic unit cell to a primitive monoclinic unit cell

    if β ≈ π / 2 && γ ≈ π / 2
        # `a` is the unique symmetry direction of the monoclinic unit cell
        return standardize(MonoclinicUnitCell(b, a, c, α; centering=primitive_centering))

    elseif α ≈ π / 2 && γ ≈ π / 2
        # `b` is the unique symmetry direction of the monoclinic unit cell
        return standardize(MonoclinicUnitCell(a, b, c, β; centering=primitive_centering))

    elseif α ≈ π / 2 && β ≈ π / 2
        # `c` is the unique symmetry direction of the monoclinic unit cell
        return standardize(MonoclinicUnitCell(a, c, b, γ; centering=primitive_centering))
    end

    throw(
        ErrorException(
            "The triclinic unit cell defined by `unit_cell` is not equivalent " *
            "to a primitive monoclinic unit cell.",
        ),
    )
end

"""
    convert_to_mI(unit_cell::TriclinicUnitCell) -> MonoclinicUnitCell

Attempt to convert the triclinic unit cell defined by `unit_cell` to an equivalent
body-centered monoclinic unit cell.

Return values
=============
- unit cell for the equivalent body-centered monoclinic unit cell if one exists

Exceptions
==========
Throws an `ErrorException` if the triclinic unit cell defined by `unit_cell` is not
equivalent to a body-centered monoclinic unit cell.
"""
function convert_to_mI(unit_cell::TriclinicUnitCell)
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
        m_basis_a, m_basis_b, m_basis_c = convert_to_mI_case_1(unit_cell)

        @debug "aP --> mI (case 1)"
        return convert_to_mI_basis_to_unit_cell(m_basis_a, m_basis_b, m_basis_c)

    catch error
        if !(error isa ErrorException) || (
            error.msg !=
            "The triclinic basis vectors defined by `unit_cell` do not " *
            "include the unique monoclinic symmetry direction, one other monoclinic " *
            "basis vector, and one body-centered lattice vector."
        )
            rethrow(error)
        end
    end

    # Case: triclinic basis does not contain unique monoclinic symmetry direction
    try
        # Compute monoclinic basis vectors
        m_basis_a, m_basis_b, m_basis_c = convert_to_mI_case_2(unit_cell)

        @debug "aP --> mI (case 2)"
        return convert_to_mI_basis_to_unit_cell(m_basis_a, m_basis_b, m_basis_c)

    catch error
        if !(error isa ErrorException) || (
            error.msg !=
            "The triclinic basis vectors defined by `unit_cell` do not " *
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
        m_basis_a, m_basis_b, m_basis_c = convert_to_mI_case_3(unit_cell)

        @debug "aP --> mI (case 3)"
        return convert_to_mI_basis_to_unit_cell(m_basis_a, m_basis_b, m_basis_c)

    catch error
        if !(error isa ErrorException) || (
            error.msg !=
            "The triclinic basis vectors defined by `unit_cell` do not " *
            "include the unique monoclinic symmetry direction and two body-centered " *
            "lattice vectors."
        )
            rethrow(error)
        end
    end

    # Case: triclinic basis contains the unique monoclinic symmetry direction
    try
        # Compute monoclinic basis vectors
        m_basis_a, m_basis_b, m_basis_c = convert_to_mI_case_4(unit_cell)

        @debug "aP --> mI (case 4)"
        return convert_to_mI_basis_to_unit_cell(m_basis_a, m_basis_b, m_basis_c)

    catch error
        if !(error isa ErrorException) || (
            error.msg !=
            "The triclinic basis vectors defined by `unit_cell` do not " *
            "include one monoclinic basis vector in the B-face and two body-centered " *
            "lattice vectors."
        )
            rethrow(error)
        end
    end

    # ------ Case: triclinic unit cell basis contains 3 body-centered lattice vectors

    try
        # Compute monoclinic basis vectors
        m_basis_a, m_basis_b, m_basis_c = convert_to_mI_case_5(unit_cell)

        @debug "aP --> mI (case 5)"
        return convert_to_mI_basis_to_unit_cell(m_basis_a, m_basis_b, m_basis_c)

    catch error
        if !(error isa ErrorException) || (
            error.msg !=
            "The triclinic basis vectors defined by `unit_cell` do not " *
            "include three body-centered lattice vectors."
        )
            rethrow(error)
        end
    end

    # --- Throw exception if unit cell is not equivalent to a body-centered monoclinic unit

    if isnothing(m_basis_a) || isnothing(m_basis_b) || isnothing(m_basis_c)
        throw(
            ErrorException(
                "The triclinic unit cell defined by `unit_cell` is not " *
                "equivalent to a body-centered monoclinic unit cell.",
            ),
        )
    end
end

function convert_to_mI_basis_to_unit_cell(
    m_basis_a::Vector{<:Real}, m_basis_b::Vector{<:Real}, m_basis_c::Vector{<:Real}
)
    # Compute monoclinic lattice constants
    m_a = norm(m_basis_a)
    m_b = norm(m_basis_b)
    m_c = norm(m_basis_c)
    m_β = angle(m_basis_a, m_basis_c)

    # Check IUCr conventions
    if m_β < π / 2
        throw(ErrorException("m_β be at least π/2. (m_β=$m_β)"))
    end

    return MonoclinicUnitCell(m_a, m_b, m_c, m_β; centering=body_centering)
end

function convert_to_mI_case_1(unit_cell::TriclinicUnitCell)
    # Notes
    # =====
    # - Case
    #   - triclinic unit cell basis contains 2 monoclinic unit cell basis vectors and
    #     1 body-centered lattice vector
    #
    #   - triclinic basis contains the unique monoclinic symmetry direction m_basis_b
    #
    #   - basis_a and basis_b are m_basis_a and m_basis_b, respectively
    #
    # - Conditions
    #
    #   - abs(a_dot_b) = 0
    #
    #   - 2 * abs(b_dot_c) = b_dot_b
    #
    # - This method adopts the same variable conventions as convert_to_mI().

    for (basis_a, basis_b, basis_c) in permutations(basis(unit_cell))
        # Initialize monoclinic basis vectors
        m_basis_a = nothing
        m_basis_b = nothing
        m_basis_c = nothing

        a_dot_a = dot(basis_a, basis_a)
        b_dot_b = dot(basis_b, basis_b)

        a_dot_b = dot(basis_a, basis_b)
        b_dot_c = dot(basis_b, basis_c)
        c_dot_a = dot(basis_c, basis_a)

        if abs(a_dot_b) < sqrt(a_dot_a) * sqrt(b_dot_b) * COS_APPROX_ZERO &&
            2 * abs(b_dot_c) ≈ b_dot_b

            # Compute monoclinic basis
            m_basis_a = basis_a
            m_basis_b = basis_b

            if 2 * abs(c_dot_a) > a_dot_a
                if c_dot_a > 0
                    m_basis_c = m_basis_a - 2 * basis_c
                else
                    m_basis_c = m_basis_a + 2 * basis_c
                end
            else
                m_basis_c = -m_basis_a + 2 * basis_c
            end
            m_basis_c -= dot(m_basis_c, m_basis_b) / b_dot_b * m_basis_b
        end

        if !isnothing(m_basis_a) && !isnothing(m_basis_b) && !isnothing(m_basis_c)
            return m_basis_a, m_basis_b, m_basis_c
        end
    end

    # --- Failed to convert the triclinic unit cell to a base-centered monoclinic unit cell

    throw(
        ErrorException(
            "The triclinic basis vectors defined by `unit_cell` do not include " *
            "the unique monoclinic symmetry direction, one other monoclinic basis " *
            "vector, and one body-centered lattice vector.",
        ),
    )
end

function convert_to_mI_case_2(unit_cell::TriclinicUnitCell)
    # --- Preparations

    # Initialize monoclinic basis vectors
    m_basis_a = nothing
    m_basis_b = nothing
    m_basis_c = nothing

    # --- Attempt to convert the triclinic unit cell to a body-centered monoclinic unit cell

    try
        # Compute monoclinic basis vectors
        m_basis_a, m_basis_b, m_basis_c = convert_to_mI_case_2a(unit_cell)
        @debug "aP --> mI (case 2a)"

    catch error
        if !(error isa ErrorException) || (
            error.msg !=
            "The triclinic basis vectors defined by `unit_cell` do not " *
            "satisfy conditions for case 2a."
        )
            rethrow(error)
        end
    end

    try
        # Compute monoclinic basis vectors
        m_basis_a, m_basis_b, m_basis_c = convert_to_mI_case_2b(unit_cell)
        @debug "aP --> mI (case 2b)"

    catch error
        if !(error isa ErrorException) || (
            error.msg !=
            "The triclinic basis vectors defined by `unit_cell` do not " *
            "satisfy conditions for case 2b."
        )
            rethrow(error)
        end
    end

    # --- Construct return value

    if isnothing(m_basis_a) || isnothing(m_basis_b) || isnothing(m_basis_c)
        throw(
            ErrorException(
                "The triclinic basis vectors defined by `unit_cell` do not " *
                "include two monoclinic basis vectors in the B-face and one " *
                "body-centered lattice vector.",
            ),
        )
    end

    return m_basis_a, m_basis_b, m_basis_c
end

function convert_to_mI_case_2a(unit_cell::TriclinicUnitCell)
    # Notes
    # =====
    # - Case
    #   - triclinic unit cell basis contains 2 monoclinic unit cell basis vectors and
    #     1 body-centered lattice vector
    #
    #   - triclinic basis contains m_basis_a and m_basis_c (not the unique monoclinic
    #     symmetry direction m_basis_b)
    #
    #   - basis_b is the vector to body-centered lattice point in monoclinic unit cell
    #
    #   - basis_a and basis_c are m_basis_a and m_basis_c, respectively
    #
    #   - m_basis_a and m_basis_c have the same sign in basis_b
    #
    # - Conditions
    #
    #   - 2 * abs(a_dot_b) = abs(a_dot_a + c_dot_a)
    #
    #   - 2 * abs(b_dot_c) = abs(c_dot_c + c_dot_a)
    #
    #   - abs(c_dot_a) != 0
    #
    # - This method adopts the same variable conventions as convert_to_mI().

    for (basis_a, basis_b, basis_c) in permutations(basis(unit_cell))
        # Initialize monoclinic basis vectors
        m_basis_a = nothing
        m_basis_b = nothing
        m_basis_c = nothing

        a_dot_a = dot(basis_a, basis_a)
        c_dot_c = dot(basis_c, basis_c)

        a_dot_b = dot(basis_a, basis_b)
        b_dot_c = dot(basis_b, basis_c)
        c_dot_a = dot(basis_c, basis_a)

        if (
            2 * abs(a_dot_b) ≈ abs(a_dot_a + c_dot_a) &&
            2 * abs(b_dot_c) ≈ abs(c_dot_c + c_dot_a) &&
            abs(c_dot_a) > sqrt(c_dot_c) * sqrt(a_dot_a) * COS_APPROX_ZERO
        )
            # Compute monoclinic basis vectors
            m_basis_a = basis_a

            if 2 * a_dot_b ≈ a_dot_a + c_dot_a
                m_basis_b = 2 * basis_b - basis_a - basis_c
            else
                m_basis_b = 2 * basis_b + basis_a + basis_c
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
            "The triclinic basis vectors defined by `unit_cell` do not satisfy " *
            "conditions for case 2a.",
        ),
    )
end

function convert_to_mI_case_2b(unit_cell::TriclinicUnitCell)
    # Notes
    # =====
    # - Case
    #   - triclinic unit cell basis contains 2 monoclinic unit cell basis vectors and
    #     1 body-centered lattice vector
    #
    #   - triclinic basis contains m_basis_a and m_basis_c (not the unique monoclinic
    #     symmetry direction m_basis_b)
    #
    #   - basis_b is the vector to body-centered lattice point in monoclinic unit cell
    #
    #   - basis_a and basis_c are m_basis_a and m_basis_c, respectively
    #
    #   - m_basis_a and m_basis_c have the opposite signs in basis_b
    #
    # - Conditions
    #
    #   - 2 * abs(a_dot_b) = abs(a_dot_a - c_dot_a)
    #
    #   - 2 * abs(b_dot_c) = abs(c_dot_c - c_dot_a)
    #
    #   - abs(c_dot_a) != 0
    #
    # - This method adopts the same variable conventions as convert_to_mI().

    for (basis_a, basis_b, basis_c) in permutations(basis(unit_cell))
        # Initialize monoclinic basis vectors
        m_basis_a = nothing
        m_basis_b = nothing
        m_basis_c = nothing

        a_dot_a = dot(basis_a, basis_a)
        c_dot_c = dot(basis_c, basis_c)

        a_dot_b = dot(basis_a, basis_b)
        b_dot_c = dot(basis_b, basis_c)
        c_dot_a = dot(basis_c, basis_a)

        if (
            2 * abs(a_dot_b) ≈ abs(a_dot_a - c_dot_a) &&
            2 * abs(b_dot_c) ≈ abs(c_dot_c - c_dot_a) &&
            abs(c_dot_a) > sqrt(c_dot_c) * sqrt(a_dot_a) * COS_APPROX_ZERO
        )
            # Compute monoclinic basis vectors
            m_basis_a = basis_a

            if 2 * a_dot_b ≈ a_dot_a - c_dot_a
                m_basis_b = 2 * basis_b - basis_a + basis_c
            else
                m_basis_b = 2 * basis_b + basis_a - basis_c
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
            "The triclinic basis vectors defined by `unit_cell` do not satisfy " *
            "conditions for case 2b.",
        ),
    )
end

function convert_to_mI_case_3(unit_cell::TriclinicUnitCell)
    # Notes
    # =====
    # - Case
    #   - triclinic unit cell basis contains 1 monoclinic unit cell basis vector and
    #     2 body-centered lattice vectors
    #
    #   - triclinic basis contains the unique monoclinic symmetry direction m_basis_b
    #
    # - Conditions
    #
    #   - abs(b_dot_c) = abs(a_dot_b) = 0.5 * b_dot_b
    #
    # - This method adopts the same variable conventions as convert_to_mI().

    for (basis_a, basis_b, basis_c) in permutations(basis(unit_cell))
        # Initialize monoclinic basis vectors
        m_basis_a = nothing
        m_basis_b = nothing
        m_basis_c = nothing

        a_dot_a = dot(basis_a, basis_a)
        b_dot_b = dot(basis_b, basis_b)
        c_dot_c = dot(basis_c, basis_c)

        a_dot_b = dot(basis_a, basis_b)
        b_dot_c = dot(basis_b, basis_c)

        if abs(b_dot_c) ≈ abs(a_dot_b) ≈ 0.5 * b_dot_b
            # Compute monoclinic basis vectors
            m_basis_b = basis_b

            m_basis_a = basis_a + basis_c
            m_basis_a -= dot(m_basis_a, m_basis_b) / b_dot_b * m_basis_b

            if a_dot_a < c_dot_c
                m_basis_c = basis_a - basis_c
            else
                m_basis_c = -basis_a + basis_c
            end
            m_basis_c -= dot(m_basis_c, m_basis_b) / b_dot_b * m_basis_b
        end

        if !isnothing(m_basis_a) && !isnothing(m_basis_b) && !isnothing(m_basis_c)
            return m_basis_a, m_basis_b, m_basis_c
        end
    end

    # --- Failed to convert the triclinic unit cell to a base-centered monoclinic unit cell

    throw(
        ErrorException(
            "The triclinic basis vectors defined by `unit_cell` do not " *
            "include the unique monoclinic symmetry direction and two body-centered " *
            "lattice vectors.",
        ),
    )
end

function convert_to_mI_case_4(unit_cell::TriclinicUnitCell)
    # --- Preparations

    # Initialize monoclinic basis vectors
    m_basis_a = nothing
    m_basis_b = nothing
    m_basis_c = nothing

    # --- Attempt to convert the triclinic unit cell to a body-centered monoclinic unit cell

    try
        # Compute monoclinic basis vectors
        m_basis_a, m_basis_b, m_basis_c = convert_to_mI_case_4a(unit_cell)
        @debug "aP --> mI (case 4a)"

    catch error
        if !(error isa ErrorException) || (
            error.msg !=
            "The triclinic basis vectors defined by `unit_cell` do not " *
            "satisfy conditions for case 4a."
        )
            rethrow(error)
        end
    end

    try
        # Compute monoclinic basis vectors
        m_basis_a, m_basis_b, m_basis_c = convert_to_mI_case_4b(unit_cell)
        @debug "aP --> mI (case 4b)"

    catch error
        if !(error isa ErrorException) || (
            error.msg !=
            "The triclinic basis vectors defined by `unit_cell` do not " *
            "satisfy conditions for case 4b."
        )
            rethrow(error)
        end
    end

    try
        # Compute monoclinic basis vectors
        m_basis_a, m_basis_b, m_basis_c = convert_to_mI_case_4c(unit_cell)

    catch error
        if !(error isa ErrorException) || (
            error.msg !=
            "The triclinic basis vectors defined by `unit_cell` do not " *
            "satisfy conditions for case 4c."
        )
            rethrow(error)
        end
    end

    # --- Construct return value

    if isnothing(m_basis_a) || isnothing(m_basis_b) || isnothing(m_basis_c)
        throw(
            ErrorException(
                "The triclinic basis vectors defined by `unit_cell` do not " *
                "include one monoclinic basis vector in the B-face and two body-centered " *
                "lattice vectors.",
            ),
        )
    end

    return m_basis_a, m_basis_b, m_basis_c
end

function convert_to_mI_case_4a(unit_cell::TriclinicUnitCell)
    # Notes
    # =====
    # - Case
    #   - triclinic unit cell basis contains 1 monoclinic unit cell basis vector and
    #     2 body-centered lattice vectors
    #
    #   - triclinic basis contains m_basis_a (not the unique monoclinic symmetry direction
    #     m_basis_b)
    #
    #   - the coefficients of m_basis_a and m_basis_c have the same relative signs in
    #     both basis_b and basis_c
    #
    # - Conditions
    #
    #   - b_dot_b = c_dot_c
    #
    #   - abs(a_dot_b) = abs(c_dot_a) != 0.5 * a_dot_a
    #
    # - This method adopts the same variable conventions as convert_to_mI().

    for (basis_a, basis_b, basis_c) in permutations(basis(unit_cell))
        # Initialize monoclinic basis vectors
        m_basis_a = nothing
        m_basis_b = nothing
        m_basis_c = nothing

        a_dot_a = dot(basis_a, basis_a)
        b_dot_b = dot(basis_b, basis_b)
        c_dot_c = dot(basis_c, basis_c)

        a_dot_b = dot(basis_a, basis_b)
        c_dot_a = dot(basis_c, basis_a)

        if (b_dot_b ≈ c_dot_c) && (abs(a_dot_b) ≈ abs(c_dot_a) ≉ 0.5 * a_dot_a)
            # Compute m_basis_a
            m_basis_a = basis_a

            # Compute m_basis_b and m_basis_c
            if a_dot_b ≈ c_dot_a
                # Case: coefficients of m_basis_a and m_basis_c are the same in basis_b
                #       and basis_c, so the sign of m_basis_b must be opposite in basis_b
                #       and basis_c

                # Compute m_basis_b
                m_basis_b = basis_b - basis_c

                # Compute m_basis_c
                if 2 * a_dot_b < a_dot_a
                    m_basis_c = basis_b + basis_c - basis_a
                else
                    m_basis_c = -(basis_b + basis_c - basis_a)
                end
            else
                # Case: coefficients of m_basis_a and m_basis_c are opposite in basis_b
                #       and basis_c, so the sign of m_basis_b must be the same in basis_b
                #       and basis_c

                # Compute m_basis_b
                m_basis_b = basis_b + basis_c

                # Compute m_basis_c
                if 2 * a_dot_b < a_dot_a
                    m_basis_c = basis_b - basis_c - basis_a
                else
                    m_basis_c = -(basis_b - basis_c - basis_a)
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
            "The triclinic basis vectors defined by `unit_cell` do not satisfy " *
            "conditions for case 4a.",
        ),
    )
end

function convert_to_mI_case_4b(unit_cell::TriclinicUnitCell)
    # Notes
    # =====
    # - Case
    #   - triclinic unit cell basis contains 1 monoclinic unit cell basis vector and
    #     2 body-centered lattice vectors
    #
    #   - triclinic basis contains m_basis_a (not the unique monoclinic symmetry direction
    #     m_basis_b)
    #     m_basis_b
    #
    #   - the coefficients of m_basis_a in basis_b and basis_c have the same sign
    #
    #   - the coefficients of m_basis_c in basis_b and basis_c have opposite sign
    #
    # - Conditions
    #
    #   - abs(a_dot_b + c_dot_a) ≈ a_dot_a
    #
    #   - abs(a_dot_b) != abs(c_dot_a)
    #
    # - This method adopts the same variable conventions as convert_to_mI().

    for (basis_a, basis_b, basis_c) in permutations(basis(unit_cell))
        # Initialize monoclinic basis vectors
        m_basis_a = nothing
        m_basis_b = nothing
        m_basis_c = nothing

        a_dot_a = dot(basis_a, basis_a)

        a_dot_b = dot(basis_a, basis_b)
        c_dot_a = dot(basis_c, basis_a)

        if (abs(a_dot_b + c_dot_a) ≈ a_dot_a) && (abs(a_dot_b) ≉ abs(c_dot_a))
            # Compute monoclinic basis vectors
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
        end

        if !isnothing(m_basis_a) && !isnothing(m_basis_b) && !isnothing(m_basis_c)
            return m_basis_a, m_basis_b, m_basis_c
        end
    end

    # --- Failed to convert the triclinic unit cell to a base-centered monoclinic unit cell

    throw(
        ErrorException(
            "The triclinic basis vectors defined by `unit_cell` do not satisfy " *
            "conditions for case 4b.",
        ),
    )
end

function convert_to_mI_case_4c(unit_cell::TriclinicUnitCell)
    # Notes
    # =====
    # - Case
    #   - triclinic unit cell basis contains 1 monoclinic unit cell basis vector and
    #     2 body-centered lattice vectors
    #
    #   - triclinic basis contains m_basis_a (not the unique monoclinic symmetry direction
    #     m_basis_b)
    #
    #   - the coefficients of m_basis_a in basis_b and basis_c have opposite sign
    #
    #   - the coefficients of m_basis_c in basis_b and basis_c have the same sign
    #
    # - Conditions
    #
    #   - abs(a_dot_b - c_dot_a) ≈ a_dot_a
    #
    #   - abs(a_dot_b) != abs(c_dot_a)
    #
    # - This method adopts the same variable conventions as convert_to_mI().

    for (basis_a, basis_b, basis_c) in permutations(basis(unit_cell))
        # Initialize monoclinic basis vectors
        m_basis_a = nothing
        m_basis_b = nothing
        m_basis_c = nothing

        a_dot_a = dot(basis_a, basis_a)

        a_dot_b = dot(basis_a, basis_b)
        c_dot_a = dot(basis_c, basis_a)

        if (abs(a_dot_b - c_dot_a) ≈ a_dot_a) && (abs(a_dot_b) ≉ abs(c_dot_a))

            # Compute monoclinic basis vectors
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
        end

        if !isnothing(m_basis_a) && !isnothing(m_basis_b) && !isnothing(m_basis_c)
            return m_basis_a, m_basis_b, m_basis_c
        end
    end

    # --- Failed to convert the triclinic unit cell to a base-centered monoclinic unit cell

    throw(
        ErrorException(
            "The triclinic basis vectors defined by `unit_cell` do not satisfy " *
            "conditions for case 4c.",
        ),
    )
end

function convert_to_mI_case_5(unit_cell::TriclinicUnitCell)
    # Notes
    # =====
    # - Case
    #   - triclinic unit cell basis contains 3 body-centered lattice vectors
    #
    #   - basis_b and basis_c are the triclinic basis vector with the same length
    #
    # - Conditions
    #
    #   - b_dot_b = c_dot_c
    #
    #   - abs(a_dot_b) != abs(c_dot_a)
    #
    # - This method adopts the same variable conventions as convert_to_mI().

    for (basis_a, basis_b, basis_c) in permutations(basis(unit_cell))
        # Initialize monoclinic basis vectors
        m_basis_a = nothing
        m_basis_b = nothing
        m_basis_c = nothing

        b_dot_b = dot(basis_b, basis_b)
        c_dot_c = dot(basis_c, basis_c)

        a_dot_b = dot(basis_a, basis_b)
        b_dot_c = dot(basis_b, basis_c)
        c_dot_a = dot(basis_c, basis_a)

        if b_dot_b ≈ c_dot_c && abs(a_dot_b) ≉ abs(c_dot_a)

            # Compute m_basis_b and precursors to m_basis_a and m_basis_c
            m_basis_d_1 = nothing
            m_basis_d_2 = nothing

            if 2 * abs(a_dot_b + c_dot_a) ≈ b_dot_b + c_dot_c + 2 * b_dot_c
                m_basis_b = basis_b + basis_c
                m_basis_d_1 = basis_b - basis_c

            elseif 2 * abs(a_dot_b - c_dot_a) ≈ b_dot_b + c_dot_c - 2 * b_dot_c
                m_basis_b = basis_b - basis_c
                m_basis_d_1 = basis_b + basis_c
            end

            m_basis_d_2 = 2 * basis_a
            m_basis_d_2 -=
                dot(m_basis_d_2, m_basis_b) / dot(m_basis_b, m_basis_b) * m_basis_b

            # Compute m_basis_a and m_basis_c
            if !isnothing(m_basis_d_1) && !isnothing(m_basis_d_2)
                m_basis_a = 0.5 * (m_basis_d_1 + m_basis_d_2)

                if dot(m_basis_d_1, m_basis_d_1) < dot(m_basis_d_2, m_basis_d_2)
                    m_basis_c = 0.5 * (m_basis_d_1 - m_basis_d_2)
                else
                    m_basis_c = 0.5 * (-m_basis_d_1 + m_basis_d_2)
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
            "The triclinic basis vectors defined by `unit_cell` do not include " *
            "three body-centered lattice vectors.",
        ),
    )
end

"""
    convert_to_mC(
        unit_cell::TriclinicUnitCell
    ) -> MonoclinicUnitCell

Attempt to convert the triclinic unit cell defined by `unit_cell` to an equivalent
base-centered monoclinic unit cell.

Return values
=============
- lattice constants for the equivalent base-centered monoclinic unit cell if one exists;
  `nothing` otherwise

  !!! warn

      Returned lattice constants are _not_ guaranteed to be standardized.

Exceptions
==========
Throws an `ErrorException` if the triclinic unit cell defined by `unit_cell` is not
equivalent to a base-centered monoclinic unit cell.
"""
function convert_to_mC(unit_cell::TriclinicUnitCell)
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
        m_basis_a, m_basis_b, m_basis_c = convert_to_mC_case_1(unit_cell)
        return convert_to_mC_basis_to_unit_cell(m_basis_a, m_basis_b, m_basis_c)

    catch error
        if !(error isa ErrorException) || (
            error.msg !=
            "The triclinic basis vectors defined by `unit_cell` do not " *
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
        m_basis_a, m_basis_b, m_basis_c = convert_to_mC_case_2(unit_cell)
        return convert_to_mC_basis_to_unit_cell(m_basis_a, m_basis_b, m_basis_c)

    catch error
        if !(error isa ErrorException) || (
            error.msg !=
            "The triclinic basis vectors defined by `unit_cell` do not " *
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
                "The triclinic unit cell defined by `unit_cell` is not " *
                "equivalent to a base-centered monoclinic unit cell.",
            ),
        )
    end
end

function convert_to_mC_basis_to_unit_cell(
    m_basis_a::Vector{<:Real}, m_basis_b::Vector{<:Real}, m_basis_c::Vector{<:Real}
)
    # Compute monoclinic lattice constants
    m_a = norm(m_basis_a)
    m_b = norm(m_basis_b)
    m_c = norm(m_basis_c)
    m_β = angle(m_basis_a, m_basis_c)

    # Check IUCr conventions
    if m_β < π / 2
        throw(ErrorException("m_β be at least π/2. (m_β=$m_β)"))
    end

    return MonoclinicUnitCell(m_a, m_b, m_c, m_β; centering=base_centering)
end

function convert_to_mC_case_1(unit_cell::TriclinicUnitCell)
    # --- Preparations

    # Initialize monoclinic basis vectors
    m_basis_a = nothing
    m_basis_b = nothing
    m_basis_c = nothing

    # --- Attempt to convert the triclinic unit cell to a base-centered monoclinic unit cell

    try
        # Compute monoclinic basis vectors
        m_basis_a, m_basis_b, m_basis_c = convert_to_mC_case_1a(unit_cell)
        @debug "aP --> mC (case 1a)"

    catch error
        if !(error isa ErrorException) || (
            error.msg !=
            "The triclinic basis vectors defined by `unit_cell` do not " *
            "satisfy conditions for case 1a."
        )
            rethrow(error)
        end
    end

    try
        # Compute monoclinic basis vectors
        m_basis_a, m_basis_b, m_basis_c = convert_to_mC_case_1b(unit_cell)
        @debug "aP --> mC (case 1b)"

    catch error
        if !(error isa ErrorException) || (
            error.msg !=
            "The triclinic basis vectors defined by `unit_cell` do not " *
            "satisfy conditions for case 1b."
        )
            rethrow(error)
        end
    end

    # --- Construct return value

    if isnothing(m_basis_a) || isnothing(m_basis_b) || isnothing(m_basis_c)
        throw(
            ErrorException(
                "The triclinic basis vectors defined by `unit_cell` do not " *
                "include the mononclinic basis vectors that are not the unique " *
                "monoclinic symmetry direction and one base-centered lattice vector.",
            ),
        )
    end

    return m_basis_a, m_basis_b, m_basis_c
end

function convert_to_mC_case_1a(unit_cell::TriclinicUnitCell)
    # Notes
    # =====
    # - Case
    #   - triclinic unit cell basis contains m_basis_a, m_basis_c, and 1 base-centered
    #     lattice vector
    #
    #   - base-centered vector does not include a contribution from m_basis_c
    #
    #   - basis_a and basis_c are m_basis_a and m_basis_c, respectively
    #
    # - Conditions
    #
    #   - 2 * abs(a_dot_b) = a_dot_a
    #
    #   - 2 * abs(b_dot_c) = abs(c_dot_a)
    #
    # - This method adopts the same variable conventions as convert_to_mC().

    # --- Attempt to convert the triclinic unit cell to a base-centered monoclinic unit cell

    for (basis_a, basis_b, basis_c) in permutations(basis(unit_cell))
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

            # Compute monoclinic basis
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
            "The triclinic basis vectors defined by `unit_cell` do not " *
            "satisfy conditions for case 1a.",
        ),
    )
end

function convert_to_mC_case_1b(unit_cell::TriclinicUnitCell)
    # Notes
    # =====
    # - Case
    #   - triclinic unit cell basis contains m_basis_a, m_basis_c, and 1 base-centered
    #     lattice vector
    #
    #   - base-centered vector include a contribution from m_basis_c
    #
    #   - basis_a and basis_c are m_basis_a and m_basis_c, respectively
    #
    # - Conditions
    #
    #   - abs(a_dot_b - 0.5 * a_dot_a) = abs(c_dot_a) or
    #     abs(a_dot_b + 0.5 * a_dot_a) = abs(c_dot_a)
    #
    #   - 2 * abs(b_dot_c - c_dot_c) = abs(c_dot_a) or
    #     2 * abs(b_dot_c + c_dot_c) = abs(c_dot_a)
    #
    # - This method adopts the same variable conventions as convert_to_mC().

    # --- Attempt to convert the triclinic unit cell to a base-centered monoclinic unit cell

    for (basis_a, basis_b, basis_c) in permutations(basis(unit_cell))
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
            "The triclinic basis vectors defined by `unit_cell` do not " *
            "satisfy conditions for case 1b.",
        ),
    )
end

function convert_to_mC_case_2(unit_cell::TriclinicUnitCell)
    # --- Preparations

    # Initialize monoclinic basis vectors
    m_basis_a = nothing
    m_basis_b = nothing
    m_basis_c = nothing

    # --- Attempt to convert the triclinic unit cell to a base-centered monoclinic unit cell

    try
        # Compute monoclinic basis vectors
        m_basis_a, m_basis_b, m_basis_c = convert_to_mC_case_2a(unit_cell)
        @debug "aP --> mC (case 2a)"

    catch error
        if !(error isa ErrorException) || (
            error.msg !=
            "The triclinic basis vectors defined by `unit_cell` do not " *
            "satisfy conditions for case 2a."
        )
            rethrow(error)
        end
    end

    try
        # Compute monoclinic basis vectors
        m_basis_a, m_basis_b, m_basis_c = convert_to_mC_case_2b(unit_cell)
        @debug "aP --> mC (case 2b)"

    catch error
        if !(error isa ErrorException) || (
            error.msg !=
            "The triclinic basis vectors defined by `unit_cell` do not " *
            "satisfy conditions for case 2b."
        )
            rethrow(error)
        end
    end

    try
        # Compute monoclinic basis vectors
        m_basis_a, m_basis_b, m_basis_c = convert_to_mC_case_2c(unit_cell)
        @debug "aP --> mC (case 2c)"

    catch error
        if !(error isa ErrorException) || (
            error.msg !=
            "The triclinic basis vectors defined by `unit_cell` do not " *
            "satisfy conditions for case 2c."
        )
            rethrow(error)
        end
    end

    try
        # Compute monoclinic basis vectors
        m_basis_a, m_basis_b, m_basis_c = convert_to_mC_case_2d(unit_cell)
        @debug "aP --> mC (case 2d)"

    catch error
        if !(error isa ErrorException) || (
            error.msg !=
            "The triclinic basis vectors defined by `unit_cell` do not " *
            "satisfy conditions for case 2d."
        )
            rethrow(error)
        end
    end

    # --- Construct return value

    if isnothing(m_basis_a) || isnothing(m_basis_b) || isnothing(m_basis_c)
        throw(
            ErrorException(
                "The triclinic basis vectors defined by `unit_cell` do not " *
                "include the mononclinic basis vector that is the intersection of the " *
                "B-face and C-face and two base-centered lattice vectors.",
            ),
        )
    end

    return m_basis_a, m_basis_b, m_basis_c
end

function convert_to_mC_case_2a(unit_cell::TriclinicUnitCell)
    # Notes
    # =====
    # - Case
    #   - triclinic unit cell basis contains m_basis_a and 2 base-centered lattice vectors
    #
    #   - the sign of the coefficient of m_basis_a in basis_b is positive and the signs of
    #     the coefficients of m_basis_b in basis_b and basis_c are the same
    #
    # - Conditions
    #
    #   - 2 * a_dot_b = a_dot_a
    #
    #   - 2 * abs(b_dot_c - b_dot_b - 0.5 * a_dot_a + a_dot_b) = abs(c_dot_a - a_dot_b)
    #
    # - This method adopts the same variable conventions as convert_to_mC().

    # --- Attempt to convert the triclinic unit cell to a base-centered monoclinic unit cell

    for (basis_a, basis_b, basis_c) in permutations(basis(unit_cell))
        # Initialize monoclinic basis vectors
        m_basis_a = nothing
        m_basis_b = nothing
        m_basis_c = nothing

        a_dot_a = dot(basis_a, basis_a)
        b_dot_b = dot(basis_b, basis_b)

        a_dot_b = dot(basis_a, basis_b)
        b_dot_c = dot(basis_b, basis_c)
        c_dot_a = dot(basis_c, basis_a)

        if 2 * a_dot_b ≈ a_dot_a &&
            2 * abs(b_dot_c - b_dot_b - 0.5 * a_dot_a + a_dot_b) ≈ abs(c_dot_a - a_dot_b)

            # Compute monoclinic basis
            m_basis_a = basis_a
            m_basis_b = 2 * basis_b - basis_a

            if a_dot_b < c_dot_a
                m_basis_c = basis_b - basis_c
            else
                m_basis_c = -basis_b + basis_c
            end
        end

        if !isnothing(m_basis_a) && !isnothing(m_basis_b) && !isnothing(m_basis_c)
            return m_basis_a, m_basis_b, m_basis_c
        end
    end

    # --- Failed to convert the triclinic unit cell to a base-centered monoclinic unit cell

    throw(
        ErrorException(
            "The triclinic basis vectors defined by `unit_cell` do not " *
            "satisfy conditions for case 2a.",
        ),
    )
end

function convert_to_mC_case_2b(unit_cell::TriclinicUnitCell)
    # Notes
    # =====
    # - Case
    #   - triclinic unit cell basis contains m_basis_a and 2 base-centered lattice vectors
    #
    #   - the sign of the coefficient of m_basis_a in basis_b is positive and the signs of
    #     the coefficients of m_basis_b in basis_b and basis_c are opposite
    #
    # - Conditions
    #
    #   - 2 * a_dot_b = a_dot_a
    #
    #   - 2 * abs(b_dot_c + b_dot_b - a_dot_b) = abs(c_dot_a - a_dot_b)
    #
    # - This method adopts the same variable conventions as convert_to_mC().

    # --- Attempt to convert the triclinic unit cell to a base-centered monoclinic unit cell

    for (basis_a, basis_b, basis_c) in permutations(basis(unit_cell))
        # Initialize monoclinic basis vectors
        m_basis_a = nothing
        m_basis_b = nothing
        m_basis_c = nothing

        a_dot_a = dot(basis_a, basis_a)
        b_dot_b = dot(basis_b, basis_b)

        a_dot_b = dot(basis_a, basis_b)
        b_dot_c = dot(basis_b, basis_c)
        c_dot_a = dot(basis_c, basis_a)

        if 2 * a_dot_b ≈ a_dot_a &&
            2 * abs(b_dot_c + b_dot_b - a_dot_b) ≈ abs(c_dot_a - a_dot_b)

            # Compute monoclinic basis
            m_basis_a = basis_a
            m_basis_b = 2 * basis_b - basis_a

            if a_dot_b < -c_dot_a
                m_basis_c = basis_b + basis_c
            else
                m_basis_c = -basis_b - basis_c
            end
        end

        if !isnothing(m_basis_a) && !isnothing(m_basis_b) && !isnothing(m_basis_c)
            return m_basis_a, m_basis_b, m_basis_c
        end
    end

    # --- Failed to convert the triclinic unit cell to a base-centered monoclinic unit cell

    throw(
        ErrorException(
            "The triclinic basis vectors defined by `unit_cell` do not " *
            "satisfy conditions for case 2b.",
        ),
    )
end

function convert_to_mC_case_2c(unit_cell::TriclinicUnitCell)
    # Notes
    # =====
    # - Case
    #   - triclinic unit cell basis contains m_basis_a and 2 base-centered lattice vectors
    #
    #   - the sign of the coefficient of m_basis_a in basis_b is negative and the signs of
    #     the coefficients of m_basis_b in basis_b and basis_c are the same
    #
    # - Conditions
    #
    #   - 2 * a_dot_b = -a_dot_a
    #
    #   - 2 * abs(b_dot_c - b_dot_b - 0.5 * a_dot_a - a_dot_b) = abs(c_dot_a - a_dot_b)
    #
    # - This method adopts the same variable conventions as convert_to_mC().

    # --- Attempt to convert the triclinic unit cell to a base-centered monoclinic unit cell

    for (basis_a, basis_b, basis_c) in permutations(basis(unit_cell))
        # Initialize monoclinic basis vectors
        m_basis_a = nothing
        m_basis_b = nothing
        m_basis_c = nothing

        a_dot_a = dot(basis_a, basis_a)
        b_dot_b = dot(basis_b, basis_b)

        a_dot_b = dot(basis_a, basis_b)
        b_dot_c = dot(basis_b, basis_c)
        c_dot_a = dot(basis_c, basis_a)

        if 2 * a_dot_b ≈ -a_dot_a &&
            2 * abs(b_dot_c - b_dot_b - 0.5 * a_dot_a - a_dot_b) ≈ abs(c_dot_a - a_dot_b)

            # Compute monoclinic basis
            m_basis_a = basis_a
            m_basis_b = 2 * basis_b + basis_a

            if a_dot_b < c_dot_a
                m_basis_c = basis_b - basis_c
            else
                m_basis_c = -basis_b + basis_c
            end
        end

        if !isnothing(m_basis_a) && !isnothing(m_basis_b) && !isnothing(m_basis_c)
            return m_basis_a, m_basis_b, m_basis_c
        end
    end

    # --- Failed to convert the triclinic unit cell to a base-centered monoclinic unit cell

    throw(
        ErrorException(
            "The triclinic basis vectors defined by `unit_cell` do not " *
            "satisfy conditions for case 2c.",
        ),
    )
end

function convert_to_mC_case_2d(unit_cell::TriclinicUnitCell)
    # Notes
    # =====
    # - Case
    #   - triclinic unit cell basis contains m_basis_a and 2 base-centered lattice vectors
    #
    #   - the sign of the coefficient of m_basis_a in basis_b is negative and the signs of
    #     the coefficients of m_basis_b in basis_b and basis_c are opposite
    #
    # - Conditions
    #
    #   - 2 * a_dot_b ≈ -a_dot_a
    #
    #   - 2 * abs(b_dot_c + b_dot_b + a_dot_b) = abs(c_dot_a - a_dot_b)
    #
    # - This method adopts the same variable conventions as convert_to_mC().

    # --- Attempt to convert the triclinic unit cell to a base-centered monoclinic unit cell

    for (basis_a, basis_b, basis_c) in permutations(basis(unit_cell))
        # Initialize monoclinic basis vectors
        m_basis_a = nothing
        m_basis_b = nothing
        m_basis_c = nothing

        a_dot_a = dot(basis_a, basis_a)
        b_dot_b = dot(basis_b, basis_b)

        a_dot_b = dot(basis_a, basis_b)
        b_dot_c = dot(basis_b, basis_c)
        c_dot_a = dot(basis_c, basis_a)

        if 2 * a_dot_b ≈ -a_dot_a &&
            2 * abs(b_dot_c + b_dot_b + a_dot_b) ≈ abs(c_dot_a - a_dot_b)

            # Compute monoclinic basis
            m_basis_a = basis_a
            m_basis_b = 2 * basis_b + basis_a

            if a_dot_b < -c_dot_a
                m_basis_c = basis_b + basis_c
            else
                m_basis_c = -basis_b - basis_c
            end
        end

        if !isnothing(m_basis_a) && !isnothing(m_basis_b) && !isnothing(m_basis_c)
            return m_basis_a, m_basis_b, m_basis_c
        end
    end

    # --- Failed to convert the triclinic unit cell to a base-centered monoclinic unit cell

    throw(
        ErrorException(
            "The triclinic basis vectors defined by `unit_cell` do not " *
            "satisfy conditions for case 2d.",
        ),
    )
end

"""
    is_triclinic_type_I_cell(unit_cell::TriclinicUnitCell) -> Bool

Determine whether the unit cell defined by `unit_cell` is a Type I or Type II cell.

A triclinic unit cell is Type I if the product of the dot products of all pairs of basis
vectors for unit cell is positive:

```math
(\\vec{a} \\cdot \\vec{b})(\\vec{b} \\cdot \\vec{c})(\\vec{c} \\cdot \\vec{a}) > 0.
```

Otherwise, the triclinic unit cell is Type II.

Return values
=============
- `true` if `unit_cell` defines a Type I cell; `false` if `unit_cell` defines a Type II
  cell
"""
function is_triclinic_type_I_cell(unit_cell::TriclinicUnitCell)
    # Extract lattice angles
    lattice_constants_ = lattice_constants(unit_cell)
    α = lattice_constants_.α
    β = lattice_constants_.β
    γ = lattice_constants_.γ

    return cos(α) * cos(β) * cos(γ) >= 0 || abs(cos(α) * cos(β) * cos(γ)) < COS_APPROX_ZERO
end

"""
    satisfies_triclinic_angle_constraints(α::Real, β::Real, γ::Real) -> Bool

Determine whether `α`, `β`, and `γ` satisfy the angle constraints for triclinic lattices:

* ``0 <  α + β + γ < 2π``
* ``0 <  α + β - γ < 2π``
* ``0 <  α - β + γ < 2π``
* ``0 < -α + β + γ < 2π``

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
