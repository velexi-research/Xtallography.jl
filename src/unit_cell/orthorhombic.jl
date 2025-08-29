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
Orthorhombic  unit cell types and functions
"""
# --- Exports

# Types
export OrthorhombicUnitCell, OrthorhombicUnitCellDelta

# --- Types

# ------ OrthorhombicUnitCell

"""
    OrthorhombicUnitCell

Lattice constants and symmetry for an orthorhombic unit cell

Fields
======
* `a`, `b`, `c`: lengths of the edges of the unit cell

* `symmetry`: unit cell symmetry
"""
const OrthorhombicUnitCell = UnitCell{Orthorhombic}

# Outer constructor
"""
    OrthorhombicUnitCell(
        a::Real, b::Real, c::Real;
        centering::Centering=primitive_centering,
        symmetry_elements::Union{Set,Vector,Nothing}=nothing
    )

Construct a OrthorhombicUnitCell object from a set of lattice constants.

Keyword Arguments
=================
- `centering`: centering of unit cell

- `symmetry_elements`: symmetry elements of crystal
"""
function OrthorhombicUnitCell(
    a::Real,
    b::Real,
    c::Real;
    centering::Centering=primitive_centering,
    symmetry_elements::Union{Set,Vector,Nothing}=nothing,
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

    # symmetry elements
    # TODO

    # --- Construct and return OrthorhombicUnitCell object

    return OrthorhombicUnitCell(
        (a=a, b=b, c=c); centering=centering, symmetry_elements=symmetry_elements
    )
end

# ------ OrthorhombicUnitCellDelta

"""
    OrthorhombicUnitCellDelta

Lattice constant deltas for an orthorhombic unit cell

Fields
======
* `Δa`, `Δb`, `Δc`: deltas of the lengths of the edges of the unit cell
"""
const OrthorhombicUnitCellDelta = UnitCellDelta{Orthorhombic}

# Outer constructors
"""
    OrthorhombicUnitCellDelta(Δa::Real, Δb::Real, Δc::Real)

Construct a OrthorhombicUnitCellDelta object from a set of lattice constant deltas.
"""
function OrthorhombicUnitCellDelta(Δa::Real, Δb::Real, Δc::Real)
    Δlattice_constants = (Δa=Δa, Δb=Δb, Δc=Δc)
    return OrthorhombicUnitCellDelta(Δlattice_constants)
end

# --- Functions/Methods

# ------ UnitCell methods

function standardize(unit_cell::OrthorhombicUnitCell)
    # --- Check arguments

    standardize_check_args(unit_cell)

    # --- Preparations

    # Extract lattice constants
    lattice_constants_ = lattice_constants(unit_cell)
    a = lattice_constants_.a
    b = lattice_constants_.b
    c = lattice_constants_.c

    # Extract centering
    centering_ = centering(unit_cell)

    # --- Standardize lattice constants

    if centering_ === base_centering
        return OrthorhombicUnitCell(sort([a, b])..., c; centering=centering_)
    end

    # all other centerings
    return OrthorhombicUnitCell(sort([a, b, c])...; centering=centering_)
end

function basis(unit_cell::OrthorhombicUnitCell)
    # Get lattice constants
    lattice_constants_ = lattice_constants(unit_cell)

    # Construct basis
    basis_a = Vector{Float64}([lattice_constants_.a, 0, 0])
    basis_b = Vector{Float64}([0, lattice_constants_.b, 0])
    basis_c = Vector{Float64}([0, 0, lattice_constants_.c])

    return basis_a, basis_b, basis_c
end

function volume(unit_cell::OrthorhombicUnitCell)
    # Get lattice constants
    lattice_constants_ = lattice_constants(unit_cell)

    # Compute volume
    return lattice_constants_.a * lattice_constants_.b * lattice_constants_.c
end

function surface_area(unit_cell::OrthorhombicUnitCell)
    # Get lattice constants
    lattice_constants_ = lattice_constants(unit_cell)
    a = lattice_constants_.a
    b = lattice_constants_.b
    c = lattice_constants_.c

    # Compute surface area
    return 2 * (a * b + b * c + c * a)
end

function conventional_cell(::Orthorhombic, unit_cell::UnitCell)
    # --- Check arguments

    conventional_cell_check_args(unit_cell)

    # --- Preparations

    # Standardize unit cell
    unit_cell = standardize(unit_cell)

    # Get lattice constants and centering
    lattice_constants_ = lattice_constants(unit_cell)
    a = lattice_constants.a
    b = lattice_constants.b
    c = lattice_constants.c

    centering_ = centering(unit_cell)

    # --- Compute IUCr conventional cell

    # Check limiting cases
    if centering_ === primitive_centering
        # Tetragonal, primitive-centering
        if a ≈ b
            @debug "oP --> tP"
            return conventional_cell(
                TetragonalUnitCell(a, c; centering=primitive_centering)
            )
        elseif b ≈ c
            @debug "oP --> tP"
            return conventional_cell(
                TetragonalUnitCell(c, a; centering=primitive_centering)
            )
        end

    elseif centering_ === body_centering
        # Tetragonal, body-centering
        if a ≈ b
            @debug "oI --> tI"
            return conventional_cell(TetragonalUnitCell(a, c; centering=body_centering))
        elseif b ≈ c
            @debug "oI --> tI"
            return conventional_cell(TetragonalUnitCell(c, a; centering=body_centering))
        end

    elseif centering_ === face_centering
        # Tetragonal, body-centering
        if a ≈ b
            @debug "oF --> tI"
            return conventional_cell(
                TetragonalUnitCell(a * SIN_PI_OVER_FOUR, c; centering=body_centering)
            )
        elseif b ≈ c
            @debug "oF --> tI"
            return conventional_cell(
                TetragonalUnitCell(c * SIN_PI_OVER_FOUR, a; centering=body_centering)
            )
        end

    elseif centering_ === base_centering
        if a ≈ b
            # Tetragonal, primitive-centering
            @debug "oC --> tP"
            return conventional_cell(
                TetragonalUnitCell(a * SIN_PI_OVER_FOUR, c; centering=primitive_centering)
            )

        elseif b ≈ 2 * a * SIN_PI_OVER_THREE
            # Hexagonal, primitive-centering
            @debug "oC --> hP"
            return conventional_cell(HexagonalUnitCell(a, c; centering=primitive_centering))
        end
    end

    # Not a limiting case, so return unit cell with standardized lattice constants
    return UnitCell(lattice_constants_; centering=centering_)
end
