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
Tetragonal unit cell types and functions
"""
# --- Exports

# Types
export TetragonalUnitCell, TetragonalUnitCellDelta

# --- Types

# ------ TetragonalUnitCell

"""
    TetragonalUnitCell

Lattice constants and symmetry for a tetragonal unit cell

Fields
======
* `a`, `c`: lengths of the edges of the unit cell

* `symmetry`: unit cell symmetry
"""
const TetragonalUnitCell = UnitCell{Tetragonal}

# Outer constructor
"""
    TetragonalUnitCell(
        a::Real, c::Real;
        centering::Centering=primitive_centering,
        symmetry_elements::Union{Set,Vector,Nothing}=nothing
    )

Construct a TetragonalUnitCell object from a set of lattice constants.

!!! note

    No constraints are imposed on `centering`. The unit cell does _not_ to be a valid
    Bravais lattice.

Keyword Arguments
=================
- `centering`: centering of unit cell

- `symmetry_elements`: symmetry elements of crystal
"""
function TetragonalUnitCell(
    a::Real,
    c::Real;
    centering::Centering=primitive_centering,
    symmetry_elements::Union{Set,Vector,Nothing}=nothing,
)

    # --- Check arguments

    # lattice constants
    if a <= 0
        throw(DomainError(a, "`a` must be positive"))
    end

    if c <= 0
        throw(DomainError(c, "`c` must be positive"))
    end

    # symmetry elements
    # TODO

    # --- Construct and return TetragonalUnitCell object

    return TetragonalUnitCell(
        (a=a, c=c); centering=centering, symmetry_elements=symmetry_elements
    )
end

# ------ TetragonalUnitCellDelta

"""
    TetragonalUnitCellDelta

Lattice constant deltas for a tetragonal unit cell

Fields
======
* `Δa`, `Δc`: deltas of the lengths of the edges of the unit cell
"""
const TetragonalUnitCellDelta = UnitCellDelta{Tetragonal}

# Outer constructors
"""
    TetragonalUnitCellDelta(Δa::Real, Δc::Real)

Construct a TetragonalUnitCellDelta object from a set of lattice constant deltas.
"""
function TetragonalUnitCellDelta(Δa::Real, Δc::Real)
    Δlattice_constants = (Δa=Δa, Δc=Δc)
    return TetragonalUnitCellDelta(Δlattice_constants)
end

# --- Functions/Methods

# ------ UnitCell methods

function basis(unit_cell::TetragonalUnitCell)
    # Construct basis
    lattice_constants_ = lattice_constants(unit_cell)
    basis_a = Vector{Float64}([lattice_constants_.a, 0, 0])
    basis_b = Vector{Float64}([0, lattice_constants_.a, 0])
    basis_c = Vector{Float64}([0, 0, lattice_constants_.c])

    return basis_a, basis_b, basis_c
end

function volume(unit_cell::TetragonalUnitCell)
    lattice_constants_ = lattice_constants(unit_cell)
    return lattice_constants_.a^2 * lattice_constants_.c
end

function surface_area(unit_cell::TetragonalUnitCell)
    lattice_constants_ = lattice_constants(unit_cell)
    return 2 * lattice_constants_.a^2 + 4 * lattice_constants_.a * lattice_constants_.c
end

function conventional_cell(::Tetragonal, unit_cell::UnitCell)
    # --- Check arguments

    conventional_cell_check_args(unit_cell)

    # --- Preparations

    # Get lattice constants and centering
    lattice_constants_ = lattice_constants(unit_cell)
    a = lattice_constants_.a
    c = lattice_constants_.c

    centering_ = centering(unit_cell)

    # --- Compute IUCr conventional cell

    # Check limiting cases
    if centering_ === primitive_centering
        if a ≈ c
            # cubic, primitive-centering, edge length `a`
            @debug "tP --> cP"
            return conventional_cell(CubicUnitCell(a; centering=primitive_centering))
        end

    elseif centering_ === body_centering
        if a ≈ c
            # cubic, body-centering, edge length `a`
            @debug "tI --> cI"
            return conventional_cell(CubicUnitCell(a; centering=body_centering))

        elseif c * SIN_PI_OVER_FOUR ≈ a
            # cubic, face-centering, edge length `c`
            @debug "tI --> cF"
            return conventional_cell(CubicUnitCell(c; centering=face_centering))
        end
    end

    # Not a limiting case, so return unit cell with original lattice constants
    return unit_cell
end
