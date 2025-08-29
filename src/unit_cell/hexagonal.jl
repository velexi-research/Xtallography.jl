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
Hexagonal unit cell types and functions
"""
# --- Exports

# Types
export HexagonalUnitCell, HexagonalUnitCellDelta

# --- Types

# ------ HexagonalUnitCell

"""
    HexagonalUnitCell

Lattice constants and symmetry for a hexagonal unit cell

Fields
======
* `a`, `c`: lengths of the edges of the unit cell

* `symmetry`: unit cell symmetry
"""
const HexagonalUnitCell = UnitCell{Hexagonal}

# Outer constructor
"""
    HexagonalUnitCell(
        a::Real, c::Real;
        centering::Centering=primitive_centering,
        symmetry_elements::Union{Set,Vector,Nothing}=nothing
    )

Construct a HexagonalUnitCell object from a set of lattice constants.

!!! note

    No constraints are imposed on `centering`. The unit cell does _not_ to be a valid
    Bravais lattice.

Keyword Arguments
=================
- `centering`: centering of unit cell

- `symmetry_elements`: symmetry elements of crystal
"""
function HexagonalUnitCell(
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

    # --- Construct and return HexagonalUnitCell object

    return HexagonalUnitCell(
        (a=a, c=c); centering=centering, symmetry_elements=symmetry_elements
    )
end

# ------ HexagonalUnitCellDelta

"""
    HexagonalUnitCellDelta

Lattice constant deltas for a hexagonal unit cell

Fields
======
* `Δa`, `Δc`: deltas of the lengths of the edges of the unit cell
"""
const HexagonalUnitCellDelta = UnitCellDelta{Hexagonal}

# Outer constructors
"""
    HexagonalUnitCellDelta(Δa::Real, Δc::Real)

Construct a HexagonalUnitCellDelta object from a set of lattice constant deltas.
"""
function HexagonalUnitCellDelta(Δa::Real, Δc::Real)
    Δlattice_constants = (Δa=Δa, Δc=Δc)
    return HexagonalUnitCellDelta(Δlattice_constants)
end

# --- Functions/Methods

# ------ UnitCell methods

function basis(unit_cell::HexagonalUnitCell)
    # Get lattice constants
    lattice_constants_ = lattice_constants(unit_cell)
    a = lattice_constants_.a
    c = lattice_constants_.c

    # Construct basis
    basis_a = Vector{Float64}([a, 0, 0])
    basis_b = Vector{Float64}([-0.5 * a, SIN_PI_OVER_THREE * a, 0])
    basis_c = Vector{Float64}([0, 0, c])

    return basis_a, basis_b, basis_c
end

function volume(unit_cell::HexagonalUnitCell)
    # Get lattice constants
    lattice_constants_ = lattice_constants(unit_cell)
    a = lattice_constants_.a
    c = lattice_constants_.c

    # Compute volume
    return a^2 * SIN_PI_OVER_THREE * c
end

function surface_area(unit_cell::HexagonalUnitCell)
    # Get lattice constants
    lattice_constants_ = lattice_constants(unit_cell)
    a = lattice_constants_.a
    c = lattice_constants_.c

    # Compute surface area
    return 4 * a * c + 2 * a^2 * SIN_PI_OVER_THREE
end
