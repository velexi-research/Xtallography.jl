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
Cubic unit cell types and functions
"""
# --- Exports

# Types
export CubicUnitCell, CubicUnitCellDelta

# --- Types

# ------ CubicUnitCell

"""
    CubicUnitCell

Lattice constant and symmetry for a cubic unit cell

Fields
======
* `a`: length of the edge of the unit cell

* `symmetry`: unit cell symmetry
"""
const CubicUnitCell = UnitCell{Cubic}

# Outer constructor
"""
    CubicUnitCell(
        a::Real;
        centering::Centering=primitive_centering,
        symmetry_elements::Union{Set,Vector,Nothing}=nothing
    )

Construct a CubicUnitCell object from a set of lattice constants.

!!! note

    No constraints are imposed on `centering`. The unit cell does _not_ to be a valid
    Bravais lattice.

Keyword Arguments
=================
- `centering`: centering of unit cell

- `symmetry_elements`: symmetry elements of crystal
"""
function CubicUnitCell(
    a::Real;
    centering::Centering=primitive_centering,
    symmetry_elements::Union{Set,Vector,Nothing}=nothing,
)

    # --- Check arguments

    # lattice constants
    if a <= 0
        throw(DomainError(a, "`a` must be positive"))
    end

    # symmetry elements
    # TODO

    # --- Construct and return new CubicUnitCell

    return CubicUnitCell((a=a,); centering=centering, symmetry_elements=symmetry_elements)
end

# ------ CubicUnitCellDelta

"""
    CubicUnitCellDelta

Lattice constant delta for a cubic unit cell

Fields
======
* `Δa`: delta of the length of the edge of the unit cell
"""
const CubicUnitCellDelta = UnitCellDelta{Cubic}

# Outer constructor
"""
    CubicUnitCellDelta(Δa::Real)

Construct a CubicUnitCellDelta object from a lattice constant delta.
"""
function CubicUnitCellDelta(Δa::Real)
    Δlattice_constants = (Δa=Δa,)
    return CubicUnitCellDelta(Δlattice_constants)
end

# --- Functions/Methods

# ------ UnitCell methods

function basis(unit_cell::CubicUnitCell)
    lattice_constants_ = lattice_constants(unit_cell)
    basis_a = Vector{Float64}([lattice_constants_.a, 0, 0])
    basis_b = Vector{Float64}([0, lattice_constants_.a, 0])
    basis_c = Vector{Float64}([0, 0, lattice_constants_.a])

    return basis_a, basis_b, basis_c
end

function volume(unit_cell::CubicUnitCell)
    return lattice_constants(unit_cell).a^3
end

function surface_area(unit_cell::CubicUnitCell)
    return 6 * lattice_constants(unit_cell).a^2
end

function is_supercell(
    unit_cell_test::CubicUnitCell,
    unit_cell_ref::CubicUnitCell;
    tol::Real=1e-3,
    max_index::Integer=3,
)
    # --- Check arguments

    if tol <= 0
        throw(DomainError(tol, "`tol` must be positive"))
    end

    # --- Compare lattice constants

    # Check that a_test is an integer multiple of a_ref
    multiplier = lattice_constants(unit_cell_test).a / lattice_constants(unit_cell_ref).a
    diff_from_int = abs(multiplier - round(multiplier))
    return diff_from_int < tol && !isapprox(multiplier, 1; atol=tol)
end
