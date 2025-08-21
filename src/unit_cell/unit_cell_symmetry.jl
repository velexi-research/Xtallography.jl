#   Copyright 2025 Velexi Corporation
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
UnitCellSymmetry type and functions
"""
# --- Exports

# ------ Types

export UnitCellSymmetry

# ------ Constants

export primitive_unit_cell_symmetry

# --- Types

using Xtallography: Centering, primitive_centering

"""
    UnitCellSymmetry

Type representing the symmetry of a unit cell
"""
struct UnitCellSymmetry
    # Fields
    centering::Centering
    symmetry_elements::Vector{<:SymmetryElement}
end

# Outer constructor
"""
    UnitCellSymmetry(
        centering::Centering;
        symmetry_elements::Union{Vector{<:SymmetryElement},Nothing}=nothing
    )


Construct a `UnitCellSymmetry` object with the specified `centering` and
`symmetry_elements`.
"""
function UnitCellSymmetry(
    centering::Centering; symmetry_elements::Union{Vector,Nothing}=nothing
)

    # Check arguments
    if isnothing(symmetry_elements)
        symmetry_elements = []
    end

    # Construct UnitCellSymmetry object
    return UnitCellSymmetry(centering, Vector{SymmetryElement}(symmetry_elements))
end

"""
    UnitCellSymmetry()

Construct a `UnitCellSymmetry` object representing a primitive unit cell with no
additional symmetry elements.
"""
function UnitCellSymmetry()
    return UnitCellSymmetry(primitive_centering)
end

"""
    primitive_unit_cell_symmetry

`UnitCellSymmetry` object representing a primitive unit cell with no additional symmetry
elements.
"""
const primitive_unit_cell_symmetry = UnitCellSymmetry()

# Methods
function Base.:(==)(x::UnitCellSymmetry, y::UnitCellSymmetry)
    return (
        x.centering === y.centering && Set(x.symmetry_elements) == Set(y.symmetry_elements)
    )
end
