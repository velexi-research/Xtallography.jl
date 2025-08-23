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
Unit cell delta type and functions
"""
# --- Exports

# Types
export UnitCellDelta

# Functions
export isapprox

# --- Types

"""
    UnitCellDelta{T<:LatticeSystem}

Delta between two unit cells (for the same lattice system)

Type Aliases
============
[`TriclinicUnitCellDelta`](@ref), [`MonoclinicUnitCellDelta`](@ref),
[`OrthorhombicUnitCellDelta`](@ref), [`TetragonalUnitCellDelta`](@ref),
[`RhombohedralUnitCellDelta`](@ref), [`HexagonalUnitCellDelta`](@ref),
[`CubicUnitCellDelta`](@ref)
"""
#abstract type UnitCellDelta{T<:LatticeSystem} end
struct UnitCellDelta{T<:LatticeSystem}
    # Fields
    Δlattice_constants::NamedTuple

    # Constructor
    function UnitCellDelta{T}(Δlattice_constants::NamedTuple) where {T<:LatticeSystem}

        # --- Check arguments

        # TODO: check that lattice constants is consistent with LatticeSystem

        # --- Return new UnitCellDelta

        return new(Δlattice_constants)
    end
end

# --- Functions/Methods

import Base.isapprox

function lattice_system(unit_cell::UnitCellDelta{T}) where {T<:LatticeSystem}
    return T()
end

function isapprox(
    x::UnitCellDelta, y::UnitCellDelta; atol::Real=0, rtol::Real=atol > 0 ? 0 : √eps()
)
    # Default isapprox() implementation to allow comparison between unit cell deltas for
    # different lattice systems.
    return false
end

function isapprox(
    x::UnitCellDelta{T}, y::UnitCellDelta{T}; atol::Real=0, rtol::Real=atol > 0 ? 0 : √eps()
) where {T<:LatticeSystem}
    return (all(
        isapprox(
            getfield(x.Δlattice_constants, name),
            getfield(y.Δlattice_constants, name);
            atol=atol,
            rtol=rtol,
        ) for name in fieldnames(typeof(x.Δlattice_constants))
    ))
end
