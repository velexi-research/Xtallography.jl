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
Tests for unit_cell/unit_cell_delta.jl
"""
# --- Imports

# Standard library
using Test

# Xtallography package
using Xtallography

# --- Tests

# Note: the following methods are tested in lattice-specific test suites
#
# * UnitCellDelta{T}(::NamedTuple)
# * lattice_system(::UnitCellDelta)
# * Δlattice_constants(::UnitCellDelta)
# * lattice_constants_deltas(::UnitCellDelta)
# * isapprox(::UnitCell)

@testset "isapprox(::UnitCellDelta): comparison between different types" begin
    # --- Tests

    # x ≈ y
    x = CubicUnitCellDelta(1 + 1e-9)
    y = CubicUnitCellDelta(1)
    @test x ≈ y

    # x ≉ y
    x = CubicUnitCellDelta(1)
    y = TetragonalUnitCellDelta(1, 2)
    @test x ≉ y
end
