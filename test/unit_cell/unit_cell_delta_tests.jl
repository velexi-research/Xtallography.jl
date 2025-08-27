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

@testset "lattice_system(::UnitCellDelta)" begin
    for lattice_system_ in
        (Triclinic, Monoclinic, Orthorhombic, Tetragonal, Hexagonal, Rhombohedral, Cubic)
        lattice_constants = (a=3,)
        unit_cell_delta = UnitCellDelta{lattice_system_}(lattice_constants)

        @test lattice_system(unit_cell_delta) === lattice_system_()
    end
end

@testset "Δlattice_constants(::UnitCellDelta)" begin
    for lattice_system_ in
        (Triclinic, Monoclinic, Orthorhombic, Tetragonal, Hexagonal, Rhombohedral, Cubic)
        Δlattice_constants_ = (a=3,)
        unit_cell_delta = UnitCellDelta{lattice_system_}(Δlattice_constants_)

        @test Δlattice_constants(unit_cell_delta) === Δlattice_constants_
    end
end

@testset "lattice_constant_deltas(::UnitCellDelta)" begin
    for lattice_system_ in
        (Triclinic, Monoclinic, Orthorhombic, Tetragonal, Hexagonal, Rhombohedral, Cubic)
        Δlattice_constants_ = (a=3,)
        unit_cell_delta = UnitCellDelta{lattice_system_}(Δlattice_constants_)

        @test lattice_constant_deltas(unit_cell_delta) === Δlattice_constants_
    end
end

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

# --- Lattice-specific tests
#
# * constructors
# * lattice_system(::UnitCellDelta)
# * Δlattice_constants(::UnitCellDelta)
# * lattice_constants_deltas(::UnitCellDelta)
# * isapprox(::UnitCell)
