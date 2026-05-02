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
Tests for HexagonalUnitCellDelta types and methods in unit_cell/hexagonal.jl
"""
# --- Imports

# Standard library
using Test

# Xtallography package
using Xtallography

# --- Tests

# ------ Constructors

@testset "HexagonalUnitCellDelta inner constructor" begin
    # --- Tests

    # Valid arguments
    Δlattice_constants_ = (Δa=1, Δc=5)
    unit_cell_delta = HexagonalUnitCellDelta(Δlattice_constants_)

    @test Δlattice_constants(unit_cell_delta) == Δlattice_constants_

    # Invalid arguments
    Δlattice_constants_ = (Δa=1,)
    expected_message = (
        "Invalid Δlattice_constants argument passed to UnitCellDelta{Hexagonal} " *
        "constructor. Expected keys: (:Δa, :Δc). " *
        "Provided keys: $(keys(Δlattice_constants_))."
    )
    @test_throws ArgumentError(expected_message) HexagonalUnitCellDelta(Δlattice_constants_)
end

@testset "HexagonalUnitCellDelta outer constructor" begin
    # --- Tests

    Δa = 1
    Δc = 5
    unit_cell_delta = HexagonalUnitCellDelta(Δa, Δc)

    @test Δlattice_constants(unit_cell_delta).Δa == Δa
    @test Δlattice_constants(unit_cell_delta).Δc == Δc
end

# ------ Methods

@testset "lattice_system(::HexagonalUnitCellDelta)" begin
    # --- Tests

    Δlattice_constants = HexagonalUnitCellDelta(1, 2)
    @test lattice_system(Δlattice_constants) === hexagonal
end

@testset "Δlattice_constants(::HexagonalUnitCellDelta)" begin
    # --- Tests

    lattice_constant_deltas = (Δa=1, Δc=3)
    unit_cell_delta = HexagonalUnitCellDelta(lattice_constant_deltas)

    @test Δlattice_constants(unit_cell_delta) == lattice_constant_deltas
end

@testset "lattice_constant_deltas(::HexagonalUnitCellDelta)" begin
    # --- Tests

    lattice_constant_deltas_ = (Δa=1, Δc=3)
    unit_cell_delta = HexagonalUnitCellDelta(lattice_constant_deltas_)

    @test lattice_constant_deltas(unit_cell_delta) == lattice_constant_deltas_
end

@testset "isapprox(::HexagonalUnitCellDelta)" begin
    # --- Preparations

    x = HexagonalUnitCellDelta(1.0, 2.0)
    y = HexagonalUnitCellDelta(1.5, 2.5)

    # --- Exercise functionality and check results

    # x ≈ (x + delta)
    @test x ≈ HexagonalUnitCellDelta(1.0 + 1e-9, 2.0)
    @test x ≈ HexagonalUnitCellDelta(1.0, 2.0 + 1e-9)

    # x !≈ y
    @test !(x ≈ y)

    # x ≈ y: atol = 1
    @test isapprox(x, y; atol=1)

    # x ≈ y: rtol = 1
    @test isapprox(x, y; rtol=1)

    # x ≈ y: atol = 0.01, rtol = 1
    @test isapprox(x, y; atol=0.01, rtol=1)

    # x ≈ y: atol = 1, rtol = 0.01
    @test isapprox(x, y; atol=1, rtol=0.01)

    # x !≈ y: atol = 0.01, rtol = 0.01
    @test !isapprox(x, y; atol=0.01, rtol=0.01)
end
