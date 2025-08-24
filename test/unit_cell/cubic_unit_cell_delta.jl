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
Tests for CubicUnitCellDelta types and methods in unit_cell/cubic.jl
"""
# --- Imports

# Standard library
using Test

# Xtallography package
using Xtallography

# --- Tests

# ------ Constructors

@testset "CubicUnitCellDelta outer constructor" begin
    # --- Tests

    Δa = 1
    unit_cell_delta = CubicUnitCellDelta(Δa)

    @test Δlattice_constants(unit_cell_delta).Δa == Δa
end

# ------ Methods

@testset "lattice_system(::CubicUnitCellDelta)" begin
    # --- Tests

    Δlattice_constants = CubicUnitCellDelta(1)
    @test lattice_system(Δlattice_constants) === cubic
end

@testset "Δlattice_constants(::CubicUnitCellDelta)" begin
    # --- Tests

    lattice_constant_deltas = (Δa=1,)
    unit_cell_delta = CubicUnitCellDelta(lattice_constant_deltas)

    @test Δlattice_constants(unit_cell_delta) == lattice_constant_deltas
end

@testset "lattice_constant_deltas(::CubicUnitCellDelta)" begin
    # --- Tests

    lattice_constant_deltas_ = (Δa=1,)
    unit_cell_delta = CubicUnitCellDelta(lattice_constant_deltas_)

    @test lattice_constant_deltas(unit_cell_delta) == lattice_constant_deltas_
end

@testset "isapprox(::CubicUnitCellDelta)" begin
    # --- Preparations

    x = CubicUnitCellDelta(1.0)
    y = CubicUnitCellDelta(2.0)

    # --- Exercise functionality and check results

    # x ≈ (x + delta)
    @test x ≈ CubicUnitCellDelta(1 + 1e-9)

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
