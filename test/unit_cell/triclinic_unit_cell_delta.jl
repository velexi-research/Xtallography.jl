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
Tests for TriclinicUnitCellDelta types and methods in unit_cell/triclinic.jl
"""
# --- Imports

# Standard library
using Test

# Xtallography package
using Xtallography

# --- Tests

# ------ Constructors

@testset "TriclinicUnitCellDelta inner constructor" begin
    # --- Tests

    # Valid arguments
    Δlattice_constants_ = (Δa=1, Δb=3, Δc=5, Δα=π / 7, Δβ=2π / 7, Δγ=3π / 7)
    unit_cell_delta = TriclinicUnitCellDelta(Δlattice_constants_)

    @test Δlattice_constants(unit_cell_delta) == Δlattice_constants_

    # Invalid arguments
    Δlattice_constants_ = (Δa=1,)
    expected_message = (
        "Invalid Δlattice_constants argument passed to UnitCellDelta{Triclinic} " *
        "constructor. Expected keys: (:Δa, :Δb, :Δc, :Δα, :Δβ, :Δγ). " *
        "Provided keys: $(keys(Δlattice_constants_))."
    )
    @test_throws ArgumentError(expected_message) TriclinicUnitCellDelta(Δlattice_constants_)
end

@testset "TriclinicUnitCellDelta outer constructor" begin
    # --- Tests

    Δa = 1
    Δb = 3
    Δc = 5
    Δα = π / 7
    Δβ = 2π / 7
    Δγ = 3π / 7
    unit_cell_delta = TriclinicUnitCellDelta(Δa, Δb, Δc, Δα, Δβ, Δγ)

    @test Δlattice_constants(unit_cell_delta).Δa == Δa
    @test Δlattice_constants(unit_cell_delta).Δb == Δb
    @test Δlattice_constants(unit_cell_delta).Δc == Δc
    @test Δlattice_constants(unit_cell_delta).Δα == Δα
    @test Δlattice_constants(unit_cell_delta).Δβ == Δβ
    @test Δlattice_constants(unit_cell_delta).Δγ == Δγ
end

# ------ Methods

@testset "lattice_system(::TriclinicConstantDelta)" begin
    # --- Tests

    Δlattice_constants = TriclinicUnitCellDelta(1, 2, 3, π / 5, 2π / 5, 3π / 5)
    @test lattice_system(Δlattice_constants) === triclinic
end

@testset "Δlattice_constants(::TriclinicUnitCellDelta)" begin
    # --- Tests

    lattice_constant_deltas = (Δa=1, Δb=2, Δc=3, Δα=π / 5, Δβ=2π / 5, Δγ=3π / 5)
    unit_cell_delta = TriclinicUnitCellDelta(lattice_constant_deltas)

    @test Δlattice_constants(unit_cell_delta) == lattice_constant_deltas
end

@testset "lattice_constant_deltas(::TriclinicUnitCellDelta)" begin
    # --- Tests

    lattice_constant_deltas_ = (Δa=1, Δb=2, Δc=3, Δα=π / 5, Δβ=2π / 5, Δγ=3π / 5)
    unit_cell_delta = TriclinicUnitCellDelta(lattice_constant_deltas_)

    @test lattice_constant_deltas(unit_cell_delta) == lattice_constant_deltas_
end

@testset "isapprox(::TriclinicUnitCellDelta)" begin
    # --- Preparations

    x = TriclinicUnitCellDelta(1.0, 2.0, 3.0, π / 5, π / 4, 2π / 5)
    y = TriclinicUnitCellDelta(1.5, 2.5, 3.5, π / 5 + 0.5, π / 4 + 0.5, 2π / 5 + 0.5)

    # --- Exercise functionality and check results

    # x ≈ (x + delta)
    @test x ≈ TriclinicUnitCellDelta(1.0 + 1e-9, 2.0, 3.0, π / 5, π / 4, 2π / 5)
    @test x ≈ TriclinicUnitCellDelta(1.0, 2.0 + 1e-9, 3.0, π / 5, π / 4, 2π / 5)
    @test x ≈ TriclinicUnitCellDelta(1.0, 2.0, 3.0 - 1e-9, π / 5, π / 4, 2π / 5)
    @test x ≈ TriclinicUnitCellDelta(1.0, 2.0, 3.0, π / 5 - 1e-9, π / 4, 2π / 5)
    @test x ≈ TriclinicUnitCellDelta(1.0, 2.0, 3.0, π / 5, π / 4 + 1e-9, 2π / 5)
    @test x ≈ TriclinicUnitCellDelta(1.0, 2.0, 3.0, π / 5, π / 4, 2π / 5 - 1e-9)

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
