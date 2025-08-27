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
Tests for OrthorhobmicUnitCellDelta types and methods in unit_cell/orthorhombic.jl
"""
# --- Imports

# Standard library
using Test

# Xtallography package
using Xtallography

# --- Tests

# ------ Constructors

@testset "OrthorhombicUnitCellDelta inner constructor" begin
    # --- Tests

    # Valid arguments
    Δlattice_constants_ = (Δa=1, Δb=3, Δc=5)
    unit_cell_delta = OrthorhombicUnitCellDelta(Δlattice_constants_)

    @test Δlattice_constants(unit_cell_delta) == Δlattice_constants_

    # Invalid arguments
    Δlattice_constants_ = (Δa=1,)
    expected_message = (
        "Invalid Δlattice_constants argument passed to UnitCellDelta{Orthorhombic} " *
        "constructor. Expected keys: (:Δa, :Δb, :Δc). " *
        "Provided keys: $(keys(Δlattice_constants_))."
    )
    @test_throws ArgumentError(expected_message) OrthorhombicUnitCellDelta(
        Δlattice_constants_
    )
end

@testset "OrthorhombicUnitCellDelta outer constructor" begin
    # --- Tests

    Δa = 1
    Δb = 3
    Δc = 5
    unit_cell_delta = OrthorhombicUnitCellDelta(Δa, Δb, Δc)

    @test Δlattice_constants(unit_cell_delta).Δa == Δa
    @test Δlattice_constants(unit_cell_delta).Δb == Δb
    @test Δlattice_constants(unit_cell_delta).Δc == Δc
end

# ------ Methods

@testset "lattice_system(::OrthorhombicUnitCellDelta)" begin
    # --- Tests

    Δlattice_constants = OrthorhombicUnitCellDelta(1, 2, 3)
    @test lattice_system(Δlattice_constants) === orthorhombic
end

@testset "Δlattice_constants(::OrthorhombicUnitCellDelta)" begin
    # --- Tests

    lattice_constant_deltas = (Δa=1, Δb=2, Δc=3)
    unit_cell_delta = OrthorhombicUnitCellDelta(lattice_constant_deltas)

    @test Δlattice_constants(unit_cell_delta) == lattice_constant_deltas
end

@testset "lattice_constant_deltas(::OrthorhombicUnitCellDelta)" begin
    # --- Tests

    lattice_constant_deltas_ = (Δa=1, Δb=2, Δc=3)
    unit_cell_delta = OrthorhombicUnitCellDelta(lattice_constant_deltas_)

    @test lattice_constant_deltas(unit_cell_delta) == lattice_constant_deltas_
end

@testset "isapprox(::OrthorhombicUnitCellDelta)" begin
    # --- Preparations

    x = OrthorhombicUnitCellDelta(1.0, 2.0, 3.0)
    y = OrthorhombicUnitCellDelta(1.5, 2.5, 3.5)

    # --- Exercise functionality and check results

    # x ≈ (x + delta)
    @test x ≈ OrthorhombicUnitCellDelta(1.0 + 1e-9, 2.0, 3.0)
    @test x ≈ OrthorhombicUnitCellDelta(1.0, 2.0 + 1e-9, 3.0)
    @test x ≈ OrthorhombicUnitCellDelta(1.0, 2.0, 3.0 - 1e-9)

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
