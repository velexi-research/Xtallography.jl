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
Tests for unit_cell/unit_cell.jl
"""
# --- Imports

# Standard library
using LinearAlgebra: qr, I
using Test

# External packages
using Combinatorics: permutations

# Xtallography package
using Xtallography

# --- Tests

# Note: the following methods are tested in lattice-specific test suites
#
# * UnitCell{T}(::NamedTuple, ::UnitCellSymmetry)
# * UnitCell{T}(::NamedTuple; ::Centering, ::Union{Set,Vector,Nothing})
# * UnitCell(unit_cell::UnitCell{T})
# * UnitCell(::Vector{<:Real}, ::Vector{<:Real}, ::Vector{<:Real}; ::Bool, ::Centering)
#
# * lattice_system(::UnitCell)
# * lattice_constants(::UnitCell)
# * symmetry(::UnitCell)
# * centering(::UnitCell)
# * symmetry_elements(::UnitCell)
# * is_bravais_lattice(::UnitCell)
# * basis(::UnitCell)
# * volume(::UnitCell)
# * surface_area(::UnitCell)
#
# * standardize(::UnitCell)
# * conventional_cell(::UnitCell)
# * reduced_cell(::UnitCell)
#
# * is_equivalent_cell(::UnitCell)
# * is_supercell(::UnitCell)
# * isapprox(::UnitCell)
#
# * -(::UnitCell,::UnitCell)

# ------ Types

@testset "UnitCell Subtypes" begin
    expected_types = [
        TriclinicUnitCell,
        MonoclinicUnitCell,
        OrthorhombicUnitCell,
        TetragonalUnitCell,
        RhombohedralUnitCell,
        HexagonalUnitCell,
        CubicUnitCell,
    ]
    for type in expected_types
        @test type <: UnitCell
    end
end

# ------ Methods

@testset "reduced_cell(): minimum sum of length squared not unique" begin
    # --- Preparations

    unit_cell = TriclinicUnitCell(
        sqrt(6),
        sqrt(8),
        sqrt(8),
        π / 3,
        acos(sqrt(3) / 6),
        acos(sqrt(3) / 4);
        centering=primitive_centering,
    )

    # --- Exercise functionality

    reduced_cell_ = reduced_cell(unit_cell)

    # --- Check results

    # Check lattice constants
    lattice_constants_ = lattice_constants(reduced_cell_)
    @test isapprox(lattice_constants_.a, 2.449; atol=0.0005)
    @test isapprox(lattice_constants_.b, 2.828; atol=0.0005)
    @test isapprox(lattice_constants_.c, 2.828; atol=0.0005)
    @test isapprox(lattice_constants_.α, 104.47 * π / 180; atol=0.0005)
    @test isapprox(lattice_constants_.β, 106.78 * π / 180; atol=0.0005)
    @test isapprox(lattice_constants_.γ, 115.66 * π / 180; atol=0.0005)
    @test centering(reduced_cell_) === primitive_centering
end

@testset "is_equivalent_unit_cell(::UnitCell): valid arguments" begin
    # --- Tests

    # equivalent unit cells with the same lattice system, default atol and rtol
    a = 1
    b = 2
    c = 3
    β = 3π / 5
    unit_cell_ref = MonoclinicUnitCell(a, b, c, β)

    c_alt = sqrt(a^2 + c^2 + 2 * a * c * cos(β))
    β_alt = π - asin(sin(β) / c_alt * c)
    unit_cell_test = MonoclinicUnitCell(a, b, c_alt, β_alt)

    @test is_equivalent_unit_cell(unit_cell_test, unit_cell_ref)

    # equivalent unit cells with different lattice systems, default atol and rtol
    b = 2
    c = 3
    β = 3π / 5
    a = -2 * c * cos(β)
    unit_cell_ref = MonoclinicUnitCell(a, b, c, β; centering=primitive_centering)
    unit_cell_test = OrthorhombicUnitCell(a, 2 * c * sin(β), b; centering=base_centering)

    @test is_equivalent_unit_cell(unit_cell_test, unit_cell_ref)

    # nonequivalent unit cells that are equivalent when atol is sufficiently large
    a = 1
    b = 2
    c = 3
    unit_cell_ref = OrthorhombicUnitCell(a, b, c; centering=primitive_centering)
    unit_cell_test = OrthorhombicUnitCell(
        a + 2, b - 1, c + 5; centering=primitive_centering
    )
    @test !is_equivalent_unit_cell(unit_cell_test, unit_cell_ref)
    @test !is_equivalent_unit_cell(unit_cell_test, unit_cell_ref; atol=4)
    @test is_equivalent_unit_cell(unit_cell_test, unit_cell_ref; atol=6)

    # nonequivalent unit cells that are equivalent when rtol is sufficiently large
    a = 1
    b = 2
    c = 3
    unit_cell_ref = OrthorhombicUnitCell(a, b, c; centering=primitive_centering)
    unit_cell_test = OrthorhombicUnitCell(
        a + 2, b - 1, c + 5; centering=primitive_centering
    )
    @test !is_equivalent_unit_cell(unit_cell_test, unit_cell_ref)
    @test !is_equivalent_unit_cell(unit_cell_test, unit_cell_ref; atol=0, rtol=0.5)
    @test is_equivalent_unit_cell(unit_cell_test, unit_cell_ref; atol=0, rtol=1)

    # nonequivalent unit cells, case #1
    b = 2
    c = 3
    β = 3π / 5
    a = -2 * c * cos(β)
    unit_cell_ref = MonoclinicUnitCell(a, b, c, β; centering=primitive_centering)
    unit_cell_test = CubicUnitCell(a; centering=face_centering)

    @test !is_equivalent_unit_cell(unit_cell_test, unit_cell_ref)

    # nonequivalent unit cells, case #2
    a = 1
    unit_cell_ref = CubicUnitCell(a; centering=primitive_centering)
    unit_cell_test = CubicUnitCell(a; centering=face_centering)

    @test !is_equivalent_unit_cell(unit_cell_test, unit_cell_ref)
end

@testset "is_equivalent_unit_cell(::UnitCell): invalid arguments" begin
    # --- Preparations

    unit_cell_ref = CubicUnitCell(1.0)
    unit_cell_test = CubicUnitCell(1.0)

    # --- Exercise functionality and check results

    # atol < 0
    expected_message = "`atol` must be nonnegative"
    @test_throws DomainError(-1, expected_message) is_equivalent_unit_cell(
        unit_cell_test, unit_cell_ref; atol=-1
    )

    # rtol < 0
    expected_message = "`rtol` must be nonnegative"
    @test_throws DomainError(-1, expected_message) is_equivalent_unit_cell(
        unit_cell_test, unit_cell_ref; rtol=-1
    )
end

@testset "is_supercell(): different UnitCell types" begin
    # --- Tests

    cubic_unit_cell = CubicUnitCell(1.0)
    orthorhombic_unit_cell = OrthorhombicUnitCell(1.0, 2.0, 3.0)

    @test !is_supercell(cubic_unit_cell, orthorhombic_unit_cell)
end

@testset "isapprox(::UnitCell): comparison between different types" begin
    # x ≉ y
    x = CubicUnitCell(1)
    y = TetragonalUnitCell(1, 2)
    @test x ≉ y
end
