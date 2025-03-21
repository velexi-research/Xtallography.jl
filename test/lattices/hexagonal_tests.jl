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
Tests for methods in lattice/hexagonal.jl (except for cell standardization methods)
"""
# --- Imports

# Standard library
using Test
using LinearAlgebra: det, norm, cross

# Xtallography package
using Xtallography

# --- Tests

# ------ Types

@testset "HexagonalLatticeConstants constructor: valid arguments" begin
    # --- Preparations

    a = 1
    c = 3

    # --- Exercise functionality

    lattice_constants = HexagonalLatticeConstants(a, c)

    # --- Check results

    @test lattice_constants.a == a
    @test lattice_constants.c == c
end

@testset "HexagonalLatticeConstants constructor: invalid arguments" begin
    # --- Preparations

    # Valid arguments
    a = 1
    c = 3

    # --- Tests

    # ------ a

    # a = 0
    expected_message = "`a` must be positive"
    @test_throws DomainError(0, expected_message) HexagonalLatticeConstants(0, c)

    # a < 0
    expected_message = "`a` must be positive"
    @test_throws DomainError(-1, expected_message) HexagonalLatticeConstants(-1.0, c)

    # ------ c

    # c = 0
    expected_message = "`c` must be positive"
    @test_throws DomainError(0, expected_message) HexagonalLatticeConstants(a, 0)

    # c < 0
    expected_message = "`c` must be positive"
    @test_throws DomainError(-1.0, expected_message) HexagonalLatticeConstants(a, -1.0)
end

@testset "HexagonalLatticeConstantDeltas constructor" begin
    # --- Tests

    Δa = 1
    Δc = 2
    Δlattice_constants = HexagonalLatticeConstantDeltas(Δa, Δc)

    @test Δlattice_constants.Δa == Δa
    @test Δlattice_constants.Δc == Δc
end

# ------ LatticeConstants functions

@testset "isapprox(::HexagonalLatticeConstants)" begin
    # --- Preparations

    x = HexagonalLatticeConstants(1.0, 2.0)
    y = HexagonalLatticeConstants(1.5, 2.5)

    # --- Exercise functionality and check results

    # x ≈ (x + delta)
    @test x ≈ HexagonalLatticeConstants(1.0 + 1e-9, 2.0)
    @test x ≈ HexagonalLatticeConstants(1.0, 2.0 + 1e-9)

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

@testset "convert(::HexagonalLatticeConstants)" begin
    # --- Tests

    x = HexagonalLatticeConstants((rand(2) .+ 0.1)...)

    x_vector = convert(Vector, x)
    @test x_vector == [x.a, x.c]
end

@testset "-(::HexagonalLatticeConstants)" begin
    # --- Tests

    x = HexagonalLatticeConstants(1, 3)
    y = HexagonalLatticeConstants(2, 10)
    @test x - y == HexagonalLatticeConstantDeltas(x.a - y.a, x.c - y.c)
end

@testset "lattice_system(::HexagonalLatticeConstants)" begin
    # --- Tests

    lattice_constants = HexagonalLatticeConstants(1, 2)
    @test lattice_system(lattice_constants) === hexagonal
end

@testset "standardize(): hexagonal" begin
    # --- Tests

    # ------ Hexagonal lattices have no lattice constants conventions for primitive
    #        centering

    lattice_constants = HexagonalLatticeConstants(1.0, 2.0)

    standardized_lattice_constants, standardized_centering = standardize(
        lattice_constants, primitive
    )

    expected_lattice_constants = lattice_constants
    @test standardized_lattice_constants ≈ expected_lattice_constants

    @test standardized_centering === primitive

    # ------ Invalid centerings

    for centering in (body_centered, face_centered, base_centered)
        expected_message =
            "Invalid Bravais lattice: " *
            "(lattice_system=Hexagonal, centering=$(nameof(typeof(centering))))"
        @test_throws ArgumentError(expected_message) standardize(
            lattice_constants, centering
        )
    end
end

# ------ LatticeConstantDeltas functions

@testset "isapprox(::HexagonalLatticeConstantDeltas)" begin
    # --- Preparations

    x = HexagonalLatticeConstantDeltas(1.0, 2.0)
    y = HexagonalLatticeConstantDeltas(1.5, 2.5)

    # --- Exercise functionality and check results

    # x ≈ (x + delta)
    @test x ≈ HexagonalLatticeConstantDeltas(1.0 + 1e-9, 2.0)
    @test x ≈ HexagonalLatticeConstantDeltas(1.0, 2.0 + 1e-9)

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

@testset "convert(::HexagonalLatticeConstantDeltas)" begin
    # --- Tests

    x = HexagonalLatticeConstantDeltas(rand(2)...)

    x_vector = convert(Vector, x)
    @test x_vector == [x.Δa, x.Δc]
end

@testset "lattice_system(::HexagonalLatticeConstantDeltas)" begin
    # --- Tests

    Δlattice_constants = HexagonalLatticeConstantDeltas(1, 2)
    @test lattice_system(Δlattice_constants) === hexagonal
end

# ------ Unit cell computations

@testset "basis(::HexagonalLatticeConstants)" begin
    # --- Preparations

    a = 2
    c = 5
    lattice_constants = HexagonalLatticeConstants(a, c)

    # --- Exercise functionality

    basis_a, basis_b, basis_c = basis(lattice_constants)

    # --- Check results

    @test basis_a ≈ [a, 0, 0]
    @test basis_b ≈ [-0.5 * a, a * sqrt(3) / 2, 0]
    @test basis_c ≈ [0, 0, c]
end

@testset "volume(::HexagonalLatticeConstants)" begin
    # --- Preparations

    # Construct basis vectors for unit cell
    a = 6
    c = 8
    lattice_constants = HexagonalLatticeConstants(a, c)

    basis_a, basis_b, basis_c = basis(lattice_constants)

    # --- Exercise functionality and check results

    @test volume(lattice_constants) ≈ abs(det(hcat(basis_a, basis_b, basis_c)))
end

@testset "surface_area(::HexagonalLatticeConstants)" begin
    # --- Preparations

    # Construct basis vectors for unit cell
    a = 6
    c = 8
    lattice_constants = HexagonalLatticeConstants(a, c)

    basis_a, basis_b, basis_c = basis(lattice_constants)

    # --- Exercise functionality and check results

    @test surface_area(lattice_constants) ≈
        2 * norm(cross(basis_a, basis_b)) +
          2 * norm(cross(basis_b, basis_c)) +
          2 * norm(cross(basis_c, basis_a))
end

@testset "reduced_cell(): hexagonal" begin
    # --- Preparations

    a = 2
    c = 5
    lattice_constants = HexagonalLatticeConstants(a, c)
    basis_a, basis_b, basis_c = basis(lattice_constants)

    # --- Exercise functionality and check results

    # primitive unit cell defined by [basis_a, basis_b, basis_c]
    unit_cell = UnitCell(lattice_constants, primitive)

    expected_reduced_cell = reduced_cell(
        UnitCell(
            LatticeConstants(basis_a, basis_b, basis_c; identify_lattice_system=false),
            primitive,
        ),
    )

    reduced_cell_ = reduced_cell(unit_cell)
    @test reduced_cell_.lattice_constants isa HexagonalLatticeConstants
    @test volume(reduced_cell_) ≈ volume(unit_cell)
    @test reduced_cell_ ≈ expected_reduced_cell

    # primitive unit cell defined by linear combination of [basis_a, basis_b, basis_c],
    # β ≈ π / 3
    unit_cell = UnitCell(LatticeConstants(basis_a + basis_b, basis_b, basis_c), primitive)

    expected_reduced_cell = reduced_cell(UnitCell(lattice_constants, primitive))

    reduced_cell_ = reduced_cell(unit_cell)
    @test reduced_cell_.lattice_constants isa HexagonalLatticeConstants
    @test volume(reduced_cell_) ≈ volume(unit_cell)
    @test reduced_cell_ ≈ expected_reduced_cell

    # primitive unit cell defined by linear combination of [basis_a, basis_b, basis_c],
    # β ≈ 2π / 3
    unit_cell = UnitCell(LatticeConstants(basis_a - basis_b, basis_b, basis_c), primitive)

    expected_reduced_cell = reduced_cell(UnitCell(lattice_constants, primitive))

    reduced_cell_ = reduced_cell(unit_cell)
    @test reduced_cell_.lattice_constants isa HexagonalLatticeConstants
    @test volume(reduced_cell_) ≈ volume(unit_cell)
    @test reduced_cell_ ≈ expected_reduced_cell
end

@testset "is_equivalent_unit_cell(::UnitCell): hexagonal" begin
    # --- Preparations

    a = 2
    c = 5
    lattice_constants = HexagonalLatticeConstants(a, c)
    basis_a, basis_b, basis_c = basis(lattice_constants)

    # --- Exercise functionality and check results

    # equivalent hexagonal and triclinic unit cells
    hexagonal_unit_cell = UnitCell(lattice_constants, primitive)
    triclinic_unit_cell = UnitCell(
        LatticeConstants(basis_a, basis_b, basis_c; identify_lattice_system=false),
        primitive,
    )
    @test is_equivalent_unit_cell(hexagonal_unit_cell, triclinic_unit_cell)

    # equivalent unit cell defined by linear combination of [basis_a, basis_b, basis_c],
    # β ≈ π / 3
    unit_cell_ref = UnitCell(lattice_constants, primitive)
    unit_cell_test = UnitCell(
        LatticeConstants(basis_a + basis_b, basis_b, basis_c), primitive
    )
    @test is_equivalent_unit_cell(unit_cell_test, unit_cell_ref)

    # equivalent unit cell defined by linear combination of [basis_a, basis_b, basis_c],
    # β ≈ 2π / 3
    unit_cell_ref = UnitCell(lattice_constants, primitive)
    unit_cell_test = UnitCell(
        LatticeConstants(basis_a - basis_b, basis_b, basis_c), primitive
    )
    @test is_equivalent_unit_cell(unit_cell_test, unit_cell_ref)
end

@testset "is_equivalent_unit_cell(::HexagonalLatticeConstants)" begin
    # --- Preparations

    a_ref = 2
    c_ref = 5
    lattice_constants_ref = HexagonalLatticeConstants(a_ref, c_ref)

    # --- Exercise functionality and check results

    # unit cells are equivalent
    lattice_constants_test = HexagonalLatticeConstants(a_ref + 1e-9, c_ref - 1e-9)
    @test is_equivalent_unit_cell(lattice_constants_test, lattice_constants_ref)

    # test unit cell and reference unit cell are unrelated
    lattice_constants_test = HexagonalLatticeConstants(2 * a_ref, c_ref)
    @test !is_equivalent_unit_cell(lattice_constants_test, lattice_constants_ref)

    # test unit cell and reference unit cell are for different lattice systems
    lattice_constants_test = CubicLatticeConstants(1)
    @test !is_equivalent_unit_cell(lattice_constants_test, lattice_constants_ref)
end
