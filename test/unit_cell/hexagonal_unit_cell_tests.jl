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
Tests for HexagonalUnitCell types and methods in unit_cell/hexagonal.jl (except for
conventional_cell() tests)
"""
# --- Imports

# Standard library
using Test
using LinearAlgebra: cross, det, norm, qr, I

# Xtallography package
using Xtallography

# Testing utilities
include("testing_utilities.jl")

# --- Tests

# ------ Constructors

@testset "HexagonalUnitCell(::NamedTuple,::UnitCellSymmetry) inner constructor" begin
    # --- Tests

    # Valid arguments
    lattice_constants_ = (a=1, c=2)
    symmetry_ = primitive_unit_cell_symmetry
    unit_cell = HexagonalUnitCell(lattice_constants_, symmetry_)

    @test lattice_constants(unit_cell) == lattice_constants_
    @test symmetry(unit_cell) == symmetry_

    # Invalid arguments
    lattice_constants_ = (a=1,)
    expected_message = (
        "Invalid lattice_constants argument passed to UnitCell{Hexagonal} " *
        "constructor. Expected keys: (:a, :c). " *
        "Provided keys: $(keys(lattice_constants_))."
    )
    @test_throws ArgumentError(expected_message) HexagonalUnitCell(
        lattice_constants_, symmetry_
    )
end

@testset "HexagonalUnitCell(::NamedTuple;::Centering,::Set) outer constructor" begin
    # --- Preparations

    lattice_constants_ = (a=1, c=5)

    # --- Tests

    # Default keyword arguments
    unit_cell = HexagonalUnitCell(lattice_constants_)

    @test lattice_constants(unit_cell) == lattice_constants_
    @test symmetry(unit_cell) == primitive_unit_cell_symmetry

    # Non-default centering keyword argument
    for centering_ in CENTERINGS
        unit_cell = HexagonalUnitCell(lattice_constants_; centering=centering_)

        @test lattice_constants(unit_cell) == lattice_constants_
        @test symmetry(unit_cell) == UnitCellSymmetry(centering_, Set{SymmetryElement}())
    end

    # Non-default symmetry_elements keyword argument
    symmetry_elements_ = [a_4_2]
    unit_cell = HexagonalUnitCell(lattice_constants_; symmetry_elements=symmetry_elements_)

    @test lattice_constants(unit_cell) == lattice_constants_
    @test symmetry(unit_cell) ==
        UnitCellSymmetry(primitive_centering, Set{SymmetryElement}(symmetry_elements_))
end

@testset "HexagonalUnitCell(::Real,::Real;::Centering,::Set) outer constructor: valid arguments" begin
    # --- Preparations

    a = 1
    c = 3

    # --- Tests

    # Default keyword arguments
    unit_cell = HexagonalUnitCell(a, c)

    @test lattice_constants(unit_cell).a == a
    @test lattice_constants(unit_cell).c == c
    @test centering(unit_cell) === primitive_centering
    @test symmetry_elements(unit_cell) isa Set
    @test isempty(symmetry_elements(unit_cell))

    # Non-default centering keyword argument
    for centering_ in CENTERINGS
        unit_cell = HexagonalUnitCell(a, c; centering=centering_)

        @test lattice_constants(unit_cell).a == a
        @test lattice_constants(unit_cell).c == c
        @test centering(unit_cell) === centering_
        @test symmetry_elements(unit_cell) isa Set
        @test isempty(symmetry_elements(unit_cell))
    end

    # Non-default symmetry_elements keyword argument
    symmetry_elements_ = [a_4_2]
    unit_cell = HexagonalUnitCell(a, c; symmetry_elements=symmetry_elements_)

    @test lattice_constants(unit_cell).a == a
    @test lattice_constants(unit_cell).c == c
    @test centering(unit_cell) === primitive_centering
    @test symmetry_elements(unit_cell) isa Set
    @test Set(symmetry_elements(unit_cell)) == Set(symmetry_elements_)
end

@testset "HexagonalUnitCell(::Real,::Real;::Centering,::Set) outer constructor: invalid arguments" begin
    # --- Preparations

    # Valid arguments
    a = 1
    c = 3

    # --- Tests

    # ------ a

    # a = 0
    expected_message = "`a` must be positive"
    @test_throws DomainError(0, expected_message) HexagonalUnitCell(0, c)

    # a < 0
    expected_message = "`a` must be positive"
    @test_throws DomainError(-1, expected_message) HexagonalUnitCell(-1.0, c)

    # ------ c

    # c = 0
    expected_message = "`c` must be positive"
    @test_throws DomainError(0, expected_message) HexagonalUnitCell(a, 0)

    # c < 0
    expected_message = "`c` must be positive"
    @test_throws DomainError(-1.0, expected_message) HexagonalUnitCell(a, -1.0)

    # ------ symmetry elements

    # TODO
end

@testset "UnitCell(::UnitCell) copy constructor" begin
    # --- Preparations

    lattice_constants_ = (a=1, c=3)
    centering = primitive_centering
    symmetry_elements_ = [a_4_2]
    reference_unit_cell = HexagonalUnitCell(
        lattice_constants_,
        UnitCellSymmetry(centering; symmetry_elements=symmetry_elements_),
    )

    # --- Tests

    unit_cell = UnitCell(reference_unit_cell)

    @test unit_cell === reference_unit_cell
end

@testset "UnitCell(::Vector,::Vector,::Vector;::Bool,::Centering) outer constructor: hexagonal" begin
    # --- Preparations

    # Generate random rotation matrix
    rotations = [I, qr(rand(3)).Q, qr(rand(3)).Q]

    # --- Tests

    # ------ a < c, β = 2π / 3

    a = 3
    c = 5

    basis_a = [a, 0, 0]
    basis_b = [-0.5 * a, 0.5 * a * sqrt(3), 0]
    basis_c = [0, 0, c]

    expected_unit_cell = standardize(HexagonalUnitCell(a, c))

    # default centering keyword argument
    test_basis_rotations_and_permutations(
        rotations, expected_unit_cell, basis_a, basis_b, basis_c
    )

    # centering keyword argument provided
    for centering in (primitive_centering,)
        test_basis_rotations_and_permutations(
            rotations, expected_unit_cell, basis_a, basis_b, basis_c; centering=centering
        )
    end

    # ------ a < c, β = π / 3

    a = 3
    c = 5

    basis_a = [a, 0, 0]
    basis_b = [0.5 * a, 0.5 * a * sqrt(3), 0]
    basis_c = [0, 0, c]

    expected_unit_cell = standardize(HexagonalUnitCell(a, c))

    # default centering keyword argument
    test_basis_rotations_and_permutations(
        rotations, expected_unit_cell, basis_a, basis_b, basis_c
    )

    # centering keyword argument provided
    for centering in (primitive_centering,)
        test_basis_rotations_and_permutations(
            rotations, expected_unit_cell, basis_a, basis_b, basis_c; centering=centering
        )
    end

    # ------ a > c, β = 2π / 3

    a = 5
    c = 3

    basis_a = [a, 0, 0]
    basis_b = [-0.5 * a, 0.5 * a * sqrt(3), 0]
    basis_c = [0, 0, c]

    expected_unit_cell = standardize(HexagonalUnitCell(a, c))

    # default centering keyword argument
    test_basis_rotations_and_permutations(
        rotations, expected_unit_cell, basis_a, basis_b, basis_c
    )

    # centering keyword argument provided
    for centering in (primitive_centering,)
        test_basis_rotations_and_permutations(
            rotations, expected_unit_cell, basis_a, basis_b, basis_c; centering=centering
        )
    end

    # ------ a > c, β = π / 3

    a = 5
    c = 3

    basis_a = [a, 0, 0]
    basis_b = [0.5 * a, 0.5 * a * sqrt(3), 0]
    basis_c = [0, 0, c]

    expected_unit_cell = standardize(HexagonalUnitCell(a, c))

    # default centering keyword argument
    test_basis_rotations_and_permutations(
        rotations, expected_unit_cell, basis_a, basis_b, basis_c
    )

    # centering keyword argument provided
    for centering in (primitive_centering,)
        test_basis_rotations_and_permutations(
            rotations, expected_unit_cell, basis_a, basis_b, basis_c; centering=centering
        )
    end
end

# --- Methods

@testset "lattice_system(::HexagonalUnitCell)" begin
    lattice_constants_ = (a=1, c=2)
    unit_cell = HexagonalUnitCell(lattice_constants_)
    @test lattice_system(unit_cell) === hexagonal
end

@testset "lattice_constants(::HexagonalUnitCell)" begin
    lattice_constants_ = (a=1, c=2)
    unit_cell = HexagonalUnitCell(lattice_constants_)
    @test lattice_constants(unit_cell) == lattice_constants_
end

@testset "symmetry(::HexagonalUnitCell)" begin
    # --- Preparations

    lattice_constants_ = (a=1, c=2)
    symmetry_elements_ = Set([a_2_1, d_perp_110])

    # --- Tests

    for centering_ in CENTERINGS
        symmetry_ = UnitCellSymmetry(;
            centering=centering_, symmetry_elements=symmetry_elements_
        )
        unit_cell = HexagonalUnitCell(lattice_constants_, symmetry_)

        @test symmetry(unit_cell) == symmetry_
    end
end

@testset "centering(::HexagonalUnitCell)" begin
    # --- Preparations

    lattice_constants_ = (a=1, c=2)

    # --- Tests

    # Default centering
    unit_cell = HexagonalUnitCell(lattice_constants_)
    @test centering(unit_cell) === primitive_centering

    # Non-default centering
    for centering_ in CENTERINGS
        unit_cell = HexagonalUnitCell(lattice_constants_; centering=centering_)
        @test centering(unit_cell) == centering_
    end
end

@testset "symmetry_elements(::HexagonalUnitCell)" begin
    # --- Preparations

    lattice_constants_ = (a=1, c=2)

    # --- Tests

    # Default symmetry elements
    unit_cell = HexagonalUnitCell(lattice_constants_)
    @test symmetry_elements(unit_cell) isa Set{SymmetryElement}
    @test isempty(symmetry_elements(unit_cell))

    # Non-default symmetry elements
    symmetry_elements_ = Set([a_2_1, d_perp_110])
    unit_cell = HexagonalUnitCell(lattice_constants_; symmetry_elements=symmetry_elements_)
    @test symmetry_elements(unit_cell) == symmetry_elements_
end

@testset "is_bravais_lattice(::HexagonalUnitCell)" begin
    # --- Preparations

    a = 1
    c = 3

    # --- Tests

    # Valid Bravais lattices
    for centering_ in (primitive_centering,)
        unit_cell = HexagonalUnitCell(a, c; centering=centering_)
        @test is_bravais_lattice(unit_cell)
    end

    # Invalid Bravais lattices
    for centering_ in (body_centering, face_centering, base_centering)
        unit_cell = HexagonalUnitCell(a, c; centering=centering_)
        @test !is_bravais_lattice(unit_cell)
    end
end

@testset "basis(::HexagonalUnitCell)" begin
    # --- Preparations

    a = 2
    c = 5
    unit_cell = HexagonalUnitCell(a, c)

    # --- Exercise functionality

    basis_a, basis_b, basis_c = basis(unit_cell)

    # --- Check results

    @test basis_a ≈ [a, 0, 0]
    @test basis_b ≈ [-0.5 * a, a * sqrt(3) / 2, 0]
    @test basis_c ≈ [0, 0, c]
end

@testset "volume(::HexagonalUnitCell)" begin
    # --- Preparations

    # Construct basis vectors for unit cell
    a = 6
    c = 8
    unit_cell = HexagonalUnitCell(a, c)

    basis_a, basis_b, basis_c = basis(unit_cell)

    # --- Exercise functionality and check results

    @test volume(unit_cell) ≈ abs(det(hcat(basis_a, basis_b, basis_c)))
end

@testset "surface_area(::HexagonalUnitCell)" begin
    # --- Preparations

    # Construct basis vectors for unit cell
    a = 6
    c = 8
    unit_cell = HexagonalUnitCell(a, c)

    basis_a, basis_b, basis_c = basis(unit_cell)

    # --- Exercise functionality and check results

    @test surface_area(unit_cell) ≈
        2 * norm(cross(basis_a, basis_b)) +
          2 * norm(cross(basis_b, basis_c)) +
          2 * norm(cross(basis_c, basis_a))
end

@testset "standardize(::HexagonalUnitCell)" begin
    # --- Tests

    # ------ Hexagonal lattices have no lattice constants conventions for primitive
    #        centering

    unit_cell = HexagonalUnitCell(1.0, 2.0)

    standardized_unit_cell = standardize(unit_cell)

    expected_unit_cell = unit_cell
    @test standardized_unit_cell ≈ expected_unit_cell

    # ------ Invalid centerings

    for centering in (body_centering, face_centering, base_centering)
        unit_cell = HexagonalUnitCell(1, 2; centering=centering)
        expected_message =
            "Invalid Bravais lattice: " *
            "(lattice_system=Hexagonal, centering=$(nameof(typeof(centering))))"
        @test_throws ArgumentError(expected_message) standardize(unit_cell)
    end
end

@testset "reduced_cell(::HexagonalUnitCell)" begin
    # --- Preparations

    a = 2
    c = 5
    reference_unit_cell = HexagonalUnitCell(a, c)
    basis_a, basis_b, basis_c = basis(reference_unit_cell)

    # --- Exercise functionality and check results

    # primitive unit cell defined by [basis_a, basis_b, basis_c]
    unit_cell = HexagonalUnitCell(
        lattice_constants(reference_unit_cell); centering=primitive_centering
    )

    expected_reduced_cell = reduced_cell(
        UnitCell(
            basis_a,
            basis_b,
            basis_c;
            identify_lattice_system=false,
            centering=primitive_centering,
        ),
    )

    reduced_cell_ = reduced_cell(unit_cell)
    @test reduced_cell_ isa HexagonalUnitCell
    @test volume(reduced_cell_) ≈ volume(unit_cell)
    @test reduced_cell_ ≈ expected_reduced_cell

    # primitive unit cell defined by linear combination of [basis_a, basis_b, basis_c],
    # β ≈ π / 3
    unit_cell = UnitCell(basis_a + basis_b, basis_b, basis_c; centering=primitive_centering)

    expected_reduced_cell = reduced_cell(reference_unit_cell)

    reduced_cell_ = reduced_cell(unit_cell)
    @test reduced_cell_ isa HexagonalUnitCell
    @test volume(reduced_cell_) ≈ volume(unit_cell)
    @test reduced_cell_ ≈ expected_reduced_cell

    # primitive unit cell defined by linear combination of [basis_a, basis_b, basis_c],
    # β ≈ 2π / 3
    unit_cell = UnitCell(basis_a - basis_b, basis_b, basis_c; centering=primitive_centering)

    expected_reduced_cell = reduced_cell(reference_unit_cell)

    reduced_cell_ = reduced_cell(unit_cell)
    @test reduced_cell_ isa HexagonalUnitCell
    @test volume(reduced_cell_) ≈ volume(unit_cell)
    @test reduced_cell_ ≈ expected_reduced_cell
end

@testset "is_equivalent_unit_cell(::UnitCell): hexagonal" begin
    # --- Preparations

    a = 2
    c = 5
    reference_unit_cell = HexagonalUnitCell(a, c)
    basis_a, basis_b, basis_c = basis(reference_unit_cell)

    # --- Exercise functionality and check results

    # equivalent hexagonal and triclinic unit cells
    hexagonal_unit_cell = reference_unit_cell
    triclinic_unit_cell = UnitCell(
        basis_a,
        basis_b,
        basis_c;
        identify_lattice_system=false,
        centering=primitive_centering,
    )
    @test is_equivalent_unit_cell(hexagonal_unit_cell, triclinic_unit_cell)

    # equivalent unit cell defined by linear combination of [basis_a, basis_b, basis_c],
    # β ≈ π / 3
    unit_cell_ref = reference_unit_cell
    unit_cell_test = UnitCell(
        basis_a + basis_b, basis_b, basis_c; centering=primitive_centering
    )
    @test is_equivalent_unit_cell(unit_cell_test, unit_cell_ref)

    # equivalent unit cell defined by linear combination of [basis_a, basis_b, basis_c],
    # β ≈ 2π / 3
    unit_cell_ref = reference_unit_cell
    unit_cell_test = UnitCell(
        basis_a - basis_b, basis_b, basis_c; centering=primitive_centering
    )
    @test is_equivalent_unit_cell(unit_cell_test, unit_cell_ref)
end

@testset "is_equivalent_unit_cell(::HexagonalUnitCell)" begin
    # --- Preparations

    a_ref = 2
    c_ref = 5
    lattice_constants_ref = HexagonalUnitCell(a_ref, c_ref)

    # --- Exercise functionality and check results

    # unit cells are equivalent
    lattice_constants_test = HexagonalUnitCell(a_ref + 1e-9, c_ref - 1e-9)
    @test is_equivalent_unit_cell(lattice_constants_test, lattice_constants_ref)

    # test unit cell and reference unit cell are unrelated
    lattice_constants_test = HexagonalUnitCell(2 * a_ref, c_ref)
    @test !is_equivalent_unit_cell(lattice_constants_test, lattice_constants_ref)

    # test unit cell and reference unit cell are for different lattice systems
    lattice_constants_test = CubicUnitCell(1)
    @test !is_equivalent_unit_cell(lattice_constants_test, lattice_constants_ref)
end

@testset ":(==)(::HexagonalUnitCell)" begin
    # --- Exercise functionality and check results

    # unit cells are equal
    @test HexagonalUnitCell(1, 0.5) == HexagonalUnitCell(1, 0.5)

    # lattice constants are not equal
    @test HexagonalUnitCell(1, 0.5) != HexagonalUnitCell(2, 0.5)

    # centerings are not equal
    @test (
        HexagonalUnitCell(1, 0.5; centering=primitive_centering) !=
        HexagonalUnitCell(1, 0.5; centering=body_centering)
    )

    # symmetry elements are not equal
    @test (
        HexagonalUnitCell(1, 0.5; symmetry_elements=Set()) !=
        HexagonalUnitCell(1, 0.5; symmetry_elements=Set([a_2_1, d_perp_110]))
    )
end

@testset "isapprox(::HexagonalUnitCell)" begin
    # --- Preparations

    x = HexagonalUnitCell(1.0, 2.0)
    y = HexagonalUnitCell(1.5, 2.5)

    # --- Exercise functionality and check results

    # x ≈ (x + delta)
    @test x ≈ HexagonalUnitCell(1.0 + 1e-9, 2.0)
    @test x ≈ HexagonalUnitCell(1.0, 2.0 + 1e-9)

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

@testset "-(::HexagonalUnitCell)" begin
    # --- Tests

    lattice_constants_x = (a=1, c=3)
    x = HexagonalUnitCell(lattice_constants_x)

    lattice_constants_y = (a=2, c=10)
    y = HexagonalUnitCell(lattice_constants_y)

    @test x - y == HexagonalUnitCellDelta(
        lattice_constants_x.a - lattice_constants_y.a,
        lattice_constants_x.c - lattice_constants_y.c,
    )
end
