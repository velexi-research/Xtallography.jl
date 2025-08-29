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
Tests for TetragonalUnitCell types and methods in unit_cell/tetragonal.jl (except for
conventional_cell() tests)
"""
# --- Imports

# Standard library
using Test
using LinearAlgebra: qr, I

# Xtallography package
using Xtallography

# Testing utilities
include("testing_utilities.jl")

# --- Tests

# ------ Types

@testset "TetragonalUnitCell(::NamedTuple,::UnitCellSymmetry) inner constructor" begin
    # --- Tests

    # Valid arguments
    lattice_constants_ = (a=1, c=5)
    symmetry_ = primitive_unit_cell_symmetry
    unit_cell = TetragonalUnitCell(lattice_constants_, symmetry_)

    @test lattice_constants(unit_cell) == lattice_constants_
    @test symmetry(unit_cell) == symmetry_

    # Invalid arguments
    lattice_constants_ = (a=1,)
    expected_message = (
        "Invalid lattice_constants argument passed to UnitCell{Tetragonal} " *
        "constructor. Expected keys: (:a, :c). " *
        "Provided keys: $(keys(lattice_constants_))."
    )
    @test_throws ArgumentError(expected_message) TetragonalUnitCell(
        lattice_constants_, symmetry_
    )
end

@testset "TetragonalUnitCell(::NamedTuple;::Centering,::Set) outer constructor" begin
    # --- Preparations

    lattice_constants_ = (a=1, c=5)

    # --- Tests

    # Default keyword arguments
    unit_cell = TetragonalUnitCell(lattice_constants_)

    @test lattice_constants(unit_cell) == lattice_constants_
    @test symmetry(unit_cell) == primitive_unit_cell_symmetry

    # Non-default centering keyword argument
    for centering_ in CENTERINGS
        unit_cell = TetragonalUnitCell(lattice_constants_; centering=centering_)

        @test lattice_constants(unit_cell) == lattice_constants_
        @test symmetry(unit_cell) == UnitCellSymmetry(centering_, Set{SymmetryElement}())
    end

    # Non-default symmetry_elements keyword argument
    symmetry_elements_ = [a_4_2]
    unit_cell = TetragonalUnitCell(lattice_constants_; symmetry_elements=symmetry_elements_)

    @test lattice_constants(unit_cell) == lattice_constants_
    @test symmetry(unit_cell) ==
        UnitCellSymmetry(primitive_centering, Set{SymmetryElement}(symmetry_elements_))
end

@testset "TetragonalUnitCell(::Real,::Real;::Centering,::Set) constructor: valid arguments" begin
    # --- Preparations

    a = 1
    c = 3

    # --- Exercise functionality and check results

    # Default keyword arguments
    unit_cell = TetragonalUnitCell(a, c)

    @test lattice_constants(unit_cell).a == a
    @test lattice_constants(unit_cell).c == c

    @test centering(unit_cell) === primitive_centering
    @test symmetry_elements(unit_cell) isa Set
    @test isempty(symmetry_elements(unit_cell))

    # Non-default centering keyword argument
    for centering_ in CENTERINGS
        unit_cell = TetragonalUnitCell(a, c; centering=centering_)

        @test lattice_constants(unit_cell).a == a
        @test lattice_constants(unit_cell).c == c
        @test centering(unit_cell) === centering_
        @test symmetry_elements(unit_cell) isa Set
        @test isempty(symmetry_elements(unit_cell))
    end

    # Non-default symmetry_elements keyword argument
    symmetry_elements_ = [a_4_2]
    unit_cell = TetragonalUnitCell(a, c; symmetry_elements=symmetry_elements_)

    @test lattice_constants(unit_cell).a == a
    @test lattice_constants(unit_cell).c == c
    @test centering(unit_cell) === primitive_centering
    @test symmetry_elements(unit_cell) isa Set
    @test symmetry_elements(unit_cell) == Set(symmetry_elements_)
end

@testset "TetragonalUnitCell(::Real,::Real;::Centering,::Set) constructor: invalid arguments" begin
    # --- Preparations

    # Valid arguments
    a = 1
    c = 3

    # --- Tests

    # ------ a

    # a = 0
    expected_message = "`a` must be positive"
    @test_throws DomainError(0, expected_message) TetragonalUnitCell(0, c)

    # a < 0
    expected_message = "`a` must be positive"
    @test_throws DomainError(-1, expected_message) TetragonalUnitCell(-1, c)

    # ------ c

    # c = 0
    expected_message = "`c` must be positive"
    @test_throws DomainError(0, expected_message) TetragonalUnitCell(a, 0)

    # c < 0
    expected_message = "`c` must be positive"
    @test_throws DomainError(-2.0, expected_message) TetragonalUnitCell(a, -2.0)

    # ------ symmetry elements

    # TODO
end

@testset "UnitCell(::UnitCell) copy constructor: tetragonal" begin
    # --- Preparations

    lattice_constants_ = (a=1, c=3)
    centering = primitive_centering
    symmetry_elements_ = [a_4_2]
    reference_unit_cell = TetragonalUnitCell(
        lattice_constants_,
        UnitCellSymmetry(centering; symmetry_elements=symmetry_elements_),
    )

    # --- Tests

    unit_cell = UnitCell(reference_unit_cell)

    @test unit_cell === reference_unit_cell
end

@testset "UnitCell(::Vector, ::Vector, ::Vector;::Bool,::Centering) constructor: tetragonal" begin
    # --- Preparations

    # Generate random rotation matrix
    rotations = [I, qr(rand(3)).Q, qr(rand(3)).Q]

    # --- Tests

    # ------ a < c

    a = 5
    c = 10

    basis_a = [a, 0, 0]
    basis_b = [0, a, 0]
    basis_c = [0, 0, c]

    expected_unit_cell = standardize(TetragonalUnitCell(a, c))

    # default centering keyword argument
    test_basis_rotations_and_permutations(
        rotations, expected_unit_cell, basis_a, basis_b, basis_c
    )

    # centering keyword argument provided
    for centering in (primitive_centering, body_centering)
        test_basis_rotations_and_permutations(
            rotations, expected_unit_cell, basis_a, basis_b, basis_c; centering=centering
        )
    end

    # ------ a > c

    a = 10
    c = 5

    basis_a = [a, 0, 0]
    basis_b = [0, a, 0]
    basis_c = [0, 0, c]

    expected_unit_cell = standardize(TetragonalUnitCell(a, c))

    # default centering keyword argument
    test_basis_rotations_and_permutations(
        rotations, expected_unit_cell, basis_a, basis_b, basis_c
    )

    # centering keyword argument provided
    for centering in (primitive_centering, body_centering)
        test_basis_rotations_and_permutations(
            rotations, expected_unit_cell, basis_a, basis_b, basis_c; centering=centering
        )
    end
end

# --- Methods

@testset "lattice_system(::TetragonalUnitCell)" begin
    lattice_constants_ = (a=1, c=2)
    unit_cell = TetragonalUnitCell(lattice_constants_)
    @test lattice_system(unit_cell) === tetragonal
end

@testset "lattice_constants(::TetragonalUnitCell)" begin
    lattice_constants_ = (a=1, c=3)
    unit_cell = TetragonalUnitCell(lattice_constants_)
    @test lattice_constants(unit_cell) == lattice_constants_
end

@testset "symmetry(::TetragonalUnitCell)" begin
    # --- Preparations

    lattice_constants_ = (a=1, c=3)
    symmetry_elements_ = Set([a_2_1, d_perp_110])

    # --- Tests

    for centering_ in CENTERINGS
        symmetry_ = UnitCellSymmetry(;
            centering=centering_, symmetry_elements=symmetry_elements_
        )
        unit_cell = TetragonalUnitCell(lattice_constants_, symmetry_)

        @test symmetry(unit_cell) == symmetry_
    end
end

@testset "centering(::TetragonalUnitCell)" begin
    # --- Preparations

    lattice_constants_ = (a=1, c=3)

    # --- Tests

    # Default centering
    unit_cell = TetragonalUnitCell(lattice_constants_)
    @test centering(unit_cell) == primitive_centering

    # Non-default centering
    for centering_ in CENTERINGS
        unit_cell = TetragonalUnitCell(lattice_constants_; centering=centering_)
        @test centering(unit_cell) === centering_
    end
end

@testset "symmetry_elements(::TetragonalUnitCell)" begin
    # --- Preparations

    lattice_constants_ = (a=1, c=3)

    # --- Tests

    # Default symmetry elements
    unit_cell = TetragonalUnitCell(lattice_constants_)
    @test symmetry_elements(unit_cell) isa Set{SymmetryElement}
    @test isempty(symmetry_elements(unit_cell))

    # Non-default symmetry elements
    symmetry_elements_ = Set([a_2_1, d_perp_110])
    unit_cell = TetragonalUnitCell(lattice_constants_; symmetry_elements=symmetry_elements_)
    @test symmetry_elements(unit_cell) == symmetry_elements_
end

@testset "is_bravais_lattice(::TetragonalUnitCell)" begin
    # --- Preparations

    a = 1
    c = 3

    # --- Tests

    # Valid Bravais lattices
    for centering_ in (primitive_centering, body_centering)
        unit_cell = TetragonalUnitCell(a, c; centering=centering_)
        @test is_bravais_lattice(unit_cell)
    end

    # Inalid Bravais lattices
    for centering_ in (face_centering, base_centering)
        unit_cell = TetragonalUnitCell(a, c; centering=centering_)
        @test !is_bravais_lattice(unit_cell)
    end
end

@testset "basis(::TetragonalUnitCell)" begin
    # --- Preparations

    a = 5
    c = 7
    lattice_constants = TetragonalUnitCell(a, c)

    # --- Exercise functionality

    basis_a, basis_b, basis_c = basis(lattice_constants)

    # --- Check results

    @test basis_a ≈ [a, 0, 0]
    @test basis_b ≈ [0, a, 0]
    @test basis_c ≈ [0, 0, c]
end

@testset "volume(::TetragonalUnitCell)" begin
    # --- Preparations

    a = 5
    c = 7
    unit_cell = TetragonalUnitCell(a, c)

    # --- Exercise functionality and check results

    @test volume(unit_cell) ≈ a^2 * c
end

@testset "surface_area(::TetragonalUnitCell)" begin
    # --- Preparations

    a = 5
    c = 7
    unit_cell = TetragonalUnitCell(a, c)

    # --- Exercise functionality and check results

    @test surface_area(unit_cell) ≈ 2 * a^2 + 4 * a * c
end

@testset "standardize(::TetragonalUnitCell)" begin
    # --- Tests

    # ------ Tetragonal lattices have no lattice constants conventions for primitive and
    #        body centerings

    a = 1.0
    c = 10.0
    reference_unit_cell = TetragonalUnitCell(a, c)

    # centering = primitive_centering
    standardized_unit_cell = standardize(
        TetragonalUnitCell(
            lattice_constants(reference_unit_cell); centering=primitive_centering
        ),
    )

    expected_unit_cell = reference_unit_cell
    @test standardized_unit_cell ≈ expected_unit_cell

    # centering = body_centering
    standardized_unit_cell = standardize(
        TetragonalUnitCell(lattice_constants(reference_unit_cell); centering=body_centering)
    )

    expected_unit_cell = TetragonalUnitCell(a, c; centering=body_centering)
    @test standardized_unit_cell ≈ expected_unit_cell

    # ------ Invalid centerings

    for centering_ in (face_centering, base_centering)
        unit_cell = TetragonalUnitCell(a, c; centering=centering_)

        expected_message =
            "Invalid Bravais lattice: " *
            "(lattice_system=Tetragonal, centering=$(nameof(typeof(centering_))))"

        @test_throws ArgumentError(expected_message) standardize(unit_cell)
    end
end

@testset "reduced_cell(::TetragonalUnitCell)" begin
    # --- Preparations

    a = 5
    c = 7
    reference_unit_cell = TetragonalUnitCell(a, c)
    basis_a, basis_b, basis_c = basis(reference_unit_cell)

    # --- Exercise functionality and check results

    # primitive unit cell
    unit_cell = TetragonalUnitCell(
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
    @test reduced_cell_ isa TetragonalUnitCell
    @test volume(reduced_cell_) ≈ volume(unit_cell)
    @test reduced_cell_ ≈ expected_reduced_cell

    # body-centered unit cell
    unit_cell = TetragonalUnitCell(
        lattice_constants(reference_unit_cell); centering=body_centering
    )

    expected_reduced_cell = reduced_cell(
        UnitCell(
            basis_a,
            basis_b,
            0.5 * (basis_a + basis_b + basis_c);
            identify_lattice_system=false,
            centering=primitive_centering,
        ),
    )

    reduced_cell_ = reduced_cell(unit_cell)
    @test reduced_cell_ isa TriclinicUnitCell
    @test volume(reduced_cell_) ≈ 0.5 * volume(unit_cell)
    @test reduced_cell_ ≈ expected_reduced_cell

    # face-centered unit cell
    unit_cell = TetragonalUnitCell(
        lattice_constants(reference_unit_cell); centering=face_centering
    )

    expected_reduced_cell = reduced_cell(
        TetragonalUnitCell(a / sqrt(2), c; centering=body_centering)
    )

    reduced_cell_ = reduced_cell(unit_cell)
    @test reduced_cell_ isa TriclinicUnitCell
    @test volume(reduced_cell_) ≈ 0.25 * volume(unit_cell)
    @test reduced_cell_ ≈ expected_reduced_cell
end

@testset "is_equivalent_unit_cell(::UnitCell, ::UnitCell)" begin
    # --- Preparations

    a = 2
    c = 5
    reference_unit_cell = TetragonalUnitCell(a, c)
    basis_a, basis_b, basis_c = basis(reference_unit_cell)

    # --- Exercise functionality and check results

    # equivalent tetragonal and triclinic unit cells
    tetragonal_unit_cell = reference_unit_cell
    triclinic_unit_cell = UnitCell(
        basis_a,
        basis_b,
        basis_c;
        identify_lattice_system=false,
        centering=primitive_centering,
    )
    @test is_equivalent_unit_cell(tetragonal_unit_cell, triclinic_unit_cell)

    # body-centered unit cell
    body_centering_unit_cell = TetragonalUnitCell(
        lattice_constants(reference_unit_cell); centering=body_centering
    )
    primitive_unit_cell = UnitCell(
        basis_a,
        basis_b,
        0.5 * (basis_a + basis_b + basis_c);
        identify_lattice_system=false,
        centering=primitive_centering,
    )
    @test is_equivalent_unit_cell(body_centering_unit_cell, primitive_unit_cell)

    # face-centered unit cell
    face_centering_unit_cell = TetragonalUnitCell(
        lattice_constants(reference_unit_cell); centering=face_centering
    )
    primitive_unit_cell = TetragonalUnitCell(a / sqrt(2), c; centering=body_centering)
    @test is_equivalent_unit_cell(face_centering_unit_cell, primitive_unit_cell)
end

@testset "is_equivalent_unit_cell(::TetragonalUnitCell)" begin
    # --- Preparations

    a_ref = 2
    c_ref = 5
    lattice_constants_ref = TetragonalUnitCell(a_ref, c_ref)

    # --- Exercise functionality and check results

    # unit cells are equivalent
    lattice_constants_test = TetragonalUnitCell(a_ref + 1e-9, c_ref - 1e-9)
    @test is_equivalent_unit_cell(lattice_constants_test, lattice_constants_ref)

    # test unit cell and reference unit cell are unrelated
    lattice_constants_test = TetragonalUnitCell(2 * a_ref, c_ref)
    @test !is_equivalent_unit_cell(lattice_constants_test, lattice_constants_ref)

    # test unit cell and reference unit cell are for different lattice systems
    lattice_constants_test = CubicUnitCell(1)
    @test !is_equivalent_unit_cell(lattice_constants_test, lattice_constants_ref)
end

@testset "isapprox(::TetragonalUnitCell)" begin
    # --- Preparations

    x = TetragonalUnitCell(1.0, 2.0)
    y = TetragonalUnitCell(1.5, 2.5)

    # --- Exercise functionality and check results

    # x ≈ (x + delta)
    @test x ≈ TetragonalUnitCell(1.0 + 1e-9, 2.0)
    @test x ≈ TetragonalUnitCell(1.0, 2.0 + 1e-9)

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

@testset "-(::TetragonalUnitCell)" begin
    # --- Tests

    lattice_constants_x = (a=1, c=5)
    x = TetragonalUnitCell(lattice_constants_x)

    lattice_constants_y = (a=2, c=2.3)
    y = TetragonalUnitCell(lattice_constants_y)

    @test x - y == TetragonalUnitCellDelta(
        lattice_constants_x.a - lattice_constants_y.a,
        lattice_constants_x.c - lattice_constants_y.c,
    )
end
