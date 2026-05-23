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
Tests for CubicUnitCell types and methods in unit_cell/cubic.jl (except for
conventional_cell() tests)
"""
# --- Imports

# Standard library
using LinearAlgebra: I
using Test

# Xtallography package
using Xtallography

# Testing utilities
include("testing_utilities.jl")

# --- Tests

# ------ Constructors

@testset "CubicUnitCell(::NamedTuple,::UnitCellSymmetry) inner constructor" begin
    # --- Tests

    # Valid arguments
    lattice_constants_ = (a=1,)
    symmetry_ = primitive_unit_cell_symmetry
    unit_cell = CubicUnitCell(lattice_constants_, symmetry_)

    @test lattice_constants(unit_cell) == lattice_constants_
    @test symmetry(unit_cell) == symmetry_

    # Invalid arguments
    lattice_constants_ = (a=1, b=3)
    expected_message = (
        "Invalid lattice_constants argument passed to UnitCell{Cubic} " *
        "constructor. Expected keys: (:a,). " *
        "Provided keys: $(keys(lattice_constants_))."
    )
    @test_throws ArgumentError(expected_message) CubicUnitCell(
        lattice_constants_, symmetry_
    )
end

@testset "CubicUnitCell(::NamedTuple;::Centering,::Set) outer constructor" begin
    # --- Preparations

    lattice_constants_ = (a=1,)

    # --- Tests

    # Default keyword arguments
    unit_cell = CubicUnitCell(lattice_constants_)

    @test lattice_constants(unit_cell) == lattice_constants_
    @test symmetry(unit_cell) == primitive_unit_cell_symmetry

    # Non-default centering keyword argument
    for centering_ in CENTERINGS
        unit_cell = CubicUnitCell(lattice_constants_; centering=centering_)

        @test lattice_constants(unit_cell) == lattice_constants_
        @test symmetry(unit_cell) == UnitCellSymmetry(centering_, Set{SymmetryElement}())
    end

    # Non-default symmetry_elements keyword argument
    symmetry_elements_ = [a_4_2]
    unit_cell = CubicUnitCell(lattice_constants_; symmetry_elements=symmetry_elements_)

    @test lattice_constants(unit_cell) == lattice_constants_
    @test symmetry(unit_cell) ==
        UnitCellSymmetry(primitive_centering, Set{SymmetryElement}(symmetry_elements_))
end

@testset "CubicUnitCell(::Real;::Centering,::Set) outer constructor: valid arguments" begin
    # --- Preparations

    a = 1

    # --- Tests

    # Default keyword arguments
    unit_cell = CubicUnitCell(a)

    @test lattice_constants(unit_cell).a == a
    @test centering(unit_cell) === primitive_centering
    @test symmetry_elements(unit_cell) isa Set
    @test isempty(symmetry_elements(unit_cell))

    # Non-default centering keyword argument
    for centering_ in CENTERINGS
        unit_cell = CubicUnitCell(a; centering=centering_)

        @test lattice_constants(unit_cell).a == a
        @test centering(unit_cell) === centering_
        @test symmetry_elements(unit_cell) isa Set
        @test isempty(symmetry_elements(unit_cell))
    end

    # Non-default symmetry_elements keyword argument
    symmetry_elements_ = [a_4_2]
    unit_cell = CubicUnitCell(a; symmetry_elements=symmetry_elements_)

    @test lattice_constants(unit_cell).a == a
    @test centering(unit_cell) === primitive_centering
    @test symmetry_elements(unit_cell) isa Set
    @test Set(symmetry_elements(unit_cell)) == Set(symmetry_elements_)
end

@testset "CubicUnitCell(::Real;::Centering,::Set) outer constructor: invalid arguments" begin
    # --- Tests

    # ------ a

    # a = 0
    expected_message = "`a` must be positive"
    @test_throws DomainError(0, expected_message) CubicUnitCell(0)

    # a < 0
    expected_message = "`a` must be positive"
    @test_throws DomainError(-1, expected_message) CubicUnitCell(-1.0)

    # ------ symmetry elements

    # TODO
end

@testset "UnitCell(::UnitCell) copy constructor: cubic" begin
    # --- Preparations

    a = 1
    reference_unit_cell = CubicUnitCell(a)

    # --- Tests

    unit_cell = UnitCell(reference_unit_cell)

    @test unit_cell === reference_unit_cell
end

@testset "UnitCell(::Vector,::Vector,::Vector;::Bool,::Centering) outer constructor: cubic" begin
    # --- Preparations

    # Generate random rotation matrix
    rotations = [I, qr(rand(3)).Q, qr(rand(3)).Q]

    # Generate basis
    a = 5
    basis_a = [a, 0, 0]
    basis_b = [0, a, 0]
    basis_c = [0, 0, a]

    # --- Tests

    # default centering keyword argument
    expected_unit_cell = standardize(CubicUnitCell(a; centering=P_centering))

    test_basis_rotations_and_permutations(
        rotations, expected_unit_cell, basis_a, basis_b, basis_c
    )

    # centering keyword argument provided
    for centering in (primitive_centering, body_centering, face_centering)
        expected_unit_cell = standardize(CubicUnitCell(a; centering=centering))

        test_basis_rotations_and_permutations(
            rotations, expected_unit_cell, basis_a, basis_b, basis_c; centering=centering
        )
    end
end

# ------ Methods

@testset "lattice_system(::CubicUnitCell)" begin
    lattice_constants_ = (a=1,)
    unit_cell = CubicUnitCell(lattice_constants_)

    @test lattice_system(unit_cell) === cubic
end

@testset "lattice_constants(::CubicUnitCell)" begin
    lattice_constants_ = (a=1,)
    unit_cell = CubicUnitCell(lattice_constants_)

    @test lattice_constants(unit_cell) == lattice_constants_
end

@testset "symmetry(::CubicUnitCell)" begin
    # --- Preparations

    lattice_constants_ = (a=1,)
    symmetry_elements_ = Set([a_2_1, d_perp_110])

    # --- Tests

    for centering_ in CENTERINGS
        symmetry_ = UnitCellSymmetry(;
            centering=centering_, symmetry_elements=symmetry_elements_
        )
        unit_cell = CubicUnitCell(lattice_constants_, symmetry_)

        @test symmetry(unit_cell) == symmetry_
    end
end

@testset "centering(::CubicUnitCell)" begin
    # --- Preparations

    lattice_constants_ = (a=1,)

    # --- Tests

    # Default centering
    unit_cell = CubicUnitCell(lattice_constants_)
    @test centering(unit_cell) === primitive_centering

    # Non-default centering
    for centering_ in CENTERINGS
        unit_cell = CubicUnitCell(lattice_constants_; centering=centering_)
        @test centering(unit_cell) == centering_
    end
end

@testset "symmetry_elements(::CubicUnitCell)" begin
    # --- Preparations

    lattice_constants_ = (a=1,)

    # --- Tests

    # Default symmetry elements
    unit_cell = CubicUnitCell(lattice_constants_)
    @test symmetry_elements(unit_cell) isa Set{SymmetryElement}
    @test isempty(symmetry_elements(unit_cell))

    # Non-default symmetry elements
    symmetry_elements_ = Set([a_2_1, d_perp_110])
    unit_cell = CubicUnitCell(lattice_constants_; symmetry_elements=symmetry_elements_)
    @test symmetry_elements(unit_cell) == symmetry_elements_
end

@testset "is_bravais_lattice(::CubicUnitCell)" begin
    # --- Preparations

    a = 1

    # --- Tests

    # Valid Bravais lattices
    for centering_ in (primitive_centering, body_centering, face_centering)
        unit_cell = CubicUnitCell(a; centering=centering_)
        @test is_bravais_lattice(unit_cell)
    end

    # Invalid Bravais lattices
    unit_cell = CubicUnitCell(a; centering=base_centering)
    @test !is_bravais_lattice(unit_cell)
end

@testset "basis(::CubicUnitCell)" begin
    # --- Preparations

    a = 5
    unit_cell = CubicUnitCell(a)

    # --- Exercise functionality

    basis_a, basis_b, basis_c = basis(unit_cell)

    # --- Check results

    @test basis_a ≈ [a, 0, 0]
    @test basis_b ≈ [0, a, 0]
    @test basis_c ≈ [0, 0, a]
end

@testset "volume(::CubicUnitCell)" begin
    # --- Preparations

    a = 5
    unit_cell = CubicUnitCell(a)

    # --- Exercise functionality and check results

    @test volume(unit_cell) ≈ a^3
end

@testset "surface_area(::CubicUnitCell)" begin
    # --- Preparations

    a = 5
    unit_cell = CubicUnitCell(a)

    # --- Exercise functionality and check results

    @test surface_area(unit_cell) ≈ 6 * a^2
end

@testset "standardize(::CubicUnitCell)" begin
    # --- Tests

    # ------ Cubic lattices have no lattice constants conventions for primitive, body, and
    #        face centerings

    for centering in (primitive_centering, body_centering, face_centering)
        unit_cell = CubicUnitCell(1.0; centering=centering)

        @test standardize(unit_cell) === unit_cell
    end

    # ------ Invalid centering

    for centering in (base_centering,)
        unit_cell = CubicUnitCell(1.0; centering=centering)
        expected_message =
            "Invalid Bravais lattice: " *
            "(lattice_system=Cubic, centering=$(nameof(typeof(centering))))"

        @test_throws ArgumentError(expected_message) standardize(unit_cell)
    end
end

@testset "reduced_cell(::CubicUnitCell)" begin
    # --- Preparations

    a = 5
    reference_unit_cell = CubicUnitCell(a; centering=primitive_centering)
    basis_a, basis_b, basis_c = basis(reference_unit_cell)

    # --- Exercise functionality and check results

    # primitive unit cell
    unit_cell = CubicUnitCell(
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
    @test reduced_cell_ isa CubicUnitCell
    @test volume(reduced_cell_) ≈ volume(reference_unit_cell)
    @test reduced_cell_ ≈ expected_reduced_cell

    # body-centered unit cell
    unit_cell = CubicUnitCell(
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
    @test reduced_cell_ isa RhombohedralUnitCell
    @test volume(reduced_cell_) ≈ 0.5 * volume(reference_unit_cell)
    @test reduced_cell_ ≈ expected_reduced_cell

    # face-centered unit cell
    unit_cell = CubicUnitCell(
        lattice_constants(reference_unit_cell); centering=face_centering
    )

    expected_reduced_cell = reduced_cell(
        UnitCell(
            0.5 * (basis_a + basis_b),
            0.5 * (basis_b - basis_c),
            0.5 * (basis_b + basis_c);
            identify_lattice_system=false,
            centering=primitive_centering,
        ),
    )

    reduced_cell_ = reduced_cell(unit_cell)
    @test reduced_cell_ isa TriclinicUnitCell
    @test volume(reduced_cell_) ≈ 0.25 * volume(reference_unit_cell)
    @test reduced_cell_ ≈ expected_reduced_cell
end

@testset "is_equivalent(::UnitCell): cubic" begin
    # --- Preparations

    a = 5
    primitive_cubic_unit_cell = CubicUnitCell(a)
    basis_a, basis_b, basis_c = basis(primitive_cubic_unit_cell)

    # --- Tests

    # equivalent cubic and primitive triclinic unit cells
    cubic_unit_cell = primitive_cubic_unit_cell
    primitive_unit_cell = UnitCell(
        basis_a,
        basis_b,
        basis_c;
        identify_lattice_system=false,
        centering=primitive_centering,
    )

    @test is_equivalent(cubic_unit_cell, primitive_unit_cell)

    # equivalent body-centered unit cell and primitive rhombohedral unit cell
    cubic_unit_cell = CubicUnitCell(
        lattice_constants(primitive_cubic_unit_cell); centering=body_centering
    )
    primitive_unit_cell = UnitCell(
        basis_a, basis_b, 0.5 * (basis_a + basis_b + basis_c); centering=primitive_centering
    )

    @test is_equivalent(cubic_unit_cell, primitive_unit_cell)

    # equivalent face-centered unit cell and primitive triclinic unit cell
    cubic_unit_cell = CubicUnitCell(
        lattice_constants(primitive_cubic_unit_cell); centering=face_centering
    )
    primitive_unit_cell = UnitCell(
        0.5 * (basis_a + basis_b),
        0.5 * (basis_b - basis_c),
        0.5 * (basis_b + basis_c);
        centering=primitive_centering,
    )
    @test is_equivalent(cubic_unit_cell, primitive_unit_cell)
end

@testset "is_equivalent(::CubicUnitCell)" begin
    # --- Preparations

    reference_unit_cell = CubicUnitCell(2.0)

    # --- Exercise functionality and check results

    # unit cells are equivalent
    unit_cell = CubicUnitCell(2.0 + 1e-9)
    @test is_equivalent(unit_cell, reference_unit_cell)

    # test unit cell is a supercell of the reference unit cell
    unit_cell = CubicUnitCell(10.0)
    @test !is_equivalent(unit_cell, reference_unit_cell)

    # test unit cell and reference unit cell are unrelated
    unit_cell = CubicUnitCell(3)
    @test !is_equivalent(unit_cell, reference_unit_cell)

    # test unit cell and reference unit cell are for different lattice systems
    unit_cell = OrthorhombicUnitCell(1, 2, 3)
    @test !is_equivalent(unit_cell, reference_unit_cell)
end

@testset "is_supercell(::CubicUnitCell): valid arguments" begin
    # --- Preparations

    reference_unit_cell = CubicUnitCell(2.5)

    # --- Exercise functionality and check results

    # ------ default keyword arguments

    # test unit cell is a supercell of reference unit cell
    unit_cell = CubicUnitCell(5.0 + 1e-9)
    @test is_supercell(unit_cell, reference_unit_cell)

    # test unit cell is not a supercell of reference unit cell
    unit_cell = CubicUnitCell(6.0)
    @test !is_supercell(unit_cell, reference_unit_cell)

    # ------ keyword arguments tests

    # `tol` large enough for multiplier deviation to pass check
    unit_cell = CubicUnitCell(6.0)
    @test is_supercell(unit_cell, reference_unit_cell; tol=0.45)

    # `tol` large enough for multiplier deviation to pass check but where multiplier
    # is close to 1 so that the test unit cell is not a proper supercell
    unit_cell = CubicUnitCell(3.0)
    @test !is_supercell(unit_cell, reference_unit_cell; tol=0.2)
end

@testset "is_supercell(::CubicUnitCell): invalid arguments" begin
    # --- Preparations

    reference_unit_cell = CubicUnitCell(2.5)
    test_unit_cell = CubicUnitCell(6.0)

    # --- Exercise functionality and check results

    # ------ `tol`

    # tol = 0
    expected_message = "`tol` must be positive"
    @test_throws DomainError(0, expected_message) is_supercell(
        test_unit_cell, reference_unit_cell; tol=0
    )

    # tol < 0
    expected_message = "`tol` must be positive"
    @test_throws DomainError(-0.1, expected_message) is_supercell(
        test_unit_cell, reference_unit_cell; tol=-0.1
    )
end

@testset ":(==)(::CubicUnitCell)" begin
    # --- Exercise functionality and check results

    # unit cells are equal
    @test CubicUnitCell(1.0) == CubicUnitCell(1.0)

    # lattice constants are not equal
    @test CubicUnitCell(1.0) != CubicUnitCell(2.0)

    # centerings are not equal
    @test (
        CubicUnitCell(1.0; centering=primitive_centering) !=
        CubicUnitCell(1.0; centering=body_centering)
    )

    # symmetry elements are not equal
    @test (
        CubicUnitCell(1.0; symmetry_elements=Set()) !=
        CubicUnitCell(1.0; symmetry_elements=Set([a_2_1, d_perp_110]))
    )
end

@testset "isapprox(::CubicUnitCell)" begin
    # --- Preparations

    x = CubicUnitCell(1.0)
    y = CubicUnitCell(2.0)

    # --- Exercise functionality and check results

    # x ≈ (x + delta)
    @test x ≈ CubicUnitCell(1 + 1e-9)

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

@testset "-(::CubicUnitCell)" begin
    # --- Tests

    x = CubicUnitCell(1)
    y = CubicUnitCell(2)

    @test x - y == CubicUnitCellDelta(lattice_constants(x).a - lattice_constants(y).a)
end
