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
Tests for RhombohedralUnitCell types and methods in unit_cell/rhombohedral.jl (except for
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

# ------ Types

@testset "RhombohedralUnitCell(::NamedTuple,::UnitCellSymmetry) inner constructor" begin
    # --- Tests

    # Valid arguments
    lattice_constants_ = (a=1, α=3π / 5)
    symmetry_ = primitive_unit_cell_symmetry
    unit_cell = RhombohedralUnitCell(lattice_constants_, symmetry_)

    @test lattice_constants(unit_cell) == lattice_constants_
    @test symmetry(unit_cell) == symmetry_

    # Invalid arguments
    lattice_constants_ = (a=1,)
    expected_message = (
        "Invalid lattice_constants argument passed to UnitCell{Rhombohedral} " *
        "constructor. Expected keys: (:a, :α). " *
        "Provided keys: $(keys(lattice_constants_))."
    )
    @test_throws ArgumentError(expected_message) RhombohedralUnitCell(
        lattice_constants_, symmetry_
    )
end

@testset "RhombohedralUnitCell(::NamedTuple;::Centering,::Set) outer constructor" begin
    # --- Preparations

    lattice_constants_ = (a=1, α=3π / 5)

    # --- Tests

    # Default keyword arguments
    unit_cell = RhombohedralUnitCell(lattice_constants_)

    @test lattice_constants(unit_cell) == lattice_constants_
    @test symmetry(unit_cell) == primitive_unit_cell_symmetry

    # Non-default centering keyword argument
    for centering_ in CENTERINGS
        unit_cell = RhombohedralUnitCell(lattice_constants_; centering=centering_)

        @test lattice_constants(unit_cell) == lattice_constants_
        @test symmetry(unit_cell) == UnitCellSymmetry(centering_, Set{SymmetryElement}())
    end

    # Non-default symmetry_elements keyword argument
    symmetry_elements_ = [a_4_2]
    unit_cell = RhombohedralUnitCell(
        lattice_constants_; symmetry_elements=symmetry_elements_
    )

    @test lattice_constants(unit_cell) == lattice_constants_
    @test symmetry(unit_cell) ==
        UnitCellSymmetry(primitive_centering, Set{SymmetryElement}(symmetry_elements_))
end

@testset "RhombohedralUnitCell(::Real,::Real;::Centering,::Set) constructor: valid arguments" begin
    # --- Preparations

    a = 1
    α = π / 4

    # --- Exercise functionality and check results

    # Default keyword arguments
    unit_cell = RhombohedralUnitCell(a, α)

    @test lattice_constants(unit_cell).a == a
    @test lattice_constants(unit_cell).α == α

    @test centering(unit_cell) === primitive_centering
    @test symmetry_elements(unit_cell) isa Set
    @test isempty(symmetry_elements(unit_cell))

    # Non-default centering keyword argument
    for centering_ in CENTERINGS
        unit_cell = RhombohedralUnitCell(a, α; centering=centering_)

        @test lattice_constants(unit_cell).a == a
        @test lattice_constants(unit_cell).α == α
        @test centering(unit_cell) === centering_
        @test symmetry_elements(unit_cell) isa Set
        @test isempty(symmetry_elements(unit_cell))
    end

    # Non-default symmetry_elements keyword argument
    symmetry_elements_ = [a_4_2]
    unit_cell = RhombohedralUnitCell(a, α; symmetry_elements=symmetry_elements_)

    @test lattice_constants(unit_cell).a == a
    @test lattice_constants(unit_cell).α == α
    @test centering(unit_cell) === primitive_centering
    @test symmetry_elements(unit_cell) isa Set
    @test symmetry_elements(unit_cell) == Set(symmetry_elements_)
end

@testset "RhombohedralUnitCell(::Real,::Real;::Centering,::Set) constructor: invalid arguments" begin
    # --- Preparations

    # Valid arguments
    a = 1
    α = π / 4

    # --- Tests

    # ------ a

    # a = 0
    expected_message = "`a` must be positive"
    @test_throws DomainError(0, expected_message) RhombohedralUnitCell(0, α)

    # a < 0
    expected_message = "`a` must be positive"
    @test_throws DomainError(-2.0, expected_message) RhombohedralUnitCell(-2.0, α)

    # ------ α

    # α = 0
    expected_message = "`α` must satisfy 0 < α < 2π / 3"
    @test_throws DomainError(0, expected_message) RhombohedralUnitCell(a, 0)

    # α < 0
    expected_message = "`α` must satisfy 0 < α < 2π / 3"
    @test_throws DomainError(-1.0, expected_message) RhombohedralUnitCell(a, -1)

    # α = 2π / 3
    expected_message = "`α` must satisfy 0 < α < 2π / 3"
    @test_throws DomainError(2π / 3, expected_message) RhombohedralUnitCell(a, 2π / 3)

    # α > 2π / 3
    expected_message = "`α` must satisfy 0 < α < 2π / 3"
    @test_throws DomainError(3π / 2, expected_message) RhombohedralUnitCell(a, 3π / 2)

    # ------ symmetry elements

    # TODO
end

@testset "UnitCell(::UnitCell) copy constructor: orthorhombic" begin
    # --- Preparations

    lattice_constants_ = (a=1, α=3π / 5)
    centering = primitive_centering
    symmetry_elements_ = [a_4_2]
    reference_unit_cell = RhombohedralUnitCell(
        lattice_constants_,
        UnitCellSymmetry(centering; symmetry_elements=symmetry_elements_),
    )

    # --- Tests

    unit_cell = UnitCell(reference_unit_cell)

    @test unit_cell === reference_unit_cell
end

@testset "UnitCell(::Vector, ::Vector, ::Vector;::Bool,::Centering) constructor: rhombohedral" begin
    # --- Preparations

    # Generate random rotation matrix
    rotations = [I, qr(rand(3)).Q, qr(rand(3)).Q]

    # --- Tests

    a = 5
    α = 2π / 5

    basis_a = [a, 0, 0]
    basis_b = [a * cos(α), a * sin(α), 0]
    basis_c = [
        a * cos(α),
        a / sin(α) * (cos(α) - cos(α)^2),
        a / sin(α) * sqrt(1 - 3 * cos(α)^2 + 2 * cos(α)^3),
    ]

    expected_unit_cell = standardize(RhombohedralUnitCell(a, α))

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

@testset "lattice_system(::RhombohedralUnitCell)" begin
    lattice_constants_ = (a=1, α=π / 5)
    unit_cell = RhombohedralUnitCell(lattice_constants_)
    @test lattice_system(unit_cell) === rhombohedral
end

@testset "lattice_constants(::RhombohedralUnitCell)" begin
    lattice_constants_ = (a=1, α=π / 5)
    unit_cell = RhombohedralUnitCell(lattice_constants_)
    @test lattice_constants(unit_cell) == lattice_constants_
end

@testset "symmetry(::RhombohedralUnitCell)" begin
    # --- Preparations

    lattice_constants_ = (a=1, α=π / 5)
    symmetry_elements_ = Set([a_2_1, d_perp_110])

    # --- Tests

    for centering_ in CENTERINGS
        symmetry_ = UnitCellSymmetry(;
            centering=centering_, symmetry_elements=symmetry_elements_
        )
        unit_cell = RhombohedralUnitCell(lattice_constants_, symmetry_)

        @test symmetry(unit_cell) == symmetry_
    end
end

@testset "centering(::RhombohedralUnitCell)" begin
    # --- Preparations

    lattice_constants_ = (a=1, α=π / 5)

    # --- Tests

    # Default centering
    unit_cell = RhombohedralUnitCell(lattice_constants_)
    @test centering(unit_cell) == primitive_centering

    # Non-default centering
    for centering_ in CENTERINGS
        unit_cell = RhombohedralUnitCell(lattice_constants_; centering=centering_)
        @test centering(unit_cell) === centering_
    end
end

@testset "symmetry_elements(::RhombohedralUnitCell)" begin
    # --- Preparations

    lattice_constants_ = (a=1, α=π / 5)

    # --- Tests

    # Default symmetry elements
    unit_cell = RhombohedralUnitCell(lattice_constants_)
    @test symmetry_elements(unit_cell) isa Set{SymmetryElement}
    @test isempty(symmetry_elements(unit_cell))

    # Non-default symmetry elements
    symmetry_elements_ = Set([a_2_1, d_perp_110])
    unit_cell = RhombohedralUnitCell(
        lattice_constants_; symmetry_elements=symmetry_elements_
    )
    @test symmetry_elements(unit_cell) == symmetry_elements_
end

@testset "is_bravais_lattice(::RhombohedralUnitCell)" begin
    # --- Preparations

    lattice_constants_ = (a=1, α=π / 5)

    # --- Tests

    # Valid Bravais lattices
    for centering_ in (primitive_centering,)
        unit_cell = RhombohedralUnitCell(lattice_constants_; centering=centering_)
        @test is_bravais_lattice(unit_cell)
    end

    for centering_ in (body_centering, face_centering, base_centering)
        unit_cell = RhombohedralUnitCell(lattice_constants_; centering=centering_)
        @test !is_bravais_lattice(unit_cell)
    end
end

@testset "basis(::RhombohedralUnitCell)" begin
    # --- Preparations

    a = 2
    α = 3π / 5
    unit_cell = RhombohedralUnitCell(a, α)

    # --- Exercise functionality

    basis_a, basis_b, basis_c = basis(unit_cell)

    # --- Check results

    @test basis_a ≈ [a, 0, 0]
    @test basis_b ≈ [a * cos(α), a * sin(α), 0]
    @test basis_c ≈ [
        a * cos(α),
        a / sin(α) * (cos(α) - cos(α)^2),
        a / sin(α) * sqrt(1 - 3 * cos(α)^2 + 2 * cos(α)^3),
    ]
end

@testset "volume(::RhombohedralUnitCell)" begin
    # --- Preparations

    # Construct basis vectors for unit cell
    a = 6
    α = 3π / 5
    unit_cell = RhombohedralUnitCell(a, α)

    basis_a, basis_b, basis_c = basis(unit_cell)

    # --- Exercise functionality and check results

    @test volume(unit_cell) ≈ abs(det(hcat(basis_a, basis_b, basis_c)))
end

@testset "surface_area(::RhombohedralUnitCell)" begin
    # --- Preparations

    # Construct basis vectors for unit cell
    a = 6
    α = 3π / 5
    unit_cell = RhombohedralUnitCell(a, α)

    basis_a, basis_b, basis_c = basis(unit_cell)

    # --- Exercise functionality and check results

    @test surface_area(unit_cell) ≈
        2 * norm(cross(basis_a, basis_b)) +
          2 * norm(cross(basis_b, basis_c)) +
          2 * norm(cross(basis_c, basis_a))
end

@testset "standardize(::RhombohedralUnitCell)" begin
    # --- Tests

    # ------ Rhombohedral lattices have no lattice constants conventions for primitive
    #        centering

    unit_cell = RhombohedralUnitCell(1.0, 2π / 5)

    standardized_unit_cell = standardize(unit_cell)

    expected_unit_cell = unit_cell
    @test standardized_unit_cell ≈ expected_unit_cell

    # ------ Invalid centering

    for centering in (body_centering, face_centering, base_centering)
        unit_cell = RhombohedralUnitCell(1.0, 2π / 5; centering=centering)
        expected_message =
            "Invalid Bravais lattice: " *
            "(lattice_system=Rhombohedral, centering=$(nameof(typeof(centering))))"

        @test_throws ArgumentError(expected_message) standardize(unit_cell)
    end
end

@testset "reduced_cell(::RhombohedralUnitCell)" begin
    # --- Preparations

    a = 2
    α = 2π / 5
    reference_unit_cell = RhombohedralUnitCell(a, α)
    basis_a, basis_b, basis_c = basis(reference_unit_cell)

    # --- Exercise functionality and check results

    # primitive unit cell
    unit_cell = RhombohedralUnitCell(
        lattice_constants(reference_unit_cell); centering=primitive_centering
    )

    expected_reduced_cell = reduced_cell(
        UnitCell(basis_a, basis_b, basis_c; centering=primitive_centering)
    )

    reduced_cell_ = reduced_cell(unit_cell)
    @test reduced_cell_ isa RhombohedralUnitCell
    @test volume(reduced_cell_) ≈ volume(unit_cell)
    @test reduced_cell_ ≈ expected_reduced_cell
end

@testset "is_equivalent_unit_cell(::UnitCell): rhombohedral" begin
    # --- Preparations

    a = 2
    α = 3π / 5
    reference_unit_cell = RhombohedralUnitCell(a, α)
    basis_a, basis_b, basis_c = basis(reference_unit_cell)

    # --- Exercise functionality and check results

    # equivalent rhombohedral and triclinic unit cells
    rhombohedral_unit_cell = reference_unit_cell
    triclinic_unit_cell = UnitCell(
        basis_a,
        basis_b,
        basis_c;
        identify_lattice_system=false,
        centering=primitive_centering,
    )
    @test is_equivalent_unit_cell(rhombohedral_unit_cell, triclinic_unit_cell)
end

@testset "is_equivalent_unit_cell(::RhombohedralUnitCell)" begin
    # --- Preparations

    a_ref = 2
    α_ref = π / 4
    lattice_constants_ref = RhombohedralUnitCell(a_ref, α_ref)

    # --- Exercise functionality and check results

    # unit cells are equivalent
    lattice_constants_test = RhombohedralUnitCell(a_ref + 1e-9, α_ref - 1e-9)
    @test is_equivalent_unit_cell(lattice_constants_test, lattice_constants_ref)

    # test unit cell and reference unit cell are unrelated
    lattice_constants_test = RhombohedralUnitCell(2 * a_ref, α_ref)
    @test !is_equivalent_unit_cell(lattice_constants_test, lattice_constants_ref)

    # test unit cell and reference unit cell are for different lattice systems
    lattice_constants_test = CubicUnitCell(1)
    @test !is_equivalent_unit_cell(lattice_constants_test, lattice_constants_ref)
end

@testset "isapprox(::RhombohedralUnitCell)" begin
    # --- Preparations

    x = RhombohedralUnitCell(1.0, π / 4)
    y = RhombohedralUnitCell(1.5, π / 5)

    # --- Exercise functionality and check results

    # x ≈ (x + delta)
    @test x ≈ RhombohedralUnitCell(1.0 + 1e-9, π / 4)
    @test x ≈ RhombohedralUnitCell(1.0, π / 4 + 1e-9)

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

@testset "-(::RhombohedralUnitCell)" begin
    # --- Tests

    lattice_constants_x = (a=1, α=π / 6)
    x = RhombohedralUnitCell(lattice_constants_x)

    lattice_constants_y = (a=2, α=π / 5)
    y = RhombohedralUnitCell(lattice_constants_y)

    @test x - y == RhombohedralUnitCellDelta(
        lattice_constants_x.a - lattice_constants_y.a,
        lattice_constants_x.α - lattice_constants_y.α,
    )
end
