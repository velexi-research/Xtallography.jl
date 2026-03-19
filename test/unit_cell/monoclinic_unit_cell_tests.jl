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
Tests for MonoclinicUnitCell types and methods in unit_cell/monoclinic.jl (except for
conventional_cell() tests)
"""
# --- Imports

# Standard library
using Test
using LinearAlgebra: cross, det, dot, norm, qr, I

# Xtallography package
using Xtallography

# Testing utilities
include("testing_utilities.jl")

# --- Tests

# ------ Types

@testset "MonoclinicUnitCell(::NamedTuple,::UnitCellSymmetry) inner constructor" begin
    # --- Tests

    # Valid arguments
    lattice_constants_ = (a=1, b=3, c=5, β=2π / 7)
    symmetry_ = primitive_unit_cell_symmetry
    unit_cell = MonoclinicUnitCell(lattice_constants_, symmetry_)

    @test lattice_constants(unit_cell) == lattice_constants_
    @test symmetry(unit_cell) == symmetry_

    # Invalid arguments
    lattice_constants_ = (a=1,)
    expected_message = (
        "Invalid lattice_constants argument passed to UnitCell{Monoclinic} " *
        "constructor. Expected keys: (:a, :b, :c, :β). " *
        "Provided keys: $(keys(lattice_constants_))."
    )
    @test_throws ArgumentError(expected_message) MonoclinicUnitCell(
        lattice_constants_, symmetry_
    )
end

@testset "MonoclinicUnitCell(::NamedTuple;::Centering,::Set) outer constructor" begin
    # --- Preparations

    lattice_constants_ = (a=1, b=3, c=5, β=2π / 7)

    # --- Tests

    # Default keyword arguments
    unit_cell = MonoclinicUnitCell(lattice_constants_)

    @test lattice_constants(unit_cell) == lattice_constants_
    @test symmetry(unit_cell) == primitive_unit_cell_symmetry

    # Non-default centering keyword argument
    for centering_ in CENTERINGS
        unit_cell = MonoclinicUnitCell(lattice_constants_; centering=centering_)

        @test lattice_constants(unit_cell) == lattice_constants_
        @test symmetry(unit_cell) == UnitCellSymmetry(centering_, Set{SymmetryElement}())
    end

    # Non-default symmetry_elements keyword argument
    symmetry_elements_ = [a_4_2]
    unit_cell = MonoclinicUnitCell(lattice_constants_; symmetry_elements=symmetry_elements_)

    @test lattice_constants(unit_cell) == lattice_constants_
    @test symmetry(unit_cell) ==
        UnitCellSymmetry(primitive_centering, Set{SymmetryElement}(symmetry_elements_))
end

@testset "MonoclinicUnitCell(::Real,::Real,::Real,::Real;::Centering,::Set) outer constructor: valid arguments" begin
    # --- Preparations

    a = 1
    b = 2
    c = 3
    β = π / 4

    # --- Exercise functionality and check results

    # Default keyword arguments
    unit_cell = MonoclinicUnitCell(a, b, c, β)

    @test lattice_constants(unit_cell).a == a
    @test lattice_constants(unit_cell).b == b
    @test lattice_constants(unit_cell).c == c
    @test lattice_constants(unit_cell).β == β

    @test centering(unit_cell) === primitive_centering
    @test symmetry_elements(unit_cell) isa Set
    @test isempty(symmetry_elements(unit_cell))

    # Non-default centering keyword argument
    for centering_ in CENTERINGS
        unit_cell = MonoclinicUnitCell(a, b, c, β; centering=centering_)

        @test lattice_constants(unit_cell).a == a
        @test lattice_constants(unit_cell).b == b
        @test lattice_constants(unit_cell).c == c
        @test lattice_constants(unit_cell).β == β
        @test centering(unit_cell) === centering_
        @test symmetry_elements(unit_cell) isa Set
        @test isempty(symmetry_elements(unit_cell))
    end

    # Non-default symmetry_elements keyword argument
    symmetry_elements_ = [a_4_2]
    unit_cell = MonoclinicUnitCell(a, b, c, β; symmetry_elements=symmetry_elements_)

    @test lattice_constants(unit_cell).a == a
    @test lattice_constants(unit_cell).b == b
    @test lattice_constants(unit_cell).c == c
    @test lattice_constants(unit_cell).β == β
    @test centering(unit_cell) === primitive_centering
    @test symmetry_elements(unit_cell) isa Set
    @test symmetry_elements(unit_cell) == Set(symmetry_elements_)
end

@testset "MonoclinicUnitCell(::Real,::Real,::Real,::Real;::Centering,::Set) outer constructor: invalid arguments" begin
    # --- Preparations

    # Valid arguments
    a = 1
    b = 2
    c = 3
    β = π / 4

    # --- Tests

    # ------ a

    # a = 0
    expected_message = "`a` must be positive"
    @test_throws DomainError(0, expected_message) MonoclinicUnitCell(0, b, c, β)

    # a < 0
    expected_message = "`a` must be positive"
    @test_throws DomainError(-1.0, expected_message) MonoclinicUnitCell(-1.0, b, c, β)
    # ------ b

    # b = 0
    expected_message = "`b` must be positive"
    @test_throws DomainError(0, expected_message) MonoclinicUnitCell(a, 0, c, β)

    # b < 0
    expected_message = "`b` must be positive"
    @test_throws DomainError(-2.0, expected_message) MonoclinicUnitCell(a, -2.0, c, β)

    # ------ c

    # c = 0
    expected_message = "`c` must be positive"
    @test_throws DomainError(0, expected_message) MonoclinicUnitCell(a, b, 0, β)

    # c < 0
    expected_message = "`c` must be positive"
    @test_throws DomainError(-3.0, expected_message) MonoclinicUnitCell(a, b, -3.0, β)

    # ------ β

    # β = 0
    expected_message = "`β` must satisfy 0 < β < π"
    @test_throws DomainError(0, expected_message) MonoclinicUnitCell(a, b, c, 0)

    # β < 0
    expected_message = "`β` must satisfy 0 < β < π"
    @test_throws DomainError(-1, expected_message) MonoclinicUnitCell(a, b, c, -1)

    # β > π
    expected_message = "`β` must satisfy 0 < β < π"
    @test_throws DomainError(π + 1, expected_message) MonoclinicUnitCell(a, b, c, π + 1)

    # β = π
    expected_message = "`β` must satisfy 0 < β < π"
    @test_throws DomainError(π, expected_message) MonoclinicUnitCell(a, b, c, π)

    # ------ symmetry elements

    # TODO
end

@testset "UnitCell(::UnitCell) copy constructor: monoclinic" begin
    # --- Preparations

    lattice_constants_ = (a=1, b=2, c=3, β=2π / 5)
    centering = primitive_centering
    symmetry_elements_ = [a_4_2]
    reference_unit_cell = MonoclinicUnitCell(
        lattice_constants_,
        UnitCellSymmetry(centering; symmetry_elements=symmetry_elements_),
    )

    # --- Tests

    unit_cell = UnitCell(reference_unit_cell)

    @test unit_cell === reference_unit_cell
end

@testset "UnitCell(::Vector,::Vector,::Vector;::Bool,::Centering) outer constructor: monoclinic" begin
    # --- Preparations

    # Generate random rotation matrix
    rotations = [I, qr(rand(3)).Q, qr(rand(3)).Q]

    # --- Tests

    a = 5
    b = 8
    c = 10
    β = 3π / 5

    basis_a = [a, 0, 0]
    basis_b = [0, b, 0]
    basis_c = [c * cos(β), 0, c * sin(β)]

    # default centering keyword argument
    expected_unit_cell = standardize(MonoclinicUnitCell(a, b, c, β))
    test_basis_rotations_and_permutations(
        rotations, expected_unit_cell, basis_a, basis_b, basis_c
    )

    # centering = primitive
    expected_unit_cell = standardize(MonoclinicUnitCell(a, b, c, β))
    test_basis_rotations_and_permutations(
        rotations,
        expected_unit_cell,
        basis_a,
        basis_b,
        basis_c;
        centering=primitive_centering,
    )

    # centering = body-centering
    expected_unit_cell = standardize(
        MonoclinicUnitCell(a, b, c, β; centering=body_centering)
    )
    test_basis_rotations_and_permutations(
        rotations, expected_unit_cell, basis_a, basis_b, basis_c; centering=body_centering
    )

    # centering = base-centering
    expected_unit_cell = standardize(
        MonoclinicUnitCell(a, b, c, β; centering=base_centering)
    )
    test_basis_rotations_and_permutations(
        rotations, expected_unit_cell, basis_a, basis_b, basis_c; centering=base_centering
    )
end

# ------ Methods

@testset "lattice_system(::MonoclinicUnitCell)" begin
    lattice_constants_ = (a=1, b=2, c=3, β=2π / 5)
    unit_cell = MonoclinicUnitCell(lattice_constants_)
    @test lattice_system(unit_cell) === monoclinic
end

@testset "lattice_constants(::MonoclinicUnitCell)" begin
    lattice_constants_ = (a=1, b=2, c=3, β=π / 2)
    unit_cell = MonoclinicUnitCell(lattice_constants_)
    @test lattice_constants(unit_cell) == lattice_constants_
end

@testset "symmetry(::MonoclinicUnitCell)" begin
    # --- Preparations

    lattice_constants_ = (a=1, b=2, c=3, β=π / 2)
    symmetry_elements_ = Set([a_2_1, d_perp_110])

    # --- Tests

    for centering_ in CENTERINGS
        symmetry_ = UnitCellSymmetry(;
            centering=centering_, symmetry_elements=symmetry_elements_
        )
        unit_cell = MonoclinicUnitCell(lattice_constants_, symmetry_)

        @test symmetry(unit_cell) == symmetry_
    end
end

@testset "centering(::MonoclinicUnitCell)" begin
    # --- Preparations

    lattice_constants_ = (a=1, b=2, c=3, β=π / 2)

    # --- Tests

    # Default centering
    unit_cell = MonoclinicUnitCell(lattice_constants_)
    @test centering(unit_cell) == primitive_centering

    # Non-default centering
    for centering_ in CENTERINGS
        unit_cell = MonoclinicUnitCell(lattice_constants_; centering=centering_)
        @test centering(unit_cell) === centering_
    end
end

@testset "symmetry_elements(::MonoclinicUnitCell)" begin
    # --- Preparations

    lattice_constants_ = (a=1, b=2, c=3, β=π / 2)

    # --- Tests

    # Default symmetry elements
    unit_cell = MonoclinicUnitCell(lattice_constants_)
    @test symmetry_elements(unit_cell) isa Set{SymmetryElement}
    @test isempty(symmetry_elements(unit_cell))

    # Non-default symmetry elements
    symmetry_elements_ = Set([a_2_1, d_perp_110])
    unit_cell = MonoclinicUnitCell(lattice_constants_; symmetry_elements=symmetry_elements_)
    @test symmetry_elements(unit_cell) == symmetry_elements_
end

@testset "is_bravais_lattice(::MonoclinicUnitCell)" begin
    # --- Preparations

    a = 1
    b = 2
    c = 3
    β = π / 2

    # --- Tests

    # Valid Bravais lattices
    for centering_ in (primitive_centering, body_centering, base_centering)
        unit_cell = MonoclinicUnitCell(a, b, c, β; centering=centering_)
        @test is_bravais_lattice(unit_cell)
    end

    # Invalid Bravais lattices
    for centering_ in (face_centering,)
        unit_cell = MonoclinicUnitCell(a, b, c, β; centering=centering_)
        @test !is_bravais_lattice(unit_cell)
    end
end

@testset "basis(::MonoclinicUnitCell)" begin
    # --- Preparations

    a = 2
    b = 3
    c = 5
    β = 3π / 5
    unit_cell = MonoclinicUnitCell(a, b, c, β)

    # --- Exercise functionality

    basis_a, basis_b, basis_c = basis(unit_cell)

    # --- Check results

    @test basis_a ≈ [a, 0, 0]
    @test basis_b ≈ [0, b, 0]
    @test basis_c ≈ [c * cos(β), 0, c * sin(β)]
end

@testset "volume(::MonoclinicUnitCell)" begin
    # --- Preparations

    # Construct basis vectors for unit cell
    a = 6
    b = 3
    c = 10
    β = 3π / 5
    unit_cell = MonoclinicUnitCell(a, b, c, β)

    basis_a, basis_b, basis_c = basis(unit_cell)

    # --- Exercise functionality and check results

    @test volume(unit_cell) ≈ abs(det(hcat(basis_a, basis_b, basis_c)))
end

@testset "surface_area(::MonoclinicUnitCell)" begin
    # --- Preparations

    # Construct basis vectors for unit cell
    a = 6
    b = 3
    c = 10
    β = 3π / 5
    unit_cell = MonoclinicUnitCell(a, b, c, β)

    basis_a, basis_b, basis_c = basis(unit_cell)

    # --- Exercise functionality and check results

    @test surface_area(unit_cell) ≈
        2 * norm(cross(basis_a, basis_b)) +
          2 * norm(cross(basis_b, basis_c)) +
          2 * norm(cross(basis_c, basis_a))
end

@testset "standardize(::MonoclinicUnitCell)" begin
    # --- Tests

    # ------ lattice constants already in standard form

    a = 6
    b = 10
    c = 8
    β = 1.1 * π / 2
    unit_cell = MonoclinicUnitCell(a, b, c, β)

    # centering = primitive_centering
    standardized_unit_cell = standardize(unit_cell)

    expected_unit_cell = MonoclinicUnitCell(a, b, c, β)
    @test standardized_unit_cell ≈ expected_unit_cell

    # centering = body_centering
    standardized_unit_cell = standardize(
        MonoclinicUnitCell(lattice_constants(unit_cell); centering=body_centering)
    )

    expected_unit_cell = MonoclinicUnitCell(a, b, c, β; centering=body_centering)
    @test standardized_unit_cell ≈ expected_unit_cell

    # ------ β ∉ [π/2, π]

    a = 6
    b = 10
    c = 8
    β = π - 1.1 * π / 2
    unit_cell = MonoclinicUnitCell(a, b, c, β)

    # centering = primitive_centering
    standardized_unit_cell = standardize(unit_cell)

    expected_unit_cell = MonoclinicUnitCell(a, b, c, π - β)
    @test standardized_unit_cell ≈ expected_unit_cell

    # centering = body_centering
    standardized_unit_cell = standardize(
        MonoclinicUnitCell(lattice_constants(unit_cell); centering=body_centering)
    )

    expected_unit_cell = MonoclinicUnitCell(a, b, c, π - β; centering=body_centering)
    @test standardized_unit_cell ≈ expected_unit_cell

    # ------ a > c

    a = 8
    b = 10
    c = 6
    β = 1.1 * π / 2
    unit_cell = MonoclinicUnitCell(a, b, c, β)

    # centering = primitive_centering
    standardized_unit_cell = standardize(unit_cell)

    expected_unit_cell = MonoclinicUnitCell(c, b, a, β)
    @test standardized_unit_cell ≈ expected_unit_cell

    # centering = body_centering
    standardized_unit_cell = standardize(
        MonoclinicUnitCell(lattice_constants(unit_cell); centering=body_centering)
    )

    expected_unit_cell = MonoclinicUnitCell(c, b, a, β; centering=body_centering)
    @test standardized_unit_cell ≈ expected_unit_cell

    # ------ base-centered unit cell, a < c

    a = 6
    b = 10
    c = 8
    β = 1.1 * π / 2
    unit_cell = MonoclinicUnitCell(a, b, c, β)

    standardized_unit_cell = standardize(
        MonoclinicUnitCell(lattice_constants(unit_cell); centering=base_centering)
    )

    a_body = sqrt(a^2 + c^2 + 2 * a * c * cos(β))
    β_body = π - asin(sin(β) / a_body * a)
    expected_unit_cell = MonoclinicUnitCell(c, b, a_body, β_body; centering=body_centering)
    @test standardized_unit_cell ≈ expected_unit_cell

    # ------ base-centered unit cell, a > c

    a = 8
    b = 10
    c = 6
    β = 1.1 * π / 2
    unit_cell = MonoclinicUnitCell(a, b, c, β)

    standardized_unit_cell = standardize(
        MonoclinicUnitCell(lattice_constants(unit_cell); centering=base_centering)
    )

    a_body = sqrt(a^2 + c^2 + 2 * a * c * cos(β))
    β_body = π - asin(sin(β) / a_body * a)
    expected_unit_cell = MonoclinicUnitCell(c, b, a_body, β_body; centering=body_centering)
    @test standardized_unit_cell ≈ expected_unit_cell

    # ------ unit cell requires single reduction in plane normal to b-axis

    a_ref = 6
    b_ref = 10
    c_ref = 8
    β_ref = 2π / 3

    # centering = primitive_centering
    a = a_ref
    b = b_ref
    c = c_ref
    β = β_ref
    unit_cell = MonoclinicUnitCell(a, b, c, β)

    standardized_unit_cell = standardize(unit_cell)

    expected_c = sqrt(a^2 + c^2 + 2 * a * c * cos(β))
    expected_β = π - asin(sin(β) / expected_c * c)
    expected_unit_cell = MonoclinicUnitCell(a, b, expected_c, expected_β)
    @test standardized_unit_cell ≈ expected_unit_cell

    # centering = body_centering
    a = a_ref
    b = b_ref
    c = sqrt((2 * a_ref)^2 + c_ref^2 + 2 * (2 * a_ref) * c_ref * cos(β_ref))
    β = π - asin(sin(β_ref) / c * c_ref)
    unit_cell = MonoclinicUnitCell(a, b, c, β; centering=body_centering)

    standardized_unit_cell = standardize(unit_cell)

    expected_unit_cell = MonoclinicUnitCell(
        a_ref, b_ref, c_ref, β_ref; centering=body_centering
    )
    @test standardized_unit_cell ≈ expected_unit_cell

    # ------ unit cell requires multiple reductions in plane normal to b-axis

    a_ref = 6
    b_ref = 10
    c_ref = 8
    β_ref = 2π / 3

    # centering = primitive_centering
    a = a_ref
    b = b_ref
    c = sqrt((5 * a_ref)^2 + c_ref^2 + 2 * (5 * a_ref) * c_ref * cos(β_ref))
    β = asin(sin(β_ref) / c * c_ref)
    unit_cell = MonoclinicUnitCell(a, b, c, β)

    standardized_unit_cell = standardize(unit_cell)

    expected_c = sqrt(a_ref^2 + c_ref^2 + 2 * a_ref * c_ref * cos(β_ref))
    expected_β = π - asin(sin(β_ref) / expected_c * c_ref)
    expected_unit_cell = MonoclinicUnitCell(a, b, expected_c, expected_β)
    @test standardized_unit_cell ≈ expected_unit_cell

    # centering = body_centering
    a = a_ref
    b = b_ref
    c = sqrt((5 * a_ref)^2 + c_ref^2 + 2 * (5 * a_ref) * c_ref * cos(β_ref))
    β = asin(sin(β_ref) / c * c_ref)
    unit_cell = MonoclinicUnitCell(a, b, c, β; centering=body_centering)

    standardized_unit_cell = standardize(unit_cell)

    expected_unit_cell = MonoclinicUnitCell(
        a, b, expected_c, expected_β; centering=body_centering
    )
    @test standardized_unit_cell ≈ expected_unit_cell

    # ------ Invalid centering

    # centering = face_centering
    centering_ = face_centering
    unit_cell = MonoclinicUnitCell(lattice_constants(unit_cell); centering=centering_)
    expected_message =
        "Invalid Bravais lattice: " *
        "(lattice_system=Monoclinic, centering=$(nameof(typeof(centering_))))"
    @test_throws ArgumentError(expected_message) standardize(unit_cell)
end

@testset "convert_to_body_centering()" begin
    # --- Preparations

    # Construct basis vectors for body-centered unit cell
    a = 6
    b = 3
    c = 10
    β = 3π / 5
    base_centered_unit_cell = MonoclinicUnitCell(a, b, c, β; centering=base_centering)

    # --- Exercise functionality and check results

    # Check conversion to body-centering
    body_centered_unit_cell = convert_to_body_centering(base_centered_unit_cell)

    a_body = sqrt(a^2 + c^2 + 2 * a * c * cos(β))
    c_body = c
    β_body = π - acos((a_body^2 + c_body^2 - a^2) / 2 / a_body / c_body)
    expected_body_centered_unit_cell = MonoclinicUnitCell(
        a_body, b, c_body, β_body; centering=body_centering
    )
    @test body_centered_unit_cell ≈ expected_body_centered_unit_cell

    # Check conversion back to base-centering
    #
    # Note: when the body-centered unit cell is converted back to base-centering, a
    #       different unit cell is selected because the original base-centered unit cell
    #       had a larger value of a^2 + c^2.
    expected_base_centered_unit_cell = MonoclinicUnitCell(
        a, b, a_body, π - asin(sin(β) / a_body * c); centering=base_centering
    )
    @test convert_to_base_centering(body_centered_unit_cell) ≈
        expected_base_centered_unit_cell
end

@testset "convert_to_base_centering()" begin
    # --- Tests

    # ------ a < c, β > π/2

    # Construct basis vectors for body-centered unit cell
    a = 6
    b = 3
    c = 10
    β = 3π / 5
    body_centered_unit_cell = MonoclinicUnitCell(a, b, c, β; centering=body_centering)

    # Check conversion to base-centering
    base_centered_unit_cell = convert_to_base_centering(body_centered_unit_cell)

    a_base = sqrt(a^2 + c^2 + 2 * a * c * cos(β))
    c_base = a
    β_base = π - acos((a_base^2 + c_base^2 - c^2) / 2 / a_base / c_base)
    expected_base_centered_unit_cell = MonoclinicUnitCell(
        a_base, b, c_base, β_base; centering=base_centering
    )
    @test base_centered_unit_cell ≈ expected_base_centered_unit_cell

    # Check conversion back to body-centering via standardize()
    standardized_unit_cell = standardize(base_centered_unit_cell)

    @test standardized_unit_cell ≈ body_centered_unit_cell

    # ------ a > c, β < π/2

    # Construct basis vectors for body-centered unit cell
    a = 8
    b = 3
    c = 4
    β = 2π / 5
    body_centered_unit_cell = MonoclinicUnitCell(a, b, c, β; centering=body_centering)

    # Check conversion to base-centering
    base_centered_unit_cell = convert_to_base_centering(body_centered_unit_cell)

    a_base = sqrt(a^2 + c^2 - 2 * a * c * cos(β))
    c_base = c
    β_base = π - acos((a_base^2 + c_base^2 - a^2) / 2 / a_base / c_base)
    expected_base_centered_unit_cell = MonoclinicUnitCell(
        a_base, b, c_base, β_base; centering=base_centering
    )
    @test base_centered_unit_cell ≈ expected_base_centered_unit_cell

    # Check conversion back to body-centering via standardize()
    standardized_unit_cell = standardize(base_centered_unit_cell)
    standardized_body_centered_unit_cell = standardize(body_centered_unit_cell)

    @test standardized_unit_cell ≈ standardized_body_centered_unit_cell
end

@testset "reduced_cell(::MonoclinicUnitCell)" begin
    # --- Preparations

    a = 2
    b = 3
    c = 5
    β = 4π / 7
    reference_unit_cell = MonoclinicUnitCell(a, b, c, β)
    basis_a, basis_b, basis_c = basis(reference_unit_cell)

    # --- Exercise functionality and check results

    # primitive unit cell defined by [basis_a, basis_b, basis_c]
    unit_cell = MonoclinicUnitCell(
        lattice_constants(reference_unit_cell); centering=primitive_centering
    )

    expected_reduced_cell = reduced_cell(unit_cell)

    reduced_cell_ = reduced_cell(unit_cell)
    @test reduced_cell_ isa MonoclinicUnitCell
    @test volume(reduced_cell_) ≈ volume(reference_unit_cell)
    @test reduced_cell_ ≈ expected_reduced_cell

    # primitive unit cell defined by linear combination of [basis_a, basis_b, basis_c]
    unit_cell = UnitCell(basis_a + basis_c, basis_b, basis_c; centering=primitive_centering)

    expected_reduced_cell = reduced_cell(unit_cell)

    reduced_cell_ = reduced_cell(unit_cell)
    @test reduced_cell_ isa MonoclinicUnitCell
    @test volume(reduced_cell_) ≈ volume(reference_unit_cell)
    @test reduced_cell_ ≈ expected_reduced_cell

    # body-centered unit cell
    unit_cell = MonoclinicUnitCell(
        lattice_constants(reference_unit_cell); centering=body_centering
    )

    expected_reduced_cell = reduced_cell(
        UnitCell(
            basis_a,
            basis_b,
            0.5 * (basis_a + basis_b + basis_c);
            centering=primitive_centering,
        ),
    )

    reduced_cell_ = reduced_cell(unit_cell)
    @test reduced_cell_ isa TriclinicUnitCell
    @test volume(reduced_cell_) ≈ 0.5 * volume(reference_unit_cell)
    @test reduced_cell_ ≈ expected_reduced_cell

    # base-centered unit cell
    unit_cell = MonoclinicUnitCell(
        lattice_constants(reference_unit_cell); centering=base_centering
    )

    expected_reduced_cell = reduced_cell(
        UnitCell(basis_a, 0.5 * (basis_a + basis_b), basis_c; centering=primitive_centering)
    )

    reduced_cell_ = reduced_cell(unit_cell)
    @test reduced_cell_ isa TriclinicUnitCell
    @test volume(reduced_cell_) ≈ 0.5 * volume(reference_unit_cell)
    @test reduced_cell_ ≈ expected_reduced_cell

    # equivalent body-centered and base-centered monoclinic unit cells
    base_centered_unit_cell = MonoclinicUnitCell(
        lattice_constants(reference_unit_cell); centering=base_centering
    )
    body_centering_basis_a = basis_a + basis_c
    body_centering_basis_b = basis_b
    body_centering_basis_c = basis_c
    body_centered_unit_cell = MonoclinicUnitCell(
        norm(body_centering_basis_a),
        norm(body_centering_basis_b),
        norm(body_centering_basis_c),
        π - acos(
            dot(body_centering_basis_a, body_centering_basis_c) /
            norm(body_centering_basis_a) / norm(body_centering_basis_c),
        );
        centering=body_centering,
    )

    @test base_centered_unit_cell isa MonoclinicUnitCell
    @test body_centered_unit_cell isa MonoclinicUnitCell

    @test reduced_cell(base_centered_unit_cell) ≈ reduced_cell(body_centered_unit_cell)

    # face-centered unit cell (equivalent to smaller body-centered unit cell)
    face_centered_unit_cell = MonoclinicUnitCell(
        lattice_constants(reference_unit_cell); centering=face_centering
    )
    body_centering_basis_a = 0.5 * (basis_a + basis_c)
    body_centering_basis_b = basis_b
    body_centering_basis_c = 0.5 * (basis_a - basis_c)
    body_centered_unit_cell = MonoclinicUnitCell(
        norm(body_centering_basis_a),
        norm(body_centering_basis_b),
        norm(body_centering_basis_c),
        acos(
            dot(body_centering_basis_a, body_centering_basis_c) /
            norm(body_centering_basis_a) / norm(body_centering_basis_c),
        );
        centering=body_centering,
    )

    reduced_face_centered_unit_cell = reduced_cell(face_centered_unit_cell)
    reduced_body_centered_unit_cell = reduced_cell(body_centered_unit_cell)

    @test face_centered_unit_cell isa MonoclinicUnitCell
    @test body_centered_unit_cell isa MonoclinicUnitCell

    @test reduced_face_centered_unit_cell ≈ reduced_body_centered_unit_cell
    @test volume(reduced_body_centered_unit_cell) ≈ 0.25 * volume(face_centered_unit_cell)
    @test volume(reduced_face_centered_unit_cell) ≈ 0.25 * volume(face_centered_unit_cell)
end

@testset "is_equivalent_unit_cell(::UnitCell): monoclinic" begin
    # --- Preparations

    a = 2
    b = 3
    c = 5
    β = 5π / 7
    reference_unit_cell = MonoclinicUnitCell(a, b, c, β)
    basis_a, basis_b, basis_c = basis(reference_unit_cell)

    # --- Tests

    # equivalent monoclinic and triclinic unit cells
    monoclinic_unit_cell = MonoclinicUnitCell(
        lattice_constants(reference_unit_cell); centering=primitive_centering
    )
    triclinic_unit_cell = UnitCell(
        basis_a,
        basis_b,
        basis_c;
        identify_lattice_system=false,
        centering=primitive_centering,
    )
    @test is_equivalent_unit_cell(monoclinic_unit_cell, triclinic_unit_cell)

    # body-centered unit cell
    body_centered_unit_cell = MonoclinicUnitCell(
        lattice_constants(reference_unit_cell); centering=body_centering
    )
    primitive_unit_cell = UnitCell(
        basis_a,
        basis_b,
        0.5 * (basis_a + basis_b + basis_c);
        identify_lattice_system=false,
        centering=primitive_centering,
    )
    @test is_equivalent_unit_cell(body_centered_unit_cell, primitive_unit_cell)

    # equivalent base-centered and body-centered monoclinic unit cells
    base_centered_unit_cell = MonoclinicUnitCell(
        lattice_constants(reference_unit_cell); centering=base_centering
    )
    body_centered_unit_cell = UnitCell(
        basis_a + basis_c, basis_b, basis_c; centering=body_centering
    )
    @test is_equivalent_unit_cell(base_centered_unit_cell, body_centered_unit_cell)
end

@testset "is_equivalent_unit_cell(::MonoclinicUnitCell)" begin
    # --- Preparations

    a_ref = 6
    b_ref = 10
    c_ref = 8
    β_ref = π / 3
    lattice_constants_ref = MonoclinicUnitCell(a_ref, b_ref, c_ref, β_ref)

    # --- Exercise functionality and check results

    # Equivalent unit cell #1
    a = a_ref
    b = b_ref
    c = sqrt(a_ref^2 + c_ref^2 - 2 * a_ref * c_ref * cos(β_ref))
    β = asin(sin(β_ref) / c * c_ref)
    lattice_constants_test = MonoclinicUnitCell(a, b, c, β)

    @test is_equivalent_unit_cell(lattice_constants_test, lattice_constants_ref)

    # Equivalent unit cell #2
    a = sqrt(a_ref^2 + c_ref^2 - 2 * a_ref * c_ref * cos(β_ref))
    b = b_ref
    c = c_ref
    β = asin(sin(β_ref) / a * a_ref)
    lattice_constants_test = MonoclinicUnitCell(a, b, c, β)

    @test is_equivalent_unit_cell(lattice_constants_test, lattice_constants_ref)

    # Equivalent unit cell #3
    a = a_ref
    b = b_ref
    c = sqrt(a_ref^2 + c_ref^2 + 2 * a_ref * c_ref * cos(β_ref))
    β = asin(sin(β_ref) / c * c_ref)
    lattice_constants_test = MonoclinicUnitCell(a, b, c, β)

    @test is_equivalent_unit_cell(lattice_constants_test, lattice_constants_ref)

    # Supercell: b multiple of b_ref
    a = a_ref
    b = 3 * b_ref
    c = c_ref
    β = β_ref
    lattice_constants_test = MonoclinicUnitCell(a, b, c, β)

    @test !is_equivalent_unit_cell(lattice_constants_test, lattice_constants_ref)

    # Supercell: c multiple of c_ref
    a = a_ref
    b = b_ref
    c = 5 * c_ref
    β = β_ref
    lattice_constants_test = MonoclinicUnitCell(a, b, c, β)

    @test !is_equivalent_unit_cell(lattice_constants_test, lattice_constants_ref)

    # Supercell: c equal to diagonal from unit cell twice as high in c_ref direction
    a = a_ref
    b = b_ref
    c = sqrt(a_ref^2 + (2 * c_ref)^2 - 2 * a_ref * (2 * c_ref) * cos(β_ref))
    β = asin(sin(β_ref) / c * (2 * c_ref))
    lattice_constants_test = MonoclinicUnitCell(a, b, c, β)

    @test !is_equivalent_unit_cell(lattice_constants_test, lattice_constants_ref)

    # test unit cell and reference unit cell are for different lattice systems
    lattice_constants_test = CubicUnitCell(1)
    @test !is_equivalent_unit_cell(lattice_constants_test, lattice_constants_ref)
end

@testset "is_supercell(::MonoclinicUnitCell): valid arguments" begin
    # --- Preparations

    a_ref = 6
    b_ref = 10
    c_ref = 8
    β_ref = π / 3
    lattice_constants_ref = MonoclinicUnitCell(a_ref, b_ref, c_ref, β_ref)

    # --- Exercise functionality

    # Supercell: b multiple of b_ref
    a = a_ref
    b = 3 * b_ref
    c = c_ref
    β = β_ref
    lattice_constants_test = MonoclinicUnitCell(a, b, c, β)

    @test is_supercell(lattice_constants_test, lattice_constants_ref)

    # Supercell: c multiple of c_ref
    a = a_ref
    b = b_ref
    c = 5 * c_ref
    β = β_ref
    lattice_constants_test = MonoclinicUnitCell(a, b, c, β)

    @test !is_supercell(lattice_constants_test, lattice_constants_ref)
    @test is_supercell(lattice_constants_test, lattice_constants_ref; max_index=5)

    # Supercell: c equal to diagonal from unit cell twice as high in c_ref direction
    a = a_ref
    b = b_ref
    c = sqrt(a_ref^2 + (2 * c_ref)^2 - 2 * a_ref * (2 * c_ref) * cos(β_ref))
    β = asin(sin(β_ref) / c * (2 * c_ref))
    lattice_constants_test = MonoclinicUnitCell(a, b, c, β)

    @test is_supercell(lattice_constants_test, lattice_constants_ref)

    # Supercell: basis formed from both diagonals of parallelogram formed by a_ref
    #            and c_ref
    a = sqrt(a_ref^2 + c_ref^2 - 2 * a_ref * c_ref * cos(β_ref))
    b = b_ref
    c = sqrt(a_ref^2 + c_ref^2 + 2 * a_ref * c_ref * cos(β_ref))
    β = asin(sin(β_ref) / a * a_ref) + asin(sin(β_ref) / c * a_ref)
    lattice_constants_test = MonoclinicUnitCell(a, b, c, β)

    @test is_supercell(lattice_constants_test, lattice_constants_ref)

    # Equivalent unit cell
    a = a_ref
    b = b_ref
    c = sqrt(a_ref^2 + c_ref^2 - 2 * a_ref * c_ref * cos(β_ref))
    β = asin(sin(β_ref) / c * c_ref)
    lattice_constants_test = MonoclinicUnitCell(a, b, c, β)

    @test !is_supercell(lattice_constants_test, lattice_constants_ref)
end

@testset "is_supercell(::MonoclinicUnitCell): invalid arguments" begin
    # --- Preparations

    a_ref = 6
    b_ref = 10
    c_ref = 8
    β_ref = π / 3
    lattice_constants_ref = MonoclinicUnitCell(a_ref, b_ref, c_ref, β_ref)

    a_test = a_ref
    b_test = 3 * b_ref
    c_test = c_ref
    β_test = β_ref
    lattice_constants_test = MonoclinicUnitCell(a_test, b_test, c_test, β_test)

    # --- Exercise functionality and check results

    # ------ `tol`

    # tol = 0
    expected_message = "`tol` must be positive"
    @test_throws DomainError(0, expected_message) is_supercell(
        lattice_constants_test, lattice_constants_ref; tol=0
    )

    # tol < 0
    expected_message = "`tol` must be positive"
    @test_throws DomainError(-0.1, expected_message) is_supercell(
        lattice_constants_test, lattice_constants_ref; tol=-0.1
    )

    # ------ `max_index`

    # max_index = 0
    expected_message = "`max_index` must be positive"
    @test_throws DomainError(0, expected_message) is_supercell(
        lattice_constants_test, lattice_constants_ref; max_index=0
    )

    # max_index < 0
    expected_message = "`max_index` must be positive"
    @test_throws DomainError(-3, expected_message) is_supercell(
        lattice_constants_test, lattice_constants_ref; max_index=-3
    )
end

@testset ":(==)(::MonoclinicUnitCell)" begin
    # --- Exercise functionality and check results

    # unit cells are equal
    @test MonoclinicUnitCell(1, 2, 3, 1.5) == MonoclinicUnitCell(1, 2, 3, 1.5)

    # lattice constants are not equal
    @test MonoclinicUnitCell(1, 2, 3, 1.5) != MonoclinicUnitCell(1, 2, 3, 1.6)

    # centerings are not equal
    @test (
        MonoclinicUnitCell(1, 2, 3, 1.5; centering=primitive_centering) !=
        MonoclinicUnitCell(1, 2, 3, 1.5; centering=body_centering)
    )

    # symmetry elements are not equal
    @test (
        MonoclinicUnitCell(1, 2, 3, 1.5; symmetry_elements=Set()) !=
        MonoclinicUnitCell(1, 2, 3, 1.5; symmetry_elements=Set([a_2_1, d_perp_110]))
    )
end

@testset "isapprox(::MonoclinicUnitCell)" begin
    # --- Preparations

    x = MonoclinicUnitCell(1.0, 2.0, 3.0, π / 5)
    y = MonoclinicUnitCell(1.5, 2.5, 3.5, π / 5 + 0.5)

    # --- Exercise functionality and check results

    # x ≈ (x + delta)
    @test x ≈ MonoclinicUnitCell(1.0 + 1e-9, 2.0, 3.0, π / 5)
    @test x ≈ MonoclinicUnitCell(1.0, 2.0 + 1e-9, 3.0, π / 5)
    @test x ≈ MonoclinicUnitCell(1.0, 2.0, 3.0 - 1e-9, π / 5)
    @test x ≈ MonoclinicUnitCell(1.0, 2.0, 3.0, π / 5 - 1e-9)

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

@testset "-(::MonoclinicUnitCell)" begin
    # --- Tests

    lattice_constants_x = (a=1, b=3, c=5, β=π / 3)
    x = MonoclinicUnitCell(lattice_constants_x)

    lattice_constants_y = (a=2, b=10, c=20, β=π / 6)
    y = MonoclinicUnitCell(lattice_constants_y)

    @test x - y == MonoclinicUnitCellDelta(
        lattice_constants_x.a - lattice_constants_y.a,
        lattice_constants_x.b - lattice_constants_y.b,
        lattice_constants_x.c - lattice_constants_y.c,
        lattice_constants_x.β - lattice_constants_y.β,
    )
end
