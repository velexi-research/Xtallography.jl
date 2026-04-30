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
Tests for TriclinicUnitCell types and methods in unit_cell/triclinic.jl (except for
conventional_cell() tests)
"""
# --- Imports

# Standard library
using LinearAlgebra: cross, det, norm, qr, I
using Test

# Xtallography package
using Xtallography

# Testing utilities
include("testing_utilities.jl")

# --- Tests

# ------ Constructors

@testset "TriclinicUnitCell(::NamedTuple,::UnitCellSymmetry) inner constructor" begin
    # --- Tests

    # Valid arguments
    lattice_constants_ = (a=1, b=3, c=5, α=π / 7, β=2π / 7, γ=3π / 7)
    symmetry_ = primitive_unit_cell_symmetry
    unit_cell = TriclinicUnitCell(lattice_constants_, symmetry_)

    @test lattice_constants(unit_cell) == lattice_constants_
    @test symmetry(unit_cell) == symmetry_

    # Invalid arguments
    lattice_constants_ = (a=1,)
    expected_message = (
        "Invalid lattice_constants argument passed to UnitCell{Triclinic} " *
        "constructor. Expected keys: (:a, :b, :c, :α, :β, :γ). " *
        "Provided keys: $(keys(lattice_constants_))."
    )
    @test_throws ArgumentError(expected_message) TriclinicUnitCell(
        lattice_constants_, symmetry_
    )
end

@testset "TriclinicUnitCell(::NamedTuple;::Centering,::Set) outer constructor" begin
    # --- Preparations

    lattice_constants_ = (a=1, b=3, c=5, α=π / 7, β=2π / 7, γ=3π / 7)

    # --- Tests

    # Default keyword arguments
    unit_cell = TriclinicUnitCell(lattice_constants_)

    @test lattice_constants(unit_cell) == lattice_constants_
    @test symmetry(unit_cell) == primitive_unit_cell_symmetry

    # Non-default centering keyword argument
    for centering_ in CENTERINGS
        unit_cell = TriclinicUnitCell(lattice_constants_; centering=centering_)

        @test lattice_constants(unit_cell) == lattice_constants_
        @test symmetry(unit_cell) == UnitCellSymmetry(centering_, Set{SymmetryElement}())
    end

    # Non-default symmetry_elements keyword argument
    symmetry_elements_ = [a_4_2]
    unit_cell = TriclinicUnitCell(lattice_constants_; symmetry_elements=symmetry_elements_)

    @test lattice_constants(unit_cell) == lattice_constants_
    @test symmetry(unit_cell) ==
        UnitCellSymmetry(primitive_centering, Set{SymmetryElement}(symmetry_elements_))
end

@testset "TriclinicUnitCell(::Real,::Real,::Real,::Real,::Real,::Real;::Centering,::Set) outer constructor: valid arguments" begin
    # --- Preparations

    a = 1
    b = 2
    c = 3
    α = π / 4
    β = π / 2
    γ = 3 * π / 4

    # --- Exercise functionality and check results

    # Default keyword arguments
    unit_cell = TriclinicUnitCell(a, b, c, α, β, γ)

    @test lattice_constants(unit_cell).a == a
    @test lattice_constants(unit_cell).b == b
    @test lattice_constants(unit_cell).c == c
    @test lattice_constants(unit_cell).α == α
    @test lattice_constants(unit_cell).β == β
    @test lattice_constants(unit_cell).γ == γ

    @test centering(unit_cell) === primitive_centering
    @test symmetry_elements(unit_cell) isa Set
    @test isempty(symmetry_elements(unit_cell))

    # Non-default centering keyword argument
    for centering_ in CENTERINGS
        unit_cell = TriclinicUnitCell(a, b, c, α, β, γ; centering=centering_)

        @test lattice_constants(unit_cell).a == a
        @test lattice_constants(unit_cell).b == b
        @test lattice_constants(unit_cell).c == c
        @test lattice_constants(unit_cell).α == α
        @test lattice_constants(unit_cell).β == β
        @test lattice_constants(unit_cell).γ == γ
        @test centering(unit_cell) === centering_
        @test symmetry_elements(unit_cell) isa Set
        @test isempty(symmetry_elements(unit_cell))
    end

    # Non-default symmetry_elements keyword argument
    symmetry_elements_ = [a_4_2]
    unit_cell = TriclinicUnitCell(a, b, c, α, β, γ; symmetry_elements=symmetry_elements_)

    @test lattice_constants(unit_cell).a == a
    @test lattice_constants(unit_cell).b == b
    @test lattice_constants(unit_cell).c == c
    @test lattice_constants(unit_cell).α == α
    @test lattice_constants(unit_cell).β == β
    @test lattice_constants(unit_cell).γ == γ
    @test centering(unit_cell) === primitive_centering
    @test symmetry_elements(unit_cell) isa Set
    @test symmetry_elements(unit_cell) == Set(symmetry_elements_)
end

@testset "TriclinicUnitCell(::Real,::Real,::Real,::Real,::Real,::Real;::Centering,::Set) outer constructor: invalid arguments" begin
    # --- Preparations

    # Valid arguments
    a = 1
    b = 2
    c = 3
    α = π / 4
    β = π / 2
    γ = 3 * π / 4

    # --- Tests

    # ------ a

    # a = 0
    expected_message = "`a` must be positive"
    @test_throws DomainError(0, expected_message) TriclinicUnitCell(0, b, c, α, β, γ)

    # a < 0
    expected_message = "`a` must be positive"
    @test_throws DomainError(-10, expected_message) TriclinicUnitCell(-10, b, c, α, β, γ)

    # ------ b

    # b = 0
    expected_message = "`b` must be positive"
    @test_throws DomainError(0, expected_message) TriclinicUnitCell(a, 0, c, α, β, γ)

    # b < 0
    expected_message = "`b` must be positive"
    @test_throws DomainError(-1.0, expected_message) TriclinicUnitCell(a, -1.0, c, α, β, γ)

    # ------ c

    # c = 0
    expected_message = "`c` must be positive"
    @test_throws DomainError(0, expected_message) TriclinicUnitCell(a, b, 0, α, β, γ)

    # c < 0
    expected_message = "`c` must be positive"
    @test_throws DomainError(-1.0, expected_message) TriclinicUnitCell(a, b, -1.0, α, β, γ)

    # ------ α

    # α = 0
    expected_message = "`α` must satisfy 0 < α < π"
    @test_throws DomainError(0, expected_message) TriclinicUnitCell(a, b, c, 0, β, γ)

    # α < 0
    expected_message = "`α` must satisfy 0 < α < π"
    @test_throws DomainError(-1, expected_message) TriclinicUnitCell(a, b, c, -1, β, γ)

    # α > π
    expected_message = "`α` must satisfy 0 < α < π"
    @test_throws DomainError(4, expected_message) TriclinicUnitCell(a, b, c, 4, β, γ)

    # α = π
    expected_message = "`α` must satisfy 0 < α < π"
    @test_throws DomainError(π, expected_message) TriclinicUnitCell(a, b, c, π, β, γ)

    # ------ β

    # β = 0
    expected_message = "`β` must satisfy 0 < β < π"
    @test_throws DomainError(0, expected_message) TriclinicUnitCell(a, b, c, α, 0, γ)

    # β < 0
    expected_message = "`β` must satisfy 0 < β < π"
    @test_throws DomainError(-1, expected_message) TriclinicUnitCell(a, b, c, α, -1, γ)

    # β = π
    expected_message = "`β` must satisfy 0 < β < π"
    @test_throws DomainError(π, expected_message) TriclinicUnitCell(a, b, c, α, π, γ)

    # β > π
    expected_message = "`β` must satisfy 0 < β < π"
    @test_throws DomainError(4, expected_message) TriclinicUnitCell(a, b, c, α, 4, γ)

    # ------ γ

    # γ = 0
    expected_message = "`γ` must satisfy 0 < γ < π"
    @test_throws DomainError(0, expected_message) TriclinicUnitCell(a, b, c, α, β, 0)

    # γ < 0
    expected_message = "`γ` must satisfy 0 < γ < π"
    @test_throws DomainError(-1, expected_message) TriclinicUnitCell(a, b, c, α, β, -1)

    # γ = π
    expected_message = "`γ` must satisfy 0 < γ < π"
    @test_throws DomainError(π, expected_message) TriclinicUnitCell(a, b, c, α, β, π)

    # γ > π
    expected_message = "`γ` must satisfy 0 < γ < π"
    @test_throws DomainError(4, expected_message) TriclinicUnitCell(a, b, c, α, β, 4)

    # ------ symmetry elements

    # TODO
end

@testset "UnitCell(::UnitCell) copy constructor: triclinic" begin
    # --- Preparations

    lattice_constants_ = (a=1, b=2, c=3, α=π / 5, β=2π / 5, γ=3π / 5)
    centering = primitive_centering
    symmetry_elements_ = [a_4_2]
    reference_unit_cell = TriclinicUnitCell(
        lattice_constants_,
        UnitCellSymmetry(centering; symmetry_elements=symmetry_elements_),
    )

    # --- Tests

    unit_cell = UnitCell(reference_unit_cell)

    @test unit_cell === reference_unit_cell
end

@testset "UnitCell(::Vector,::Vector,::Vector;::Bool,::Centering) outer constructor: triclinic" begin
    # --- Preparations

    # Generate random rotation matrix
    rotations = [I, qr(rand(3)).Q, qr(rand(3)).Q]

    # --- Tests

    # ------ General triclinic unit cell

    a = 5
    b = 8
    c = 10
    α = 2π / 5
    β = 3π / 5
    γ = 4π / 5

    basis_a = [a, 0, 0]
    basis_b = [b * cos(γ), b * sin(γ), 0]
    basis_c = [
        c * cos(β),
        c / sin(γ) * (cos(α) - cos(β) * cos(γ)),
        2 * c / sin(γ) * sqrt(
            sin(0.5 * (α + β + γ)) *
            sin(0.5 * (α - β + γ)) *
            sin(0.5 * (α + β - γ)) *
            sin(0.5 * (-α + β + γ)),
        ),
    ]

    expected_unit_cell = standardize(TriclinicUnitCell(a, b, c, α, β, γ))

    # default centering keyword argument
    test_basis_rotations_and_permutations(
        rotations, expected_unit_cell, basis_a, basis_b, basis_c
    )

    # centering keyword argument provided
    for centering_ in (primitive_centering,)
        test_basis_rotations_and_permutations(
            rotations, expected_unit_cell, basis_a, basis_b, basis_c; centering=centering_
        )
    end

    # ------ Almost rhomobohedral: all three angles equal, edge lengths unequal

    a = 5
    b = 6
    c = 7
    α = 2π / 5

    basis_a = [a, 0, 0]
    basis_b = [b * cos(α), b * sin(α), 0]
    basis_c = [
        c * cos(α),
        c / sin(α) * (cos(α) - cos(α)^2),
        c / sin(α) * sqrt(1 - 3 * cos(α)^2 + 2 * cos(α)^3),
    ]

    expected_unit_cell = standardize(TriclinicUnitCell(a, b, c, α, α, α))

    # default centering keyword argument
    test_basis_rotations_and_permutations(
        rotations, expected_unit_cell, basis_a, basis_b, basis_c
    )

    # centering keyword argument provided
    for centering_ in (primitive_centering,)
        test_basis_rotations_and_permutations(
            rotations, expected_unit_cell, basis_a, basis_b, basis_c; centering=centering_
        )
    end
end

# ------ Methods

@testset "lattice_system(::TriclinicUnitCell)" begin
    lattice_constants_ = (a=1, b=2, c=3, α=π / 5, β=2π / 5, γ=3π / 5)
    unit_cell = TriclinicUnitCell(lattice_constants_)
    @test lattice_system(unit_cell) === triclinic
end

@testset "lattice_constants(::TriclinicUnitCell)" begin
    lattice_constants_ = (a=1, b=2, c=3, α=π / 4, β=π / 2, γ=3 * π / 4)
    unit_cell = TriclinicUnitCell(lattice_constants_)
    @test lattice_constants(unit_cell) == lattice_constants_
end

@testset "symmetry(::TriclinicUnitCell)" begin
    # --- Preparations

    lattice_constants_ = (a=1, b=2, c=3, α=π / 4, β=π / 2, γ=3 * π / 4)
    symmetry_elements_ = Set([a_2_1, d_perp_110])

    # --- Tests

    for centering_ in CENTERINGS
        symmetry_ = UnitCellSymmetry(;
            centering=centering_, symmetry_elements=symmetry_elements_
        )
        unit_cell = TriclinicUnitCell(lattice_constants_, symmetry_)

        @test symmetry(unit_cell) == symmetry_
    end
end

@testset "centering(::TriclinicUnitCell)" begin
    # --- Preparations

    lattice_constants_ = (a=1, b=2, c=3, α=π / 4, β=π / 2, γ=3 * π / 4)

    # --- Tests

    # Default centering
    unit_cell = TriclinicUnitCell(lattice_constants_)
    @test centering(unit_cell) == primitive_centering

    # Non-default centering
    for centering_ in CENTERINGS
        unit_cell = TriclinicUnitCell(lattice_constants_; centering=centering_)
        @test centering(unit_cell) === centering_
    end
end

@testset "symmetry_elements(::TriclinicUnitCell)" begin
    # --- Preparations

    lattice_constants_ = (a=1, b=2, c=3, α=π / 4, β=π / 2, γ=3 * π / 4)

    # --- Tests

    # Default symmetry elements
    unit_cell = TriclinicUnitCell(lattice_constants_)
    @test symmetry_elements(unit_cell) isa Set{SymmetryElement}
    @test isempty(symmetry_elements(unit_cell))

    # Non-default symmetry elements
    symmetry_elements_ = Set([a_2_1, d_perp_110])
    unit_cell = TriclinicUnitCell(lattice_constants_; symmetry_elements=symmetry_elements_)
    @test symmetry_elements(unit_cell) == symmetry_elements_
end

@testset "is_bravais_lattice(::TriclinicUnitCell)" begin
    # --- Preparations

    a = 1
    b = 2
    c = 3
    α = π / 4
    β = π / 2
    γ = 3 * π / 4

    # --- Tests

    # Valid Bravais lattices
    unit_cell = TriclinicUnitCell(a, b, c, α, β, γ; centering=primitive_centering)
    @test is_bravais_lattice(unit_cell)

    # Invalid Bravais lattices
    for centering_ in (body_centering, face_centering, base_centering)
        unit_cell = TriclinicUnitCell(a, b, c, α, β, γ; centering=centering_)
        @test !is_bravais_lattice(unit_cell)
    end
end

@testset "basis(::TriclinicUnitCell)" begin
    # --- Preparations

    a = 2
    b = 3
    c = 5
    α = π / 3
    β = π / 4
    γ = π / 5
    unit_cell = TriclinicUnitCell(a, b, c, α, β, γ)

    # --- Exercise functionality

    basis_a, basis_b, basis_c = basis(unit_cell)

    # --- Check results

    @test basis_a ≈ [a, 0, 0]
    @test basis_b ≈ [b * cos(γ), b * sin(γ), 0]

    V = volume(unit_cell)
    @test basis_c ≈
        [c * cos(β), c / sin(γ) * (cos(α) - cos(β) * cos(γ)), V / sin(γ) / a / b]
end

@testset "volume(::TricinicUnitCell)" begin
    # --- Preparations

    # Construct basis vectors for unit cell
    a = 2
    b = 3
    c = 5
    α = π / 3
    β = π / 4
    γ = π / 5
    unit_cell = TriclinicUnitCell(a, b, c, α, β, γ)

    basis_a, basis_b, basis_c = basis(unit_cell)

    # --- Exercise functionality and check results

    @test volume(unit_cell) ≈ abs(det(hcat(basis_a, basis_b, basis_c)))
end

@testset "surface_area(::TricinicUnitCell)" begin
    # --- Preparations

    # Construct basis vectors for unit cell
    a = 2
    b = 3
    c = 5
    α = π / 3
    β = π / 4
    γ = π / 5
    unit_cell = TriclinicUnitCell(a, b, c, α, β, γ)

    basis_a, basis_b, basis_c = basis(unit_cell)

    # --- Exercise functionality and check results

    @test surface_area(unit_cell) ≈
        2 * norm(cross(basis_a, basis_b)) +
          2 * norm(cross(basis_b, basis_c)) +
          2 * norm(cross(basis_c, basis_a))
end

@testset "standardize(::TriclinicUnitCell): Type I cell" begin
    # --- Tests

    # ------ lattice constants already in standard form

    a = 1.0
    b = 5.0
    c = 10.0
    α = 1π / 10
    β = 2π / 10
    γ = 3π / 10
    unit_cell = TriclinicUnitCell(a, b, c, α, β, γ)
    @test is_triclinic_type_I_cell(unit_cell)

    # Exercise functionality and check results
    standardized_unit_cell = standardize(unit_cell)

    expected_unit_cell = TriclinicUnitCell(a, b, c, α, β, γ)
    @test standardized_unit_cell ≈ expected_unit_cell

    # ------ edge lengths not sorted

    a = 1.0
    b = 10.0
    c = 5.0
    α = 1π / 10
    β = 2π / 10
    γ = 3π / 10
    unit_cell = TriclinicUnitCell(a, b, c, α, β, γ)
    @test is_triclinic_type_I_cell(unit_cell)

    # Exercise functionality and check results
    standardized_unit_cell = standardize(unit_cell)

    expected_unit_cell = TriclinicUnitCell(a, c, b, α, γ, β)
    @test standardized_unit_cell ≈ expected_unit_cell

    # ------ origin not at "homogeneous corner"

    a = 1.0
    b = 5.0
    c = 10.0
    α = 1π / 10
    β = 8π / 10
    γ = 7π / 10
    unit_cell = TriclinicUnitCell(a, b, c, α, β, γ)
    @test is_triclinic_type_I_cell(unit_cell)

    # Exercise functionality and check results
    standardized_unit_cell = standardize(unit_cell)

    expected_unit_cell = TriclinicUnitCell(a, b, c, α, π - β, π - γ)
    @test standardized_unit_cell ≈ expected_unit_cell

    # ------ edge lengths not sorted and origin not at "homogeneous corner"

    a = 10.0
    b = 1.0
    c = 5.0
    α = 1π / 10
    β = 8π / 10
    γ = 7π / 10
    unit_cell = TriclinicUnitCell(a, b, c, α, β, γ)
    @test is_triclinic_type_I_cell(unit_cell)

    # Exercise functionality and check results
    standardized_unit_cell = standardize(unit_cell)

    expected_unit_cell = TriclinicUnitCell(b, c, a, π - β, π - γ, α)
    @test standardized_unit_cell ≈ expected_unit_cell

    # ------ a ≈ b ≈ c, only need to swap α and β

    a = b = c = 1.0
    α = 2π / 10
    β = 1π / 10
    γ = 3π / 10
    unit_cell = TriclinicUnitCell(a, b, c, α, β, γ)
    @test is_triclinic_type_I_cell(unit_cell)

    # Exercise functionality and check results
    standardized_unit_cell = standardize(unit_cell)

    expected_unit_cell = TriclinicUnitCell(a, b, c, β, α, γ)
    @test standardized_unit_cell ≈ expected_unit_cell

    # ------ a ≈ b ≈ c, only need to swap β and γ

    a = b = c = 1.0
    α = 1π / 10
    β = 3π / 10
    γ = 2π / 10
    unit_cell = TriclinicUnitCell(a, b, c, α, β, γ)
    @test is_triclinic_type_I_cell(unit_cell)

    # Exercise functionality and check results
    standardized_unit_cell = standardize(unit_cell)

    expected_unit_cell = TriclinicUnitCell(a, b, c, α, γ, β)
    @test standardized_unit_cell ≈ expected_unit_cell

    # ------ a ≈ b ≈ c, only need to reverse order of angles

    a = b = c = 1.0
    α = 3π / 10
    β = 2π / 10
    γ = 1π / 10
    unit_cell = TriclinicUnitCell(a, b, c, α, β, γ)
    @test is_triclinic_type_I_cell(unit_cell)

    # Exercise functionality and check results
    standardized_unit_cell = standardize(unit_cell)

    expected_unit_cell = TriclinicUnitCell(a, b, c, γ, β, α)
    @test standardized_unit_cell ≈ expected_unit_cell

    # ------ a ≈ b (after sorting a, b, c)

    a = b = 1.0
    c = 5
    α = 4π / 10
    β = 8π / 10
    γ = 7π / 10
    unit_cell = TriclinicUnitCell(a, b, c, α, β, γ)
    @test is_triclinic_type_I_cell(unit_cell)

    # Exercise functionality and check results
    standardized_unit_cell = standardize(unit_cell)

    expected_unit_cell = TriclinicUnitCell(a, b, c, π - β, α, π - γ)
    @test standardized_unit_cell ≈ expected_unit_cell

    # ------ b ≈ c (after sorting a, b, c)

    a = b = 5.0
    c = 1
    α = 8π / 10
    β = 7π / 10
    γ = 4π / 10
    unit_cell = TriclinicUnitCell(a, b, c, α, β, γ)
    @test is_triclinic_type_I_cell(unit_cell)

    # Exercise functionality and check results
    standardized_unit_cell = standardize(unit_cell)

    expected_unit_cell = TriclinicUnitCell(c, a, b, γ, π - α, π - β)
    @test standardized_unit_cell ≈ expected_unit_cell
end

@testset "standardize(::TriclinicUnitCell): Type II cell" begin
    # --- Tests

    # ------ lattice constants already in standard form

    a = 1.0
    b = 5.0
    c = 10.0
    α = 6π / 10
    β = 7π / 10
    γ = 8π / 10
    unit_cell = TriclinicUnitCell(a, b, c, α, β, γ)
    @test !is_triclinic_type_I_cell(unit_cell)

    # Exercise functionality and check results
    standardized_unit_cell = standardize(unit_cell)

    expected_unit_cell = unit_cell
    @test standardized_unit_cell ≈ expected_unit_cell

    # ------ edge lengths not sorted

    a = 1.0
    b = 10.0
    c = 5.0
    α = 6π / 10
    β = 7π / 10
    γ = 8π / 10
    unit_cell = TriclinicUnitCell(a, b, c, α, β, γ)
    @test !is_triclinic_type_I_cell(unit_cell)

    # Exercise functionality and check results
    standardized_unit_cell = standardize(unit_cell)

    expected_unit_cell = TriclinicUnitCell(a, c, b, α, γ, β)
    @test standardized_unit_cell ≈ expected_unit_cell

    # ------ origin not at "homogeneous corner"

    a = 1.0
    b = 5.0
    c = 10.0
    α = 4π / 10
    β = 3π / 10
    γ = 8π / 10
    unit_cell = TriclinicUnitCell(a, b, c, α, β, γ)
    @test !is_triclinic_type_I_cell(unit_cell)

    # Exercise functionality and check results
    standardized_unit_cell = standardize(unit_cell)

    expected_unit_cell = TriclinicUnitCell(a, b, c, π - α, π - β, γ)
    @test standardized_unit_cell ≈ expected_unit_cell

    # ------ edge lengths not sorted and origin not at "homogeneous corner"

    a = 10.0
    b = 1.0
    c = 5.0
    α = 4π / 10
    β = 3π / 10
    γ = 8π / 10
    unit_cell = TriclinicUnitCell(a, b, c, α, β, γ)
    @test !is_triclinic_type_I_cell(unit_cell)

    # Exercise functionality and check results
    standardized_unit_cell = standardize(unit_cell)

    expected_unit_cell = TriclinicUnitCell(b, c, a, π - β, γ, π - α)
    @test standardized_unit_cell ≈ expected_unit_cell

    # ------ a ≈ b ≈ c, only need to swap α and β

    a = b = c = 1.0
    α = 8π / 10
    β = 7π / 10
    γ = 9π / 10
    unit_cell = TriclinicUnitCell(a, b, c, α, β, γ)
    @test !is_triclinic_type_I_cell(unit_cell)

    # Exercise functionality and check results
    standardized_unit_cell = standardize(unit_cell)

    expected_unit_cell = TriclinicUnitCell(a, b, c, β, α, γ)
    @test standardized_unit_cell ≈ expected_unit_cell

    # ------ a ≈ b ≈ c, only need to swap β and γ

    a = b = c = 1.0
    α = 7π / 10
    β = 9π / 10
    γ = 8π / 10
    unit_cell = TriclinicUnitCell(a, b, c, α, β, γ)
    @test !is_triclinic_type_I_cell(unit_cell)

    # Exercise functionality and check results
    standardized_unit_cell = standardize(unit_cell)

    expected_unit_cell = TriclinicUnitCell(a, b, c, α, γ, β)
    @test standardized_unit_cell ≈ expected_unit_cell

    # ------ a ≈ b ≈ c, only need to reverse order of angles

    a = b = c = 1.0
    α = 9π / 10
    β = 8π / 10
    γ = 7π / 10
    unit_cell = TriclinicUnitCell(a, b, c, α, β, γ)
    @test !is_triclinic_type_I_cell(unit_cell)

    # Exercise functionality and check results
    standardized_unit_cell = standardize(unit_cell)

    expected_unit_cell = TriclinicUnitCell(a, b, c, γ, β, α)
    @test standardized_unit_cell ≈ expected_unit_cell

    # ------ a ≈ b (after sorting a, b, c)

    a = b = 1.0
    c = 5
    α = 7π / 10
    β = 4π / 10
    γ = 4π / 10
    unit_cell = TriclinicUnitCell(a, b, c, α, β, γ)
    @test !is_triclinic_type_I_cell(unit_cell)

    # Exercise functionality and check results
    standardized_unit_cell = standardize(unit_cell)

    expected_unit_cell = TriclinicUnitCell(a, b, c, π - β, α, π - γ)
    @test standardized_unit_cell ≈ expected_unit_cell

    # ------ b ≈ c (after sorting a, b, c)

    a = b = 5.0
    c = 1
    α = 8π / 10
    β = 4π / 10
    γ = 3π / 10
    unit_cell = TriclinicUnitCell(a, b, c, α, β, γ)
    @test !is_triclinic_type_I_cell(unit_cell)

    # Exercise functionality and check results
    standardized_unit_cell = standardize(unit_cell)

    expected_unit_cell = TriclinicUnitCell(c, a, b, π - γ, π - β, α)
    @test standardized_unit_cell ≈ expected_unit_cell
end

@testset "standardize(::TriclinicUnitCell): invalid arguments" begin
    # --- Tests

    for centering_ in (body_centering, face_centering, base_centering)
        unit_cell = TriclinicUnitCell(1, 2, 3, 2π / 5, 3π / 5, 4π / 5; centering=centering_)

        expected_message =
            "Invalid Bravais lattice: " *
            "(lattice_system=Triclinic, centering=$(nameof(typeof(centering_))))"
        @test_throws ArgumentError(expected_message) standardize(unit_cell)
    end
end

@testset "reduced_cell(::TriclinicUnitCell)" begin
    # --- Preparations

    a = 5
    b = 8
    c = 10
    α = 2π / 5
    β = 3π / 5
    γ = 4π / 5
    reference_unit_cell = TriclinicUnitCell(a, b, c, α, β, γ)
    basis_a, basis_b, basis_c = basis(reference_unit_cell)

    # --- Exercise functionality and check results

    # primitive unit cell defined by [basis_a, basis_b, basis_c]
    unit_cell = TriclinicUnitCell(
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
    @test reduced_cell_ isa TriclinicUnitCell
    @test volume(reduced_cell_) ≈ volume(reference_unit_cell)
    @test reduced_cell_ ≈ expected_reduced_cell

    # primitive unit cell defined by linear combination of [basis_a, basis_b, basis_c]
    unit_cell = UnitCell(
        basis_a + basis_b + basis_c,
        basis_b,
        basis_c;
        identify_lattice_system=true,
        centering=primitive_centering,
    )

    expected_reduced_cell = reduced_cell(reference_unit_cell)

    reduced_cell_ = reduced_cell(unit_cell)
    @test reduced_cell_ isa TriclinicUnitCell
    @test volume(reduced_cell_) ≈ volume(reference_unit_cell)
    @test reduced_cell_ ≈ expected_reduced_cell

    #=
    # TODO: Delaunay reduction does not terminate
    a = 5
    b = 8
    c = 10
    α = π / 5
    β = 2π / 5
    γ = 3π / 5
    unit_cell = TriclinicUnitCell(a, b, c, α, β, γ)
    basis_a, basis_b, basis_c = basis(unit_cell)
    unit_cell = UnitCell(
        UnitCell(basis_a + basis_b + basis_c, basis_b, basis_c), primitive_centering
    )

    expected_reduced_cell = reduced_cell(unit_cell)

    reduced_cell_ = reduced_cell(unit_cell)
    =#
end

@testset "is_equivalent(::UnitCell): triclinic" begin
    # --- Preparations

    a = 2
    b = 3
    c = 5
    α = 3π / 7
    β = 4π / 7
    γ = 5π / 7
    unit_cell = TriclinicUnitCell(a, b, c, α, β, γ)
    basis_a, basis_b, basis_c = basis(unit_cell)

    # --- Tests

    # identical triclinic unit cells
    unit_cell_ref = unit_cell
    unit_cell_test = UnitCell(basis_a, basis_b, basis_c; identify_lattice_system=false)
    @test is_equivalent(unit_cell_test, unit_cell_ref)

    # primitive unit cell defined by linear combination of [basis_a, basis_b, basis_c]
    unit_cell_ref = reduced_cell(unit_cell)
    unit_cell_test = UnitCell(basis_a + basis_b + basis_c, basis_b, basis_c)
    @test is_equivalent(unit_cell_test, unit_cell_ref)

    # body-centered unit cell
    body_centered_unit_cell = TriclinicUnitCell(
        lattice_constants(unit_cell); centering=body_centering
    )
    primitive_unit_cell = UnitCell(
        basis_a, basis_b, 0.5 * (basis_a + basis_b + basis_c); identify_lattice_system=false
    )
    @test is_equivalent(body_centered_unit_cell, primitive_unit_cell)
end

@testset "is_equivalent(::TriclinicUnitCell)" begin
    # --- Preparations

    a_ref = 6
    b_ref = 10
    c_ref = 8
    α_ref = π / 3
    β_ref = π / 4
    γ_ref = π / 5
    reference_unit_cell = TriclinicUnitCell(a_ref, b_ref, c_ref, α_ref, β_ref, γ_ref)
    basis_a, basis_b, basis_c = basis(reference_unit_cell)

    # --- Exercise functionality and check results

    # equivalent primitive unit cell defined by linear combination of
    # [basis_a, basis_b, basis_c]
    unit_cell = UnitCell(basis_a + basis_b + basis_c, basis_b, basis_c)
    @test is_equivalent(unit_cell, reference_unit_cell)

    # equivalent primitive unit cell defined by a different linear combination of
    # [basis_a, basis_b, basis_c]
    unit_cell = UnitCell(basis_a + basis_b + basis_c, basis_b + basis_c, basis_c)
    @test is_equivalent(unit_cell, reference_unit_cell)

    # test unit cell and reference unit cell are for different lattice systems
    unit_cell = CubicUnitCell(1)
    @test !is_equivalent(unit_cell, reference_unit_cell)
end

@testset ":(==)(::TriclinicUnitCell)" begin
    # --- Exercise functionality and check results

    # unit cells are equal
    @test TriclinicUnitCell(1, 2, 3, 1.5, 1.5, 1.5) ==
        TriclinicUnitCell(1, 2, 3, 1.5, 1.5, 1.5)

    # lattice constants are not equal
    @test TriclinicUnitCell(1, 2, 3, 1.5, 1.5, 1.5) !=
        TriclinicUnitCell(1, 2, 3, 1.5, 1.5, 1.6)

    # centerings are not equal
    @test (
        TriclinicUnitCell(1, 2, 3, 1.5, 1.5, 1.5; centering=primitive_centering) !=
        TriclinicUnitCell(1, 2, 3, 1.5, 1.5, 1.5; centering=body_centering)
    )

    # symmetry elements are not equal
    @test (
        TriclinicUnitCell(1, 2, 3, 1.5, 1.5, 1.5; symmetry_elements=Set()) !=
        TriclinicUnitCell(
            1, 2, 3, 1.5, 1.5, 1.5; symmetry_elements=Set([a_2_1, d_perp_110])
        )
    )
end

@testset "isapprox(::TriclinicUnitCell)" begin
    # --- Preparations

    x = TriclinicUnitCell(1.0, 2.0, 3.0, π / 5, π / 4, 2π / 5)
    y = TriclinicUnitCell(1.5, 2.5, 3.5, π / 5 + 0.5, π / 4 + 0.5, 2π / 5 + 0.5)

    # --- Exercise functionality and check results

    # x ≈ (x + delta)
    @test x ≈ TriclinicUnitCell(1.0 + 1e-9, 2.0, 3.0, π / 5, π / 4, 2π / 5)
    @test x ≈ TriclinicUnitCell(1.0, 2.0 + 1e-9, 3.0, π / 5, π / 4, 2π / 5)
    @test x ≈ TriclinicUnitCell(1.0, 2.0, 3.0 - 1e-9, π / 5, π / 4, 2π / 5)
    @test x ≈ TriclinicUnitCell(1.0, 2.0, 3.0, π / 5 - 1e-9, π / 4, 2π / 5)
    @test x ≈ TriclinicUnitCell(1.0, 2.0, 3.0, π / 5, π / 4 + 1e-9, 2π / 5)
    @test x ≈ TriclinicUnitCell(1.0, 2.0, 3.0, π / 5, π / 4, 2π / 5 - 1e-9)

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

@testset "-(::TriclinicUnitCell)" begin
    # --- Preparations

    lattice_constants_x = (a=1, b=5, c=3, α=π / 7, β=2π / 7, γ=3π / 7)
    x = TriclinicUnitCell(lattice_constants_x)

    lattice_constants_y = (a=2, b=2.3, c=3, α=5π / 9, β=6π / 9, γ=7π / 9)
    y = TriclinicUnitCell(lattice_constants_y)

    # --- Tests

    @test x - y == TriclinicUnitCellDelta(
        lattice_constants_x.a - lattice_constants_y.a,
        lattice_constants_x.b - lattice_constants_y.b,
        lattice_constants_x.c - lattice_constants_y.c,
        lattice_constants_x.α - lattice_constants_y.α,
        lattice_constants_x.β - lattice_constants_y.β,
        lattice_constants_x.γ - lattice_constants_y.γ,
    )
end

@testset "is_triclinic_type_I_cell()" begin
    # --- Tests

    # ------ Type I

    # no angles equal to π/2
    a = 1
    b = 1
    c = 1
    α = 2π / 10
    β = 3π / 10
    γ = 4π / 10
    @test is_triclinic_type_I_cell(TriclinicUnitCell(a, b, c, α, β, γ))

    # α equal to π/2
    a = 1
    b = 1
    c = 1
    α = π / 2
    β = 3π / 10
    γ = 4π / 10
    @test is_triclinic_type_I_cell(TriclinicUnitCell(a, b, c, α, β, γ))

    # α equal to π/2 + ϵ (so that cos(α) < 0)
    a = 1
    b = 1
    c = 1
    α = π / 2 + 1e-9
    β = 3π / 10
    γ = 4π / 10
    @test is_triclinic_type_I_cell(TriclinicUnitCell(a, b, c, α, β, γ))

    # ------ Type II

    # Type II: no angles equal to π/2
    a = 1
    b = 1
    c = 1
    α = 6π / 10
    β = 7π / 10
    γ = 8π / 10
    @test !is_triclinic_type_I_cell(TriclinicUnitCell(a, b, c, α, β, γ))

    # γ equal to π/2
    a = 1
    b = 1
    c = 1
    α = 8π / 10
    β = 7π / 10
    γ = π / 2
    @test is_triclinic_type_I_cell(TriclinicUnitCell(a, b, c, α, β, γ))

    # β equal to π/2 + ϵ (so that cos(β) < 0)
    a = 1
    b = 1
    c = 1
    α = 7π / 10
    β = π / 2 + 1e-9
    γ = 8π / 10
    @test is_triclinic_type_I_cell(TriclinicUnitCell(a, b, c, α, β, γ))
end

@testset "satisfies_triclinic_angle_constraints()" begin
    # --- Tests

    # valid angles for a triclinic unit cell
    @test satisfies_triclinic_angle_constraints(π / 4, π / 5, π / 6)

    # invalid angles for a triclinic unit cell
    @test !satisfies_triclinic_angle_constraints(3π / 4, 4π / 5, 5π / 6)
end
