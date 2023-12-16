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
Tests for unit cell standardization methods for triclinic lattices
"""
# --- Imports

# Standard library
using Test
using LinearAlgebra: norm, dot

# External packages
using Combinatorics: combinations, permutations
using AngleBetweenVectors: angle

# XtallographyUtils package
using XtallographyUtils

# Notes
# =====
# These tests adopt the following variable conventions.
#
# - Unless otherwise noted, lattice constants and basis vectors refer to the triclinic
#   (not monoclinic) unit cell.
#
# - Lattice constants and basis vectors for the monoclinic unit cell are indicated by
#   the "m_" prefix.

# --- Tests

# ------ standardize()

@testset "standardize(): Type I cell" begin
    # --- Tests

    # ------ lattice constants already in standard form

    a = 1.0
    b = 5.0
    c = 10.0
    α = 1π / 10
    β = 2π / 10
    γ = 3π / 10
    lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)
    @test is_triclinic_type_I_cell(lattice_constants)

    # Exercise functionality and check results
    standardized_lattice_constants, standardized_centering = standardize(
        lattice_constants, XtallographyUtils.PRIMITIVE
    )

    expected_lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)
    @test standardized_lattice_constants ≈ expected_lattice_constants

    @test standardized_centering == XtallographyUtils.PRIMITIVE

    # ------ edge lengths not sorted

    a = 1.0
    b = 10.0
    c = 5.0
    α = 1π / 10
    β = 2π / 10
    γ = 3π / 10
    lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)
    @test is_triclinic_type_I_cell(lattice_constants)

    # Exercise functionality and check results
    standardized_lattice_constants, standardized_centering = standardize(
        lattice_constants, XtallographyUtils.PRIMITIVE
    )

    expected_lattice_constants = TriclinicLatticeConstants(a, c, b, α, γ, β)
    @test standardized_lattice_constants ≈ expected_lattice_constants

    @test standardized_centering == XtallographyUtils.PRIMITIVE

    # ------ origin not at "homogeneous corner"

    a = 1.0
    b = 5.0
    c = 10.0
    α = 1π / 10
    β = 8π / 10
    γ = 7π / 10
    lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)
    @test is_triclinic_type_I_cell(lattice_constants)

    # Exercise functionality and check results
    standardized_lattice_constants, standardized_centering = standardize(
        lattice_constants, XtallographyUtils.PRIMITIVE
    )

    expected_lattice_constants = TriclinicLatticeConstants(a, b, c, α, π - β, π - γ)
    @test standardized_lattice_constants ≈ expected_lattice_constants

    @test standardized_centering == XtallographyUtils.PRIMITIVE

    # ------ edge lengths not sorted and origin not at "homogeneous corner"

    a = 10.0
    b = 1.0
    c = 5.0
    α = 1π / 10
    β = 8π / 10
    γ = 7π / 10
    lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)
    @test is_triclinic_type_I_cell(lattice_constants)

    # Exercise functionality and check results
    standardized_lattice_constants, standardized_centering = standardize(
        lattice_constants, XtallographyUtils.PRIMITIVE
    )

    expected_lattice_constants = TriclinicLatticeConstants(b, c, a, π - β, π - γ, α)
    @test standardized_lattice_constants ≈ expected_lattice_constants

    @test standardized_centering == XtallographyUtils.PRIMITIVE

    # ------ a ≈ b ≈ c

    a = b = c = 1.0
    α = 1π / 10
    β = 8π / 10
    γ = 7π / 10
    lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)
    @test is_triclinic_type_I_cell(lattice_constants)

    # Exercise functionality and check results
    standardized_lattice_constants, standardized_centering = standardize(
        lattice_constants, XtallographyUtils.PRIMITIVE
    )

    expected_lattice_constants = TriclinicLatticeConstants(a, b, c, α, π - β, π - γ)
    @test standardized_lattice_constants ≈ expected_lattice_constants

    @test standardized_centering == XtallographyUtils.PRIMITIVE

    # ------ a ≈ b (after sorting a, b, c)

    a = b = 1.0
    c = 5
    α = 4π / 10
    β = 8π / 10
    γ = 7π / 10
    lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)
    @test is_triclinic_type_I_cell(lattice_constants)

    # Exercise functionality and check results
    standardized_lattice_constants, standardized_centering = standardize(
        lattice_constants, XtallographyUtils.PRIMITIVE
    )

    expected_lattice_constants = TriclinicLatticeConstants(a, b, c, π - β, α, π - γ)
    @test standardized_lattice_constants ≈ expected_lattice_constants

    @test standardized_centering == XtallographyUtils.PRIMITIVE

    # ------ b ≈ c (after sorting a, b, c)

    a = b = 5.0
    c = 1
    α = 8π / 10
    β = 7π / 10
    γ = 4π / 10
    lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)
    @test is_triclinic_type_I_cell(lattice_constants)

    # Exercise functionality and check results
    standardized_lattice_constants, standardized_centering = standardize(
        lattice_constants, XtallographyUtils.PRIMITIVE
    )

    expected_lattice_constants = TriclinicLatticeConstants(c, a, b, γ, π - α, π - β)
    @test standardized_lattice_constants ≈ expected_lattice_constants

    @test standardized_centering == XtallographyUtils.PRIMITIVE
end

@testset "standardize(): Type II cell" begin
    # --- Tests

    # ------ lattice constants already in standard form

    a = 1.0
    b = 5.0
    c = 10.0
    α = 6π / 10
    β = 7π / 10
    γ = 8π / 10
    lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)
    @test !is_triclinic_type_I_cell(lattice_constants)

    # Exercise functionality and check results
    standardized_lattice_constants, standardized_centering = standardize(
        lattice_constants, XtallographyUtils.PRIMITIVE
    )

    expected_lattice_constants = lattice_constants
    @test standardized_lattice_constants ≈ expected_lattice_constants

    @test standardized_centering == XtallographyUtils.PRIMITIVE

    # ------ edge lengths not sorted

    a = 1.0
    b = 10.0
    c = 5.0
    α = 6π / 10
    β = 7π / 10
    γ = 8π / 10
    lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)
    @test !is_triclinic_type_I_cell(lattice_constants)

    # Exercise functionality and check results
    standardized_lattice_constants, standardized_centering = standardize(
        lattice_constants, XtallographyUtils.PRIMITIVE
    )

    expected_lattice_constants = TriclinicLatticeConstants(a, c, b, α, γ, β)
    @test standardized_lattice_constants ≈ expected_lattice_constants

    @test standardized_centering == XtallographyUtils.PRIMITIVE

    # ------ origin not at "homogeneous corner"

    a = 1.0
    b = 5.0
    c = 10.0
    α = 4π / 10
    β = 3π / 10
    γ = 8π / 10
    lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)
    @test !is_triclinic_type_I_cell(lattice_constants)

    # Exercise functionality and check results
    standardized_lattice_constants, standardized_centering = standardize(
        lattice_constants, XtallographyUtils.PRIMITIVE
    )

    expected_lattice_constants = TriclinicLatticeConstants(a, b, c, π - α, π - β, γ)
    @test standardized_lattice_constants ≈ expected_lattice_constants

    @test standardized_centering == XtallographyUtils.PRIMITIVE

    # ------ edge lengths not sorted and origin not at "homogeneous corner"

    a = 10.0
    b = 1.0
    c = 5.0
    α = 4π / 10
    β = 3π / 10
    γ = 8π / 10
    lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)
    @test !is_triclinic_type_I_cell(lattice_constants)

    # Exercise functionality and check results
    standardized_lattice_constants, standardized_centering = standardize(
        lattice_constants, XtallographyUtils.PRIMITIVE
    )

    expected_lattice_constants = TriclinicLatticeConstants(b, c, a, π - β, γ, π - α)
    @test standardized_lattice_constants ≈ expected_lattice_constants

    @test standardized_centering == XtallographyUtils.PRIMITIVE

    # ------ a ≈ b ≈ c

    a = b = c = 1.0
    α = 4π / 10
    β = 7π / 10
    γ = 2π / 10
    lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)
    @test !is_triclinic_type_I_cell(lattice_constants)

    # Exercise functionality and check results
    standardized_lattice_constants, standardized_centering = standardize(
        lattice_constants, XtallographyUtils.PRIMITIVE
    )

    expected_lattice_constants = TriclinicLatticeConstants(a, b, c, π - α, β, π - γ)
    @test standardized_lattice_constants ≈ expected_lattice_constants

    @test standardized_centering == XtallographyUtils.PRIMITIVE

    # ------ a ≈ b (after sorting a, b, c)

    a = b = 1.0
    c = 5
    α = 7π / 10
    β = 4π / 10
    γ = 4π / 10
    lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)
    @test !is_triclinic_type_I_cell(lattice_constants)

    # Exercise functionality and check results
    standardized_lattice_constants, standardized_centering = standardize(
        lattice_constants, XtallographyUtils.PRIMITIVE
    )

    expected_lattice_constants = TriclinicLatticeConstants(a, b, c, π - β, α, π - γ)
    @test standardized_lattice_constants ≈ expected_lattice_constants

    @test standardized_centering == XtallographyUtils.PRIMITIVE

    # ------ b ≈ c (after sorting a, b, c)

    a = b = 5.0
    c = 1
    α = 8π / 10
    β = 4π / 10
    γ = 3π / 10
    lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)
    @test !is_triclinic_type_I_cell(lattice_constants)

    # Exercise functionality and check results
    standardized_lattice_constants, standardized_centering = standardize(
        lattice_constants, XtallographyUtils.PRIMITIVE
    )

    expected_lattice_constants = TriclinicLatticeConstants(c, a, b, π - γ, π - β, α)
    @test standardized_lattice_constants ≈ expected_lattice_constants

    @test standardized_centering == XtallographyUtils.PRIMITIVE
end

@testset "standardize(): invalid arguments" begin
    # --- Preparations

    lattice_constants = TriclinicLatticeConstants(1, 2, 3, 2π / 5, 3π / 5, 4π / 5)

    # --- Tests

    for centering in
        (XtallographyUtils.BODY, XtallographyUtils.FACE, XtallographyUtils.BASE)
        local error = nothing
        local error_message = ""
        try
            standardize(lattice_constants, centering)
        catch error
            bt = catch_backtrace()
            error_message = sprint(showerror, error, bt)
        end

        @test error isa ArgumentError

        expected_error =
            "ArgumentError: " *
            "Invalid Bravais lattice: (lattice_system=Triclinic, centering=$centering)"

        @test startswith(error_message, expected_error)
    end
end

# ------ iucr_conventional_cell()

@testset "iucr_conventional_cell(): limiting cases" begin
    # Note: this test set only tests basic functionality. Comprehensive coverage of
    #       different cases for triclinic bases is covered in the convert_to_mP()
    #       convert_to_mI(), and convert_to_mS() test sets.

    # --- Tests

    # ------ triclinic unit cell is not equivalent to a monoclinic unit cell

    # Construct lattice constants for triclinic unit cell
    a = 1.0
    b = 2.0
    c = 3.0
    α = π / 2
    β = 4π / 7
    γ = 3π / 5
    lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)

    # Exercise functionality
    iucr_unit_cell = iucr_conventional_cell(
        UnitCell(lattice_constants, XtallographyUtils.PRIMITIVE)
    )

    # Check results
    @test iucr_unit_cell.lattice_constants ≈ standardize(lattice_constants)
    @test iucr_unit_cell.centering == XtallographyUtils.PRIMITIVE

    # ------ triclinic unit cell is equivalent to a primitive monoclinic unit cell

    # Construct lattice constants for triclinic unit cell
    a = 1.0
    b = 2.0
    c = 3.0
    α = π / 4
    β = π / 2
    γ = π / 2
    lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)

    # Exercise functionality
    iucr_unit_cell = iucr_conventional_cell(
        UnitCell(lattice_constants, XtallographyUtils.PRIMITIVE)
    )

    # Check results
    @test iucr_unit_cell.lattice_constants ≈
        standardize(MonoclinicLatticeConstants(b, a, c, π - α))
    @test iucr_unit_cell.centering == XtallographyUtils.PRIMITIVE

    # ------ triclinic unit cell is equivalent to a body-centered monoclinic unit cell

    # Construct basis for monoclinic unit cell
    m_a = 1.0
    m_b = 2.0
    m_c = 3.0
    m_β = 3π / 5

    monoclinic_lattice_constants = MonoclinicLatticeConstants(m_a, m_b, m_c, m_β)

    m_basis_a, m_basis_b, m_basis_c = basis(monoclinic_lattice_constants)

    # Compute basis for triclinic unit cell
    basis_a = m_basis_a
    basis_b = m_basis_b
    basis_c = 0.5 * (basis_a + basis_b + m_basis_c)

    # Construct lattice constants for triclinic unit cell
    a = norm(basis_a)
    b = norm(basis_b)
    c = norm(basis_c)
    α = angle(basis_b, basis_c)
    β = angle(basis_c, basis_a)
    γ = angle(basis_a, basis_b)

    lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)

    # Exercise functionality
    iucr_unit_cell = iucr_conventional_cell(
        UnitCell(lattice_constants, XtallographyUtils.PRIMITIVE)
    )

    # Check results
    expected_lattice_constants, _ = standardize(
        MonoclinicLatticeConstants(m_a, m_b, m_c, m_β), XtallographyUtils.BODY
    )
    @test iucr_unit_cell.lattice_constants ≈ expected_lattice_constants
    @test iucr_unit_cell.centering == XtallographyUtils.BODY
end

@testset "iucr_conventional_cell(): limiting cases - base-centered" begin
    # Note: this test checks that all base-centered limiting cases are covered by
    #       convert_to_mI() and convert_to_mS().

    # --- Tests

    # ------ triclinic basis contains m_basis_b, m_basis_c, and one base-centered lattice
    #        vector

    # Construct basis for monoclinic unit cell
    m_a = 1.0
    m_b = 2.0
    m_c = 3.0
    m_β = 3π / 5

    monoclinic_lattice_constants = MonoclinicLatticeConstants(m_a, m_b, m_c, m_β)

    m_basis_a, m_basis_b, m_basis_c = basis(monoclinic_lattice_constants)

    # Compute basis for triclinic unit cell
    basis_a = 0.5 * (m_basis_a + m_basis_b)
    basis_b = m_basis_b
    basis_c = m_basis_c
    @test is_basis(basis_a, basis_b, basis_c)

    # Construct lattice constants for triclinic unit cell
    a = norm(basis_a)
    b = norm(basis_b)
    c = norm(basis_c)
    α = angle(basis_b, basis_c)
    β = angle(basis_c, basis_a)
    γ = angle(basis_a, basis_b)

    lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)

    # Exercise functionality
    iucr_unit_cell = iucr_conventional_cell(
        UnitCell(lattice_constants, XtallographyUtils.PRIMITIVE)
    )

    # Check results
    expected_lattice_constants, _ = standardize(
        MonoclinicLatticeConstants(m_a, m_b, m_c, m_β), XtallographyUtils.BASE
    )
    @test iucr_unit_cell.lattice_constants ≈ expected_lattice_constants
    @test iucr_unit_cell.centering == XtallographyUtils.BODY

    # ------ triclinic basis contains m_basis_a, m_basis_c, and one base-centered lattice
    #        vector. Case #1: m_basis_c makes no contribution to basis_b

    # Construct basis for monoclinic unit cell
    m_a = 1.0
    m_b = 2.0
    m_c = 3.0
    m_β = 3π / 5

    monoclinic_lattice_constants = MonoclinicLatticeConstants(m_a, m_b, m_c, m_β)

    m_basis_a, m_basis_b, m_basis_c = basis(monoclinic_lattice_constants)

    # Compute basis for triclinic unit cell
    basis_a = m_basis_a
    basis_b = 0.5 * (m_basis_a + m_basis_b)
    basis_c = m_basis_c
    @test is_basis(basis_a, basis_b, basis_c)

    # Construct lattice constants for triclinic unit cell
    a = norm(basis_a)
    b = norm(basis_b)
    c = norm(basis_c)
    α = angle(basis_b, basis_c)
    β = angle(basis_c, basis_a)
    γ = angle(basis_a, basis_b)

    lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)

    # Exercise functionality
    iucr_unit_cell = iucr_conventional_cell(
        UnitCell(lattice_constants, XtallographyUtils.PRIMITIVE)
    )

    # Check results
    expected_lattice_constants, _ = standardize(
        MonoclinicLatticeConstants(m_a, m_b, m_c, m_β), XtallographyUtils.BASE
    )
    @test iucr_unit_cell.lattice_constants ≈ expected_lattice_constants
    @test iucr_unit_cell.centering == XtallographyUtils.BODY

    # ------ triclinic basis contains m_basis_a, m_basis_c, and one base-centered lattice
    #        vector. Case #2: m_basis_c contributes to basis_b

    # Construct basis for monoclinic unit cell
    m_a = 1.0
    m_b = 2.0
    m_c = 3.0
    m_β = 3π / 5

    monoclinic_lattice_constants = MonoclinicLatticeConstants(m_a, m_b, m_c, m_β)

    m_basis_a, m_basis_b, m_basis_c = basis(monoclinic_lattice_constants)

    # Compute basis for triclinic unit cell
    basis_a = m_basis_a
    basis_b = 0.5 * (m_basis_a + m_basis_b) + m_basis_c
    basis_c = m_basis_c
    @test is_basis(basis_a, basis_b, basis_c)

    # Construct lattice constants for triclinic unit cell
    a = norm(basis_a)
    b = norm(basis_b)
    c = norm(basis_c)
    α = angle(basis_b, basis_c)
    β = angle(basis_c, basis_a)
    γ = angle(basis_a, basis_b)

    lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)

    # Exercise functionality
    iucr_unit_cell = iucr_conventional_cell(
        UnitCell(lattice_constants, XtallographyUtils.PRIMITIVE)
    )

    # Check results
    expected_lattice_constants, _ = standardize(
        MonoclinicLatticeConstants(m_a, m_b, m_c, m_β), XtallographyUtils.BASE
    )
    @test iucr_unit_cell.lattice_constants ≈ expected_lattice_constants
    @test iucr_unit_cell.centering == XtallographyUtils.BODY

    # ------ triclinic basis m_basis_a and two base-centered lattice vectors.
    #        Case #1: relative signs of m_basis_a and m_basis_b are the same in basis_b
    #        and basis_c

    # Construct basis for monoclinic unit cell
    m_a = 1.0
    m_b = 2.0
    m_c = 3.0
    m_β = 3π / 5

    monoclinic_lattice_constants = MonoclinicLatticeConstants(m_a, m_b, m_c, m_β)

    m_basis_a, m_basis_b, m_basis_c = basis(monoclinic_lattice_constants)

    # Compute basis for triclinic unit cell
    basis_a = basis_a
    basis_b = 0.5 * (m_basis_a - m_basis_b)
    basis_c = 0.5 * (m_basis_a - m_basis_b) - m_basis_c
    @test is_basis(basis_a, basis_b, basis_c)

    # Construct lattice constants for triclinic unit cell
    a = norm(basis_a)
    b = norm(basis_b)
    c = norm(basis_c)
    α = angle(basis_b, basis_c)
    β = angle(basis_c, basis_a)
    γ = angle(basis_a, basis_b)

    lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)

    # Exercise functionality
    iucr_unit_cell = iucr_conventional_cell(
        UnitCell(lattice_constants, XtallographyUtils.PRIMITIVE)
    )

    # Check results
    expected_lattice_constants, _ = standardize(
        MonoclinicLatticeConstants(m_a, m_b, m_c, m_β), XtallographyUtils.BASE
    )
    @test iucr_unit_cell.lattice_constants ≈ expected_lattice_constants
    @test iucr_unit_cell.centering == XtallographyUtils.BODY

    # ------ triclinic basis m_basis_a and two base-centered lattice vectors.
    #        Case #2: relative signs of m_basis_a and m_basis_b are opposite in basis_b
    #        and basis_c

    # Construct basis for monoclinic unit cell
    m_a = 1.0
    m_b = 2.0
    m_c = 3.0
    m_β = 3π / 5

    monoclinic_lattice_constants = MonoclinicLatticeConstants(m_a, m_b, m_c, m_β)

    m_basis_a, m_basis_b, m_basis_c = basis(monoclinic_lattice_constants)

    # Compute basis for triclinic unit cell
    basis_a = basis_a
    basis_b = 0.5 * (m_basis_a + m_basis_b)
    basis_c = 0.5 * (m_basis_a - m_basis_b) - m_basis_c
    @test is_basis(basis_a, basis_b, basis_c)

    # Construct lattice constants for triclinic unit cell
    a = norm(basis_a)
    b = norm(basis_b)
    c = norm(basis_c)
    α = angle(basis_b, basis_c)
    β = angle(basis_c, basis_a)
    γ = angle(basis_a, basis_b)

    lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)

    # Exercise functionality
    iucr_unit_cell = iucr_conventional_cell(
        UnitCell(lattice_constants, XtallographyUtils.PRIMITIVE)
    )

    # Check results
    expected_lattice_constants, _ = standardize(
        MonoclinicLatticeConstants(m_a, m_b, m_c, m_β), XtallographyUtils.BASE
    )
    @test iucr_unit_cell.lattice_constants ≈ expected_lattice_constants
    @test iucr_unit_cell.centering == XtallographyUtils.BODY

    # ------ triclinic basis m_basis_b and two base-centered lattice vectors

    # Construct basis for monoclinic unit cell
    m_a = 1.0
    m_b = 2.0
    m_c = 3.0
    m_β = 3π / 5

    monoclinic_lattice_constants = MonoclinicLatticeConstants(m_a, m_b, m_c, m_β)

    m_basis_a, m_basis_b, m_basis_c = basis(monoclinic_lattice_constants)

    # Compute basis for triclinic unit cell
    basis_a = 0.5 * (m_basis_a + m_basis_b)
    basis_b = m_basis_b
    basis_c = 0.5 * (m_basis_a + m_basis_b) + m_basis_c
    @test is_basis(basis_a, basis_b, basis_c)

    # Construct lattice constants for triclinic unit cell
    a = norm(basis_a)
    b = norm(basis_b)
    c = norm(basis_c)
    α = angle(basis_b, basis_c)
    β = angle(basis_c, basis_a)
    γ = angle(basis_a, basis_b)

    lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)

    # Exercise functionality
    iucr_unit_cell = iucr_conventional_cell(
        UnitCell(lattice_constants, XtallographyUtils.PRIMITIVE)
    )

    # Check results
    expected_lattice_constants, _ = standardize(
        MonoclinicLatticeConstants(m_a, m_b, m_c, m_β), XtallographyUtils.BASE
    )
    @test iucr_unit_cell.lattice_constants ≈ expected_lattice_constants
    @test iucr_unit_cell.centering == XtallographyUtils.BODY

    # ------ triclinic basis m_basis_c and two base-centered lattice vectors

    # Construct basis for monoclinic unit cell
    m_a = 1.0
    m_b = 2.0
    m_c = 3.0
    m_β = 3π / 5

    monoclinic_lattice_constants = MonoclinicLatticeConstants(m_a, m_b, m_c, m_β)

    m_basis_a, m_basis_b, m_basis_c = basis(monoclinic_lattice_constants)

    # Compute basis for triclinic unit cell
    basis_a = 0.5 * (m_basis_a + m_basis_b)
    basis_b = 0.5 * (m_basis_a - m_basis_b) + m_basis_c
    basis_c = m_basis_c
    @test is_basis(basis_a, basis_b, basis_c)

    # Construct lattice constants for triclinic unit cell
    a = norm(basis_a)
    b = norm(basis_b)
    c = norm(basis_c)
    α = angle(basis_b, basis_c)
    β = angle(basis_c, basis_a)
    γ = angle(basis_a, basis_b)

    lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)

    # Exercise functionality
    iucr_unit_cell = iucr_conventional_cell(
        UnitCell(lattice_constants, XtallographyUtils.PRIMITIVE)
    )

    # Check results
    expected_lattice_constants, _ = standardize(
        MonoclinicLatticeConstants(m_a, m_b, m_c, m_β), XtallographyUtils.BASE
    )
    @test iucr_unit_cell.lattice_constants ≈ expected_lattice_constants
    @test iucr_unit_cell.centering == XtallographyUtils.BODY

    # ------ triclinic basis three base-centered lattice vectors

    # Construct basis for monoclinic unit cell
    m_a = 1.0
    m_b = 2.0
    m_c = 3.0
    m_β = 3π / 5

    monoclinic_lattice_constants = MonoclinicLatticeConstants(m_a, m_b, m_c, m_β)

    m_basis_a, m_basis_b, m_basis_c = basis(monoclinic_lattice_constants)

    # Compute basis for triclinic unit cell
    basis_a = 0.5 * (m_basis_a + m_basis_b)
    basis_b = 0.5 * (m_basis_a - m_basis_b)
    basis_c = 0.5 * (m_basis_a - m_basis_b) + m_basis_c
    @test is_basis(basis_a, basis_b, basis_c)

    # Construct lattice constants for triclinic unit cell
    a = norm(basis_a)
    b = norm(basis_b)
    c = norm(basis_c)
    α = angle(basis_b, basis_c)
    β = angle(basis_c, basis_a)
    γ = angle(basis_a, basis_b)

    lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)

    # Exercise functionality
    iucr_unit_cell = iucr_conventional_cell(
        UnitCell(lattice_constants, XtallographyUtils.PRIMITIVE)
    )

    # Check results
    expected_lattice_constants, _ = standardize(
        MonoclinicLatticeConstants(m_a, m_b, m_c, m_β), XtallographyUtils.BASE
    )
    @test iucr_unit_cell.lattice_constants ≈ expected_lattice_constants
    @test iucr_unit_cell.centering == XtallographyUtils.BODY
end

@testset "iucr_conventional_cell(): invalid arguments" begin
    # --- Preparations

    # Construct lattice constants for triclinic unit cell
    a = 1.0
    b = 2.0
    c = 3.0
    α = 2π / 5
    β = 3π / 5
    γ = 4π / 5
    lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)

    # --- Tests

    for centering in
        (XtallographyUtils.BODY, XtallographyUtils.FACE, XtallographyUtils.BASE)
        local error = nothing
        local error_message = ""
        try
            iucr_conventional_cell(UnitCell(lattice_constants, centering))
        catch error
            bt = catch_backtrace()
            error_message = sprint(showerror, error, bt)
        end

        @test error isa ArgumentError

        expected_error =
            "ArgumentError: " *
            "Invalid Bravais lattice: (lattice_system=Triclinic, centering=$centering)"

        @test startswith(error_message, expected_error)
    end
end

# ------ convert_to_mP()

@testset "convert_to_mP()" begin
    # --- Tests

    # β ≈ π / 2 && γ ≈ π / 2
    a = 1.0
    b = 2.0
    c = 3.0
    α = π / 4
    β = π / 2
    γ = π / 2
    triclinic_lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)

    monoclinic_lattice_constants = convert_to_mP(triclinic_lattice_constants)

    @test monoclinic_lattice_constants ≈
        standardize(MonoclinicLatticeConstants(b, a, c, π - α))

    # α ≈ π / 2 && γ ≈ π / 2
    a = 3.0
    b = 2.0
    c = 1.0
    α = π / 2
    β = π / 3
    γ = π / 2
    triclinic_lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)

    monoclinic_lattice_constants = convert_to_mP(triclinic_lattice_constants)

    @test monoclinic_lattice_constants ≈
        standardize(MonoclinicLatticeConstants(c, b, a, π - β))

    # α ≈ π / 2 && β ≈ π / 2
    a = 1.0
    b = 2.0
    c = 3.0
    α = π / 2
    β = π / 2
    γ = 3π / 5
    triclinic_lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)

    monoclinic_lattice_constants = convert_to_mP(triclinic_lattice_constants)

    @test monoclinic_lattice_constants ≈ standardize(MonoclinicLatticeConstants(a, c, b, γ))

    # triclinic unit cell is not equivalent to a primitive monoclinic unit cell
    a = 1.0
    b = 2.0
    c = 3.0
    α = π / 2
    β = 4π / 7
    γ = 3π / 5
    triclinic_lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)

    local error = nothing
    local error_message = ""
    try
        convert_to_mP(triclinic_lattice_constants)

    catch error
        bt = catch_backtrace()
        error_message = sprint(showerror, error, bt)
    end

    @test error isa ErrorException

    expected_error =
        "The triclinic unit cell defined by `lattice_constants` is not equivalent to " *
        "a primitive monoclinic unit cell."
    @test startswith(error_message, expected_error)
end

# ------ convert_to_mI()

# Helper functions
function convert_to_mI_test_all_triclinic_basis_permutations(
    basis_a::Vector{<:Real},
    basis_b::Vector{<:Real},
    basis_c::Vector{<:Real},
    expected_monoclinic_lattice_constants::MonoclinicLatticeConstants,
)
    # --- Preparations

    # Construct expected unit cell
    expected_monoclinic_unit_cell = UnitCell(
        expected_monoclinic_lattice_constants, XtallographyUtils.BODY
    )

    # --- Test convert_to_mI() for each permutation of the basis vectors

    for permuted_basis in permutations([basis_a, basis_b, basis_c])
        # Extract permuted basis vectors
        permuted_basis_a = permuted_basis[1]
        permuted_basis_b = permuted_basis[2]
        permuted_basis_c = permuted_basis[3]

        # Construct lattice constants for triclinic unit cell
        a = norm(permuted_basis_a)
        b = norm(permuted_basis_b)
        c = norm(permuted_basis_c)
        α = angle(permuted_basis_b, permuted_basis_c)
        β = angle(permuted_basis_c, permuted_basis_a)
        γ = angle(permuted_basis_a, permuted_basis_b)

        lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)

        # Exercise functionality
        monoclinic_lattice_constants = convert_to_mI(lattice_constants)
        monoclinic_unit_cell = UnitCell(
            monoclinic_lattice_constants, XtallographyUtils.BODY
        )

        # Check results
        @test is_equivalent_unit_cell(monoclinic_unit_cell, expected_monoclinic_unit_cell)
    end
end

@testset "convert_to_mI(): {m_a, m_b} in aP basis, (m_c/m_a) |cos m_β| < 1" begin
    # --- Preparations

    # Construct basis for monoclinic unit cell
    m_a = 1.0
    m_b = 2.0
    m_c = 3.0
    m_β = 3π / 5

    # Check test conditions
    @test abs(m_c / m_a * cos(m_β)) < 1

    expected_monoclinic_lattice_constants = MonoclinicLatticeConstants(m_a, m_b, m_c, m_β)

    m_basis_a, m_basis_b, m_basis_c = basis(expected_monoclinic_lattice_constants)

    # --- Tests

    # Test all cases for body-centered basis vector
    basis_a = m_basis_a
    basis_b = m_basis_b
    for signs in Iterators.product((-1, 1), (-1, 1), (-1, 1))
        basis_c = 0.5 * (signs[1] * basis_a + signs[2] * basis_b + signs[3] * m_basis_c)

        convert_to_mI_test_all_triclinic_basis_permutations(
            basis_a, basis_b, basis_c, expected_monoclinic_lattice_constants
        )
    end
end

@testset "convert_to_mI(): {m_a, m_b} in aP basis, (m_c/m_a) |cos m_β| > 1" begin
    # --- Preparations

    # Construct basis for monoclinic unit cell
    m_a = 1.0
    m_b = 2.0
    m_c = 10.0
    m_β = 3π / 5

    # Check test conditions
    @test abs(m_c / m_a * cos(m_β)) > 1

    expected_monoclinic_lattice_constants = MonoclinicLatticeConstants(m_a, m_b, m_c, m_β)

    m_basis_a, m_basis_b, m_basis_c = basis(expected_monoclinic_lattice_constants)

    # --- Tests

    # Test all cases for body-centered basis vector
    basis_a = m_basis_a
    basis_b = m_basis_b
    for signs in Iterators.product((-1, 1), (-1, 1), (-1, 1))
        basis_c = 0.5 * (signs[1] * basis_a + signs[2] * basis_b + signs[3] * m_basis_c)

        convert_to_mI_test_all_triclinic_basis_permutations(
            basis_a, basis_b, basis_c, expected_monoclinic_lattice_constants
        )
    end
end

@testset "convert_to_mI(): {m_b, m_c} in aP basis, (m_c/m_a) |cos m_β| < 1" begin
    # --- Preparations

    # Construct basis for monoclinic unit cell
    m_a = 1.0
    m_b = 2.0
    m_c = 3.0
    m_β = 3π / 5

    # Check test conditions
    @test abs(m_c / m_a * cos(m_β)) < 1

    expected_monoclinic_lattice_constants = MonoclinicLatticeConstants(m_a, m_b, m_c, m_β)

    m_basis_a, m_basis_b, m_basis_c = basis(expected_monoclinic_lattice_constants)

    # --- Tests

    # Test all cases for body-centered basis vector
    basis_a = m_basis_b
    basis_b = m_basis_c
    for signs in Iterators.product((-1, 1), (-1, 1), (-1, 1))
        basis_c = 0.5 * (signs[1] * basis_a + signs[2] * basis_b + signs[3] * m_basis_a)

        convert_to_mI_test_all_triclinic_basis_permutations(
            basis_a, basis_b, basis_c, expected_monoclinic_lattice_constants
        )
    end
end

@testset "convert_to_mI(): {m_b, m_c} in aP basis, (m_c/m_a) |cos m_β| > 1" begin
    # --- Preparations

    # Construct basis for monoclinic unit cell
    m_a = 1.0
    m_b = 2.0
    m_c = 10.0
    m_β = 3π / 5

    # Check test conditions
    @test abs(m_c / m_a * cos(m_β)) > 1

    expected_monoclinic_lattice_constants = MonoclinicLatticeConstants(m_a, m_b, m_c, m_β)

    m_basis_a, m_basis_b, m_basis_c = basis(expected_monoclinic_lattice_constants)

    # --- Tests

    # Test all cases for body-centered basis vector
    basis_a = m_basis_b
    basis_b = m_basis_c
    for signs in Iterators.product((-1, 1), (-1, 1), (-1, 1))
        basis_c = 0.5 * (signs[1] * basis_a + signs[2] * basis_b + signs[3] * m_basis_a)

        convert_to_mI_test_all_triclinic_basis_permutations(
            basis_a, basis_b, basis_c, expected_monoclinic_lattice_constants
        )
    end
end

@testset "convert_to_mI(): {m_a, m_c} in aP basis, (m_c/m_a) |cos m_β| < 1" begin
    # --- Preparations

    # Construct basis for monoclinic unit cell
    m_a = 1.0
    m_b = 2.0
    m_c = 3.0
    m_β = 3π / 5

    # Check test conditions
    @test abs(m_c / m_a * cos(m_β)) < 1

    expected_monoclinic_lattice_constants = MonoclinicLatticeConstants(m_a, m_b, m_c, m_β)

    m_basis_a, m_basis_b, m_basis_c = basis(expected_monoclinic_lattice_constants)

    # --- Tests

    # Test all cases for body-centered basis vector
    basis_a = m_basis_a
    basis_c = m_basis_c
    for signs in Iterators.product((-1, 1), (-1, 1), (-1, 1))
        basis_b = 0.5 * (signs[1] * basis_a + signs[2] * m_basis_b + signs[3] * basis_c)

        convert_to_mI_test_all_triclinic_basis_permutations(
            basis_a, basis_b, basis_c, expected_monoclinic_lattice_constants
        )
    end
end

@testset "convert_to_mI(): {m_a, m_c} in aP basis, (m_c/m_a) |cos m_β| > 1" begin
    # --- Preparations

    # Construct basis for monoclinic unit cell
    m_a = 1.0
    m_b = 2.0
    m_c = 10.0
    m_β = 3π / 5

    # Check test conditions
    @test abs(m_c / m_a * cos(m_β)) > 1

    expected_monoclinic_lattice_constants = MonoclinicLatticeConstants(m_a, m_b, m_c, m_β)

    m_basis_a, m_basis_b, m_basis_c = basis(expected_monoclinic_lattice_constants)

    # --- Tests

    # Test all cases for body-centered basis vector
    basis_a = m_basis_a
    basis_c = m_basis_c
    for signs in Iterators.product((-1, 1), (-1, 1), (-1, 1))
        basis_b = 0.5 * (signs[1] * basis_a + signs[2] * m_basis_b + signs[3] * basis_c)

        convert_to_mI_test_all_triclinic_basis_permutations(
            basis_a, basis_b, basis_c, expected_monoclinic_lattice_constants
        )
    end
end

@testset "convert_to_mI(): only m_a in aP basis, (m_c/m_a) |cos m_β| < 1" begin
    # --- Preparations

    # Construct basis for monoclinic unit cell
    m_a = 1.0
    m_b = 2.0
    m_c = 3.0
    m_β = 3π / 5

    # Check test conditions
    @test abs(m_c / m_a * cos(m_β)) < 1

    expected_monoclinic_lattice_constants = MonoclinicLatticeConstants(m_a, m_b, m_c, m_β)

    m_basis_a, m_basis_b, m_basis_c = basis(expected_monoclinic_lattice_constants)

    # --- Tests

    # Test all cases for body-centered basis vector
    basis_a = m_basis_a
    signs = collect(Iterators.product((-1, 1), (-1, 1), (-1, 1)))
    for (signs_b, signs_c) in combinations(signs, 2)
        basis_b =
            0.5 * (signs_b[1] * basis_a + signs_b[2] * m_basis_b + signs_b[3] * m_basis_c)
        basis_c =
            0.5 * (signs_c[1] * basis_a + signs_c[2] * m_basis_b + signs_c[3] * m_basis_c)

        # Check that basis vectors are linearly independent
        if !is_basis(basis_a, basis_b, basis_c)
            continue
        end

        convert_to_mI_test_all_triclinic_basis_permutations(
            basis_a, basis_b, basis_c, expected_monoclinic_lattice_constants
        )
    end
end

@testset "convert_to_mI(): only m_a in aP basis, (m_c/m_a) |cos m_β| > 1" begin
    # --- Preparations

    # Construct basis for monoclinic unit cell
    m_a = 1.0
    m_b = 2.0
    m_c = 10.0
    m_β = 3π / 5

    # Check test conditions
    @test abs(m_c / m_a * cos(m_β)) > 1

    expected_monoclinic_lattice_constants = MonoclinicLatticeConstants(m_a, m_b, m_c, m_β)

    m_basis_a, m_basis_b, m_basis_c = basis(expected_monoclinic_lattice_constants)

    # --- Tests

    # Test all cases for body-centered basis vector
    basis_a = m_basis_a
    signs = collect(Iterators.product((-1, 1), (-1, 1), (-1, 1)))
    for (signs_b, signs_c) in combinations(signs, 2)
        basis_b =
            0.5 * (signs_b[1] * basis_a + signs_b[2] * m_basis_b + signs_b[3] * m_basis_c)
        basis_c =
            0.5 * (signs_c[1] * basis_a + signs_c[2] * m_basis_b + signs_c[3] * m_basis_c)

        # Check that basis vectors are linearly independent
        if !is_basis(basis_a, basis_b, basis_c)
            continue
        end

        convert_to_mI_test_all_triclinic_basis_permutations(
            basis_a, basis_b, basis_c, expected_monoclinic_lattice_constants
        )
    end
end

@testset "convert_to_mI(): only m_b in aP basis, (m_c/m_a) |cos m_β| < 1" begin
    # --- Preparations

    # Construct basis for monoclinic unit cell
    m_a = 1.0
    m_b = 2.0
    m_c = 3.0
    m_β = 3π / 5

    # Check test conditions
    @test abs(m_c / m_a * cos(m_β)) < 1

    expected_monoclinic_lattice_constants = MonoclinicLatticeConstants(m_a, m_b, m_c, m_β)

    m_basis_a, m_basis_b, m_basis_c = basis(expected_monoclinic_lattice_constants)

    # --- Tests

    # Test all cases for body-centered basis vector
    basis_b = m_basis_b
    signs = collect(Iterators.product((-1, 1), (-1, 1), (-1, 1)))
    for (signs_a, signs_c) in combinations(signs, 2)
        basis_a =
            0.5 * (signs_a[1] * m_basis_a + signs_a[2] * basis_b + signs_a[3] * m_basis_c)
        basis_c =
            0.5 * (signs_c[1] * m_basis_a + signs_c[2] * basis_b + signs_c[3] * m_basis_c)

        # Check that basis vectors are linearly independent
        if !is_basis(basis_a, basis_b, basis_c)
            continue
        end

        convert_to_mI_test_all_triclinic_basis_permutations(
            basis_a, basis_b, basis_c, expected_monoclinic_lattice_constants
        )
    end
end

@testset "convert_to_mI(): only m_b in aP basis, (m_c/m_a) |cos m_β| > 1" begin
    # --- Preparations

    # Construct basis for monoclinic unit cell
    m_a = 1.0
    m_b = 2.0
    m_c = 10.0
    m_β = 3π / 5

    # Check test conditions
    @test abs(m_c / m_a * cos(m_β)) > 1

    expected_monoclinic_lattice_constants = MonoclinicLatticeConstants(m_a, m_b, m_c, m_β)

    m_basis_a, m_basis_b, m_basis_c = basis(expected_monoclinic_lattice_constants)

    # --- Tests

    # Test all cases for body-centered basis vector
    basis_b = m_basis_b
    signs = collect(Iterators.product((-1, 1), (-1, 1), (-1, 1)))
    for (signs_a, signs_c) in combinations(signs, 2)
        basis_a =
            0.5 * (signs_a[1] * m_basis_a + signs_a[2] * basis_b + signs_a[3] * m_basis_c)
        basis_c =
            0.5 * (signs_c[1] * m_basis_a + signs_c[2] * basis_b + signs_c[3] * m_basis_c)

        # Check that basis vectors are linearly independent
        if !is_basis(basis_a, basis_b, basis_c)
            continue
        end

        convert_to_mI_test_all_triclinic_basis_permutations(
            basis_a, basis_b, basis_c, expected_monoclinic_lattice_constants
        )
    end
end

@testset "convert_to_mI(): only m_c in aP basis, (m_c/m_a) |cos m_β| < 1" begin
    # --- Preparations

    # Construct basis for monoclinic unit cell
    m_a = 1.0
    m_b = 2.0
    m_c = 3.0
    m_β = 3π / 5

    # Check test conditions
    @test abs(m_c / m_a * cos(m_β)) < 1

    expected_monoclinic_lattice_constants = MonoclinicLatticeConstants(m_a, m_b, m_c, m_β)

    m_basis_a, m_basis_b, m_basis_c = basis(expected_monoclinic_lattice_constants)

    # --- Tests

    # Test all cases for body-centered basis vector
    basis_c = m_basis_c
    signs = collect(Iterators.product((-1, 1), (-1, 1), (-1, 1)))
    for (signs_a, signs_b) in combinations(signs, 2)
        basis_a =
            0.5 * (signs_a[1] * m_basis_a + signs_a[2] * m_basis_b + signs_a[3] * basis_c)
        basis_b =
            0.5 * (signs_b[1] * m_basis_a + signs_b[2] * m_basis_b + signs_b[3] * basis_c)

        # Check that basis vectors are linearly independent
        if !is_basis(basis_a, basis_b, basis_c)
            continue
        end

        convert_to_mI_test_all_triclinic_basis_permutations(
            basis_a, basis_b, basis_c, expected_monoclinic_lattice_constants
        )
    end
end

@testset "convert_to_mI(): only m_c in aP basis, (m_c/m_a) |cos m_β| > 1" begin
    # --- Preparations

    # Construct basis for monoclinic unit cell
    m_a = 1.0
    m_b = 2.0
    m_c = 10.0
    m_β = 3π / 5

    # Check test conditions
    @test abs(m_c / m_a * cos(m_β)) > 1

    expected_monoclinic_lattice_constants = MonoclinicLatticeConstants(m_a, m_b, m_c, m_β)

    m_basis_a, m_basis_b, m_basis_c = basis(expected_monoclinic_lattice_constants)

    # --- Tests

    # Test all cases for body-centered basis vector
    basis_c = m_basis_c
    signs = collect(Iterators.product((-1, 1), (-1, 1), (-1, 1)))
    for (signs_a, signs_b) in combinations(signs, 2)
        basis_a =
            0.5 * (signs_a[1] * m_basis_a + signs_a[2] * m_basis_b + signs_a[3] * basis_c)
        basis_b =
            0.5 * (signs_b[1] * m_basis_a + signs_b[2] * m_basis_b + signs_b[3] * basis_c)

        # Check that basis vectors are linearly independent
        if !is_basis(basis_a, basis_b, basis_c)
            continue
        end

        convert_to_mI_test_all_triclinic_basis_permutations(
            basis_a, basis_b, basis_c, expected_monoclinic_lattice_constants
        )
    end
end

@testset "convert_to_mI(): no mI basis vectors in aP basis, (m_c/m_a) |cos m_β| < 1" begin
    # --- Preparations

    # Construct basis for monoclinic unit cell
    m_a = 1.0
    m_b = 2.0
    m_c = 3.0
    m_β = 3π / 5

    # Check test conditions
    @test abs(m_c / m_a * cos(m_β)) < 1

    expected_monoclinic_lattice_constants = MonoclinicLatticeConstants(m_a, m_b, m_c, m_β)

    m_basis_a, m_basis_b, m_basis_c = basis(expected_monoclinic_lattice_constants)

    # --- Tests

    # Test all cases for body-centered basis vector
    signs = collect(Iterators.product((-1, 1), (-1, 1), (-1, 1)))
    for (signs_a, signs_b, signs_c) in combinations(signs, 3)
        basis_a =
            0.5 * (signs_a[1] * m_basis_a + signs_a[2] * m_basis_b + signs_a[3] * m_basis_c)
        basis_b =
            0.5 * (signs_b[1] * m_basis_a + signs_b[2] * m_basis_b + signs_b[3] * m_basis_c)
        basis_c =
            0.5 * (signs_c[1] * m_basis_a + signs_c[2] * m_basis_b + signs_c[3] * m_basis_c)

        # Check that basis vectors are linearly independent
        if !is_basis(basis_a, basis_b, basis_c)
            continue
        end

        convert_to_mI_test_all_triclinic_basis_permutations(
            basis_a, basis_b, basis_c, expected_monoclinic_lattice_constants
        )
    end
end

@testset "convert_to_mI(): no mI basis vectors in aP basis, (m_c/m_a) |cos m_β| > 1" begin
    # --- Preparations

    # Construct basis for monoclinic unit cell
    m_a = 1.0
    m_b = 2.0
    m_c = 10.0
    m_β = 3π / 5

    # Check test conditions
    @test abs(m_c / m_a * cos(m_β)) > 1

    expected_monoclinic_lattice_constants = MonoclinicLatticeConstants(m_a, m_b, m_c, m_β)

    m_basis_a, m_basis_b, m_basis_c = basis(expected_monoclinic_lattice_constants)

    # --- Tests

    # Test all cases for body-centered basis vector
    signs = collect(Iterators.product((-1, 1), (-1, 1), (-1, 1)))
    for (signs_a, signs_b, signs_c) in combinations(signs, 3)
        basis_a =
            0.5 * (signs_a[1] * m_basis_a + signs_a[2] * m_basis_b + signs_a[3] * m_basis_c)
        basis_b =
            0.5 * (signs_b[1] * m_basis_a + signs_b[2] * m_basis_b + signs_b[3] * m_basis_c)
        basis_c =
            0.5 * (signs_c[1] * m_basis_a + signs_c[2] * m_basis_b + signs_c[3] * m_basis_c)

        # Check that basis vectors are linearly independent
        if !is_basis(basis_a, basis_b, basis_c)
            continue
        end

        convert_to_mI_test_all_triclinic_basis_permutations(
            basis_a, basis_b, basis_c, expected_monoclinic_lattice_constants
        )
    end
end

@testset "convert_to_mI(): aP cell not equivalent to a mI cell" begin
    # --- Tests

    # triclinic unit cell is not equivalent to a body-centered monoclinic unit cell
    a = 1.0
    b = 2.0
    c = 3.0
    α = π / 2
    β = 4π / 7
    γ = 3π / 5
    lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)

    local error = nothing
    local error_message = ""
    try
        convert_to_mI(lattice_constants)

    catch error
        bt = catch_backtrace()
        error_message = sprint(showerror, error, bt)
    end

    @test error isa ErrorException

    expected_error =
        "The triclinic unit cell defined by `lattice_constants` is not equivalent to " *
        "a body-centered monoclinic unit cell."
    @test startswith(error_message, expected_error)
end

# ------ convert_to_mS()

# Helper functions
function convert_to_mS_test_all_triclinic_basis_permutations(
    basis_a::Vector{<:Real},
    basis_b::Vector{<:Real},
    basis_c::Vector{<:Real},
    expected_monoclinic_lattice_constants::MonoclinicLatticeConstants,
)
    # --- Preparations

    # Construct expected unit cell
    expected_monoclinic_unit_cell = UnitCell(
        expected_monoclinic_lattice_constants, XtallographyUtils.BASE
    )

    # --- Test convert_to_mS() for each permutation of the basis vectors

    # for permuted_basis in permutations([basis_a, basis_b, basis_c])
    for permuted_basis in ([basis_a, basis_b, basis_c],)
        # Extract permuted basis vectors
        permuted_basis_a = permuted_basis[1]
        permuted_basis_b = permuted_basis[2]
        permuted_basis_c = permuted_basis[3]

        # Construct lattice constants for triclinic unit cell
        a = norm(permuted_basis_a)
        b = norm(permuted_basis_b)
        c = norm(permuted_basis_c)
        α = angle(permuted_basis_b, permuted_basis_c)
        β = angle(permuted_basis_c, permuted_basis_a)
        γ = angle(permuted_basis_a, permuted_basis_b)

        lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)

        # Exercise functionality
        monoclinic_lattice_constants = convert_to_mS(lattice_constants)
        monoclinic_unit_cell = UnitCell(
            monoclinic_lattice_constants, XtallographyUtils.BASE
        )

        # Check results
        @test is_equivalent_unit_cell(monoclinic_unit_cell, expected_monoclinic_unit_cell)
    end
end

@testset "convert_to_mS_case_1a(): {m_a, m_c} in aP basis, |m_a| < |m_c|" begin
    # --- Preparations

    # Construct basis for monoclinic unit cell
    m_a = 1.0
    m_b = 2.0
    m_c = 3.0
    m_β = 3π / 5

    expected_monoclinic_lattice_constants = MonoclinicLatticeConstants(m_a, m_b, m_c, m_β)

    m_basis_a, m_basis_b, m_basis_c = basis(expected_monoclinic_lattice_constants)

    # --- Tests

    # Test all cases for base-centered basis vector
    basis_a = m_basis_a
    basis_c = m_basis_c
    for signs in Iterators.product((-1, 1), (-1, 1))
        basis_b = 0.5 * (signs[1] * m_basis_a + signs[2] * m_basis_b)

        convert_to_mS_test_all_triclinic_basis_permutations(
            basis_a, basis_b, basis_c, expected_monoclinic_lattice_constants
        )
    end
end

@testset "convert_to_mS_case_1a(): {m_a, m_c} in aP basis, |m_a| > |m_c|" begin
    # --- Preparations

    # Construct basis for monoclinic unit cell
    m_a = 3.0
    m_b = 2.0
    m_c = 1.0
    m_β = 3π / 5

    expected_monoclinic_lattice_constants = MonoclinicLatticeConstants(m_a, m_b, m_c, m_β)

    m_basis_a, m_basis_b, m_basis_c = basis(expected_monoclinic_lattice_constants)

    # --- Tests

    # Test all cases for base-centered basis vector
    basis_a = m_basis_a
    basis_c = m_basis_c
    for signs in Iterators.product((-1, 1), (-1, 1))
        basis_b = 0.5 * (signs[1] * m_basis_a + signs[2] * m_basis_b)

        convert_to_mS_test_all_triclinic_basis_permutations(
            basis_a, basis_b, basis_c, expected_monoclinic_lattice_constants
        )
    end
end

@testset "convert_to_mS_case_1b(): {m_a, m_c} in aP basis, |m_a| < |m_c|" begin
    # --- Preparations

    # Construct basis for monoclinic unit cell
    m_a = 1.0
    m_b = 2.0
    m_c = 3.0
    m_β = 3π / 5

    expected_monoclinic_lattice_constants = MonoclinicLatticeConstants(m_a, m_b, m_c, m_β)

    m_basis_a, m_basis_b, m_basis_c = basis(expected_monoclinic_lattice_constants)

    # --- Tests

    # Test all cases for base-centered basis vector
    basis_a = m_basis_a
    basis_c = m_basis_c
    for signs in Iterators.product((-1, 1), (-1, 1), (-1, 1))
        basis_b = 0.5 * (signs[1] * m_basis_a + signs[2] * m_basis_b) + signs[3] * m_basis_c

        convert_to_mS_test_all_triclinic_basis_permutations(
            basis_a, basis_b, basis_c, expected_monoclinic_lattice_constants
        )
    end
end

@testset "convert_to_mS_case_1b(): {m_a, m_c} in aP basis, |m_a| > |m_c|" begin
    # --- Preparations

    # Construct basis for monoclinic unit cell
    m_a = 3.0
    m_b = 2.0
    m_c = 1.0
    m_β = 3π / 5

    expected_monoclinic_lattice_constants = MonoclinicLatticeConstants(m_a, m_b, m_c, m_β)

    m_basis_a, m_basis_b, m_basis_c = basis(expected_monoclinic_lattice_constants)

    # --- Tests

    # Test all cases for base-centered basis vector
    basis_a = m_basis_a
    basis_c = m_basis_c
    for signs in Iterators.product((-1, 1), (-1, 1), (-1, 1))
        basis_b = 0.5 * (signs[1] * m_basis_a + signs[2] * m_basis_b) + signs[3] * m_basis_c

        convert_to_mS_test_all_triclinic_basis_permutations(
            basis_a, basis_b, basis_c, expected_monoclinic_lattice_constants
        )
    end
end

@testset "convert_to_mS_case_2a(): only m_a in aP basis, |m_a| < |m_c|" begin
    println("convert_to_mS_case_2a()")
    # --- Preparations

    # Construct basis for monoclinic unit cell
    m_a = 1.0
    m_b = 2.0
    m_c = 3.0
    m_β = 3π / 5

    expected_monoclinic_lattice_constants, _ = standardize(
        MonoclinicLatticeConstants(m_a, m_b, m_c, m_β), XtallographyUtils.BASE
    )

    m_basis_a, m_basis_b, m_basis_c = basis(expected_monoclinic_lattice_constants)

    # --- Tests

    # Test all cases for base-centered basis vector
    basis_a = m_basis_a
    for signs in Iterators.product((-1, 1), (-1, 1))
        basis_b = 0.5 * (signs[1] * m_basis_a + signs[1] * m_basis_b)
        basis_c = 0.5 * (signs[1] * m_basis_a + signs[1] * m_basis_b) + signs[2] * m_basis_c

        # Check that basis vectors are linearly independent
        if !is_basis(basis_a, basis_b, basis_c)
            continue
        end

        convert_to_mS_test_all_triclinic_basis_permutations(
            basis_a, basis_b, basis_c, expected_monoclinic_lattice_constants
        )
    end
end

@testset "convert_to_mS_case_2a(): only m_a in aP basis, |m_a| > |m_c|" begin
    println("convert_to_mS_case_2a()")
    # --- Preparations

    # Construct basis for monoclinic unit cell
    m_a = 3.0
    m_b = 2.0
    m_c = 1.0
    m_β = 3π / 5

    expected_monoclinic_lattice_constants = MonoclinicLatticeConstants(m_a, m_b, m_c, m_β)

    m_basis_a, m_basis_b, m_basis_c = basis(expected_monoclinic_lattice_constants)

    # --- Tests

    # Test all cases for base-centered basis vector
    basis_a = m_basis_a
    for signs in Iterators.product((-1, 1), (-1, 1))
        basis_b = 0.5 * (signs[1] * m_basis_a + signs[1] * m_basis_b)
        basis_c = 0.5 * (signs[1] * m_basis_a + signs[1] * m_basis_b) + signs[2] * m_basis_c

        # Check that basis vectors are linearly independent
        if !is_basis(basis_a, basis_b, basis_c)
            continue
        end

        convert_to_mS_test_all_triclinic_basis_permutations(
            basis_a, basis_b, basis_c, expected_monoclinic_lattice_constants
        )
    end
end

@testset "convert_to_mS_case_2b(): only m_a in aP basis, |m_a| < |m_c|" begin
    println("convert_to_mS_case_2b()")
    # --- Preparations

    # Construct basis for monoclinic unit cell
    m_a = 1.0
    m_b = 2.0
    m_c = 3.0
    m_β = 3π / 5

    expected_monoclinic_lattice_constants, _ = standardize(
        MonoclinicLatticeConstants(m_a, m_b, m_c, m_β), XtallographyUtils.BASE
    )

    m_basis_a, m_basis_b, m_basis_c = basis(expected_monoclinic_lattice_constants)

    # --- Tests

    # Test all cases for base-centered basis vector
    basis_a = m_basis_a
    for signs in Iterators.product((-1, 1), (-1, 1))
        basis_b = 0.5 * (signs[1] * m_basis_a - signs[2] * m_basis_b)
        basis_c = 0.5 * (signs[1] * m_basis_a + signs[2] * m_basis_b) + signs[1] * m_basis_c

        # Check that basis vectors are linearly independent
        if !is_basis(basis_a, basis_b, basis_c)
            continue
        end

        convert_to_mS_test_all_triclinic_basis_permutations(
            basis_a, basis_b, basis_c, expected_monoclinic_lattice_constants
        )
    end
end

@testset "convert_to_mS_case_2b(): only m_a in aP basis, |m_a| > |m_c|" begin
    println("convert_to_mS_case_2b()")
    # --- Preparations

    # Construct basis for monoclinic unit cell
    m_a = 3.0
    m_b = 2.0
    m_c = 1.0
    m_β = 3π / 5

    expected_monoclinic_lattice_constants = MonoclinicLatticeConstants(m_a, m_b, m_c, m_β)

    m_basis_a, m_basis_b, m_basis_c = basis(expected_monoclinic_lattice_constants)

    # --- Tests

    # Test all cases for base-centered basis vector
    basis_a = m_basis_a
    for signs in Iterators.product((-1, 1), (-1, 1))
        basis_b = 0.5 * (signs[1] * m_basis_a - signs[2] * m_basis_b)
        basis_c = 0.5 * (signs[1] * m_basis_a + signs[2] * m_basis_b) + signs[1] * m_basis_c

        # Check that basis vectors are linearly independent
        if !is_basis(basis_a, basis_b, basis_c)
            continue
        end

        convert_to_mS_test_all_triclinic_basis_permutations(
            basis_a, basis_b, basis_c, expected_monoclinic_lattice_constants
        )
    end
end

#=
@testset "convert_to_mS(): {m_b, m_c} in aP basis, |m_a| < |m_c|" begin
    # --- Preparations

    # Construct basis for monoclinic unit cell
    m_a = 1.0
    m_b = 2.0
    m_c = 3.0
    m_β = 3π / 5

    expected_monoclinic_lattice_constants = MonoclinicLatticeConstants(m_a, m_b, m_c, m_β)

    m_basis_a, m_basis_b, m_basis_c = basis(expected_monoclinic_lattice_constants)

    # --- Tests

    # Test all cases for base-centered basis vector
    basis_b = m_basis_b
    basis_c = m_basis_c
    for signs in Iterators.product((-1, 1), (-1, 1))
        basis_a = 0.5 * (signs[1] * m_basis_a + signs[2] * m_basis_b)

        convert_to_mS_test_all_triclinic_basis_permutations(
            basis_a, basis_b, basis_c, expected_monoclinic_lattice_constants
        )
    end
end

@testset "convert_to_mS(): {m_b, m_c} in aP basis, |m_a| > |m_c|" begin
    # --- Preparations

    # Construct basis for monoclinic unit cell
    m_a = 3.0
    m_b = 2.0
    m_c = 1.0
    m_β = 3π / 5

    expected_monoclinic_lattice_constants = MonoclinicLatticeConstants(m_a, m_b, m_c, m_β)

    m_basis_a, m_basis_b, m_basis_c = basis(expected_monoclinic_lattice_constants)

    # --- Tests

    # Test all cases for base-centered basis vector
    basis_b = m_basis_b
    basis_c = m_basis_c
    for signs in Iterators.product((-1, 1), (-1, 1))
        basis_a = 0.5 * (signs[1] * m_basis_a + signs[2] * m_basis_b)

        convert_to_mS_test_all_triclinic_basis_permutations(
            basis_a, basis_b, basis_c, expected_monoclinic_lattice_constants
        )
    end
end

@testset "convert_to_mS(): only m_c in aP basis, |m_a| < |m_c|" begin
    # --- Preparations

    # Construct basis for monoclinic unit cell
    m_a = 1.0
    m_b = 2.0
    m_c = 3.0
    m_β = 3π / 5

    expected_monoclinic_lattice_constants, _ = standardize(
        MonoclinicLatticeConstants(m_a, m_b, m_c, m_β), XtallographyUtils.BASE
    )

    m_basis_a, m_basis_b, m_basis_c = basis(expected_monoclinic_lattice_constants)

    # --- Tests

    # Test all cases for base-centered basis vector
    basis_c = m_basis_c
    for signs in Iterators.product((-1, 1), (-1, 1))
        basis_a = 0.5 * (signs[1] * m_basis_a + signs[2] * m_basis_b)
        basis_b = 0.5 * (signs[1] * m_basis_a - signs[2] * m_basis_b)

        # Check that basis vectors are linearly independent
        if !is_basis(basis_a, basis_b, basis_c)
            continue
        end

        convert_to_mS_test_all_triclinic_basis_permutations(
            basis_a, basis_b, basis_c, expected_monoclinic_lattice_constants
        )
    end
end

@testset "convert_to_mS(): only m_c in aP basis, |m_a| > |m_c|" begin
    # --- Preparations

    # Construct basis for monoclinic unit cell
    m_a = 3.0
    m_b = 2.0
    m_c = 1.0
    m_β = 3π / 5

    expected_monoclinic_lattice_constants = MonoclinicLatticeConstants(m_a, m_b, m_c, m_β)

    m_basis_a, m_basis_b, m_basis_c = basis(expected_monoclinic_lattice_constants)

    # --- Tests

    # Test all cases for base-centered basis vector
    basis_c = m_basis_c
    for signs in Iterators.product((-1, 1), (-1, 1))
        basis_a = 0.5 * (signs[1] * m_basis_a + signs[2] * m_basis_b)
        basis_b = 0.5 * (signs[1] * m_basis_a - signs[2] * m_basis_b)

        # Check that basis vectors are linearly independent
        if !is_basis(basis_a, basis_b, basis_c)
            continue
        end

        convert_to_mS_test_all_triclinic_basis_permutations(
            basis_a, basis_b, basis_c, expected_monoclinic_lattice_constants
        )
    end
end

@testset "convert_to_mS(): no mS basis vectors in aP basis, |m_a| < |m_c|" begin
    # --- Preparations

    # Construct basis for monoclinic unit cell
    m_a = 1.0
    m_b = 2.0
    m_c = 3.0
    m_β = 3π / 5

    expected_monoclinic_lattice_constants, _ = standardize(
        MonoclinicLatticeConstants(m_a, m_b, m_c, m_β), XtallographyUtils.BASE
    )

    m_basis_a, m_basis_b, m_basis_c = basis(expected_monoclinic_lattice_constants)

    # --- Tests

    # Test all cases for body-centered basis vector
    for signs in Iterators.product((-1, 1), (-1, 1), (-1, 1))
        basis_a = 0.5 * (signs[1] * m_basis_a + signs[2] * m_basis_b)
        basis_b = 0.5 * (signs[1] * m_basis_a - signs[2] * m_basis_b)
        basis_c = 0.5 * (signs[1] * m_basis_a + signs[2] * m_basis_b) + signs[3] * m_basis_c

        # Check that basis vectors are linearly independent
        if !is_basis(basis_a, basis_b, basis_c)
            continue
        end

        convert_to_mS_test_all_triclinic_basis_permutations(
            basis_a, basis_b, basis_c, expected_monoclinic_lattice_constants
        )
    end
end

@testset "convert_to_mS(): no mS basis vectors in aP basis, |m_a| > |m_c|" begin
    # --- Preparations

    # Construct basis for monoclinic unit cell
    m_a = 3.0
    m_b = 2.0
    m_c = 1.0
    m_β = 3π / 5

    # expected_monoclinic_lattice_constants = MonoclinicLatticeConstants(m_a, m_b, m_c, m_β)
    expected_monoclinic_lattice_constants, _ = standardize(
        MonoclinicLatticeConstants(m_a, m_b, m_c, m_β), XtallographyUtils.BASE
    )

    m_basis_a, m_basis_b, m_basis_c = basis(expected_monoclinic_lattice_constants)

    # --- Tests

    # Test all cases for body-centered basis vector
    for signs in Iterators.product((-1, 1), (-1, 1), (-1, 1))
        basis_a = 0.5 * (signs[1] * m_basis_a + signs[2] * m_basis_b)
        basis_b = 0.5 * (signs[1] * m_basis_a - signs[2] * m_basis_b)
        basis_c = 0.5 * (signs[1] * m_basis_a + signs[2] * m_basis_b) + signs[3] * m_basis_c

        # Check that basis vectors are linearly independent
        if !is_basis(basis_a, basis_b, basis_c)
            continue
        end

        convert_to_mS_test_all_triclinic_basis_permutations(
            basis_a, basis_b, basis_c, expected_monoclinic_lattice_constants
        )
    end
end

@testset "convert_to_mS(): aP cell not equivalent to a mS cell" begin
    # --- Tests

    # triclinic unit cell is not equivalent to a body-centered monoclinic unit cell
    a = 1.0
    b = 2.0
    c = 3.0
    α = π / 2
    β = 4π / 7
    γ = 3π / 5
    lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)

    local error = nothing
    local error_message = ""
    try
        convert_to_mS(lattice_constants)

    catch error
        bt = catch_backtrace()
        error_message = sprint(showerror, error, bt)
    end

    @test error isa ErrorException

    expected_error =
        "The triclinic unit cell defined by `lattice_constants` is not equivalent to " *
        "a base-centered monoclinic unit cell."
    @test startswith(error_message, expected_error)
end
=#

# ------ reduced_cell()

@testset "reduced_cell()" begin
    # --- Preparations

    a = 5
    b = 8
    c = 10
    α = 2π / 5
    β = 3π / 5
    γ = 4π / 5
    lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)
    basis_a, basis_b, basis_c = basis(lattice_constants)

    # --- Exercise functionality and check results

    # primitive unit cell defined by [basis_a, basis_b, basis_c]
    unit_cell = UnitCell(lattice_constants, XtallographyUtils.PRIMITIVE)

    expected_reduced_cell = reduced_cell(
        UnitCell(lattice_constants, XtallographyUtils.PRIMITIVE)
    )

    reduced_cell_ = reduced_cell(unit_cell)
    @test reduced_cell_.lattice_constants isa TriclinicLatticeConstants
    @test volume(reduced_cell_) ≈ volume(unit_cell)
    @test reduced_cell_ ≈ expected_reduced_cell

    # primitive unit cell defined by linear combination of [basis_a, basis_b, basis_c]
    unit_cell = UnitCell(
        LatticeConstants(basis_a + basis_b + basis_c, basis_b, basis_c),
        XtallographyUtils.PRIMITIVE,
    )

    expected_reduced_cell = reduced_cell(
        UnitCell(lattice_constants, XtallographyUtils.PRIMITIVE)
    )

    reduced_cell_ = reduced_cell(unit_cell)
    @test reduced_cell_.lattice_constants isa TriclinicLatticeConstants
    @test volume(reduced_cell_) ≈ volume(unit_cell)
    @test reduced_cell_ ≈ expected_reduced_cell

    #=
    # TODO: Delaunay reduction does not terminate
    a = 5
    b = 8
    c = 10
    α = π / 5
    β = 2π / 5
    γ = 3π / 5
    lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)
    basis_a, basis_b, basis_c = basis(lattice_constants)
    unit_cell = UnitCell(
        LatticeConstants(basis_a + basis_b + basis_c, basis_b, basis_c), XtallographyUtils.PRIMITIVE
    )

    expected_reduced_cell = reduced_cell(UnitCell(lattice_constants, XtallographyUtils.PRIMITIVE))

    reduced_cell_ = reduced_cell(unit_cell)
    =#
end

# ------ is_triclinic_type_I_cell()

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
    @test is_triclinic_type_I_cell(TriclinicLatticeConstants(a, b, c, α, β, γ))

    # α equal to π/2
    a = 1
    b = 1
    c = 1
    α = π / 2
    β = 3π / 10
    γ = 4π / 10
    @test is_triclinic_type_I_cell(TriclinicLatticeConstants(a, b, c, α, β, γ))

    # α equal to π/2 + ϵ (so that cos(α) < 0)
    a = 1
    b = 1
    c = 1
    α = π / 2 + 1e-9
    β = 3π / 10
    γ = 4π / 10
    @test is_triclinic_type_I_cell(TriclinicLatticeConstants(a, b, c, α, β, γ))

    # ------ Type II

    # Type II: no angles equal to π/2
    a = 1
    b = 1
    c = 1
    α = 6π / 10
    β = 7π / 10
    γ = 8π / 10
    @test !is_triclinic_type_I_cell(TriclinicLatticeConstants(a, b, c, α, β, γ))

    # γ equal to π/2
    a = 1
    b = 1
    c = 1
    α = 8π / 10
    β = 7π / 10
    γ = π / 2
    @test is_triclinic_type_I_cell(TriclinicLatticeConstants(a, b, c, α, β, γ))

    # β equal to π/2 + ϵ (so that cos(β) < 0)
    a = 1
    b = 1
    c = 1
    α = 7π / 10
    β = π / 2 + 1e-9
    γ = 8π / 10
    @test is_triclinic_type_I_cell(TriclinicLatticeConstants(a, b, c, α, β, γ))
end
