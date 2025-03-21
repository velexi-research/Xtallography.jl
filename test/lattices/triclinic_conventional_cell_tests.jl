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
using LinearAlgebra: norm, dot
using Test

# External packages
using Combinatorics: combinations, permutations
using AngleBetweenVectors: angle

# Xtallography package
using Xtallography

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

@testset "conventional_cell():triclinic: limiting cases" begin
    # Note: this test set only tests basic functionality. Comprehensive coverage of
    #       different cases for triclinic bases is covered in the convert_to_mP()
    #       convert_to_mI(), and convert_to_mC() test sets.

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
    iucr_unit_cell = conventional_cell(UnitCell(lattice_constants, primitive))

    # Check results
    @test iucr_unit_cell.lattice_constants ≈ standardize(lattice_constants)
    @test iucr_unit_cell.centering === primitive

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
    iucr_unit_cell = conventional_cell(UnitCell(lattice_constants, primitive))

    # Check results
    @test iucr_unit_cell.lattice_constants ≈
        standardize(MonoclinicLatticeConstants(b, a, c, π - α))
    @test iucr_unit_cell.centering === primitive

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
    iucr_unit_cell = conventional_cell(UnitCell(lattice_constants, primitive))

    # Check results
    expected_lattice_constants, _ = standardize(
        MonoclinicLatticeConstants(m_a, m_b, m_c, m_β), body_centered
    )
    @test iucr_unit_cell.lattice_constants ≈ expected_lattice_constants
    @test iucr_unit_cell.centering === body_centered
end

@testset "conventional_cell():triclinic: limiting cases - base-centered" begin
    # Note: this test checks that all base-centered limiting cases are covered by
    #       convert_to_mI() and convert_to_mC().

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
    iucr_unit_cell = conventional_cell(UnitCell(lattice_constants, primitive))

    # Check results
    expected_lattice_constants, _ = standardize(
        MonoclinicLatticeConstants(m_a, m_b, m_c, m_β), base_centered
    )
    @test iucr_unit_cell.lattice_constants ≈ expected_lattice_constants
    @test iucr_unit_cell.centering === body_centered

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
    iucr_unit_cell = conventional_cell(UnitCell(lattice_constants, primitive))

    # Check results
    expected_lattice_constants, _ = standardize(
        MonoclinicLatticeConstants(m_a, m_b, m_c, m_β), base_centered
    )
    @test iucr_unit_cell.lattice_constants ≈ expected_lattice_constants
    @test iucr_unit_cell.centering === body_centered

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
    iucr_unit_cell = conventional_cell(UnitCell(lattice_constants, primitive))

    # Check results
    expected_lattice_constants, _ = standardize(
        MonoclinicLatticeConstants(m_a, m_b, m_c, m_β), base_centered
    )
    @test iucr_unit_cell.lattice_constants ≈ expected_lattice_constants
    @test iucr_unit_cell.centering === body_centered

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
    iucr_unit_cell = conventional_cell(UnitCell(lattice_constants, primitive))

    # Check results
    expected_lattice_constants, _ = standardize(
        MonoclinicLatticeConstants(m_a, m_b, m_c, m_β), base_centered
    )
    @test iucr_unit_cell.lattice_constants ≈ expected_lattice_constants
    @test iucr_unit_cell.centering === body_centered

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
    iucr_unit_cell = conventional_cell(UnitCell(lattice_constants, primitive))

    # Check results
    expected_lattice_constants, _ = standardize(
        MonoclinicLatticeConstants(m_a, m_b, m_c, m_β), base_centered
    )
    @test iucr_unit_cell.lattice_constants ≈ expected_lattice_constants
    @test iucr_unit_cell.centering === body_centered

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
    iucr_unit_cell = conventional_cell(UnitCell(lattice_constants, primitive))

    # Check results
    expected_lattice_constants, _ = standardize(
        MonoclinicLatticeConstants(m_a, m_b, m_c, m_β), base_centered
    )
    @test iucr_unit_cell.lattice_constants ≈ expected_lattice_constants
    @test iucr_unit_cell.centering === body_centered

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
    iucr_unit_cell = conventional_cell(UnitCell(lattice_constants, primitive))

    # Check results
    expected_lattice_constants, _ = standardize(
        MonoclinicLatticeConstants(m_a, m_b, m_c, m_β), base_centered
    )
    @test iucr_unit_cell.lattice_constants ≈ expected_lattice_constants
    @test iucr_unit_cell.centering === body_centered

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
    iucr_unit_cell = conventional_cell(UnitCell(lattice_constants, primitive))

    # Check results
    expected_lattice_constants, _ = standardize(
        MonoclinicLatticeConstants(m_a, m_b, m_c, m_β), base_centered
    )
    @test iucr_unit_cell.lattice_constants ≈ expected_lattice_constants
    @test iucr_unit_cell.centering === body_centered
end

@testset "conventional_cell():triclinic: invalid arguments" begin
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

    # ------ Invalid centering

    for centering in (body_centered, face_centered, base_centered)
        expected_message =
            "Invalid Bravais lattice: " *
            "(lattice_system=Triclinic, centering=$(nameof(typeof(centering))))"

        @test_throws ArgumentError(expected_message) conventional_cell(
            UnitCell(lattice_constants, centering)
        )
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
        expected_monoclinic_lattice_constants, body_centered
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
        monoclinic_unit_cell = UnitCell(monoclinic_lattice_constants, body_centered)

        # Check results
        @test is_equivalent_unit_cell(monoclinic_unit_cell, expected_monoclinic_unit_cell)
    end
end

@testset "convert_to_mI_case_1(): {m_a, m_b} in aP basis, (m_c/m_a) |cos m_β| < 1" begin
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

@testset "convert_to_mI_case_1(): {m_a, m_b} in aP basis, (m_c/m_a) |cos m_β| > 1" begin
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

@testset "convert_to_mI_case_1(): {m_b, m_c} in aP basis, (m_c/m_a) |cos m_β| < 1" begin
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

@testset "convert_to_mI_case_1(): {m_b, m_c} in aP basis, (m_c/m_a) |cos m_β| > 1" begin
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

@testset "convert_to_mI_case_2(): {m_a, m_c} in aP basis, (m_c/m_a) |cos m_β| < 1" begin
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
    for signs in Iterators.product((-1, 1), (-1, 1), (-1, 1), (-1, 1))
        basis_b = 0.5 * (signs[1] * basis_a + signs[2] * m_basis_b + signs[3] * m_basis_c)
        basis_c = signs[4] * m_basis_c

        convert_to_mI_test_all_triclinic_basis_permutations(
            basis_a, basis_b, basis_c, expected_monoclinic_lattice_constants
        )
    end
end

@testset "convert_to_mI_case_2(): {m_a, m_c} in aP basis, (m_c/m_a) |cos m_β| > 1" begin
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
    for signs in Iterators.product((-1, 1), (-1, 1), (-1, 1), (-1, 1))
        basis_b = 0.5 * (signs[1] * basis_a + signs[2] * m_basis_b + signs[3] * m_basis_c)
        basis_c = signs[4] * m_basis_c

        convert_to_mI_test_all_triclinic_basis_permutations(
            basis_a, basis_b, basis_c, expected_monoclinic_lattice_constants
        )
    end
end

@testset "convert_to_mI_case_4(): only m_a in aP basis, (m_c/m_a) |cos m_β| < 1" begin
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

@testset "convert_to_mI_case_4(): only m_a in aP basis, (m_c/m_a) |cos m_β| > 1" begin
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

@testset "convert_to_mI_case_3(): only m_b in aP basis, (m_c/m_a) |cos m_β| < 1" begin
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

@testset "convert_to_mI_case_3(): only m_b in aP basis, (m_c/m_a) |cos m_β| > 1" begin
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

@testset "convert_to_mI_case_4(): only m_c in aP basis, (m_c/m_a) |cos m_β| < 1" begin
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

@testset "convert_to_mI_case_4(): only m_c in aP basis, (m_c/m_a) |cos m_β| > 1" begin
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

@testset "convert_to_mI_case_5(): no mI basis vectors in aP basis, (m_c/m_a) |cos m_β| < 1" begin
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

@testset "convert_to_mI_case_5(): no mI basis vectors in aP basis, (m_c/m_a) |cos m_β| > 1" begin
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

# ------ convert_to_mC()

# Helper functions
function convert_to_mC_test_all_triclinic_basis_permutations(
    basis_a::Vector{<:Real},
    basis_b::Vector{<:Real},
    basis_c::Vector{<:Real},
    expected_monoclinic_lattice_constants::MonoclinicLatticeConstants,
)
    # --- Preparations

    # Construct expected unit cell
    expected_monoclinic_unit_cell = UnitCell(
        expected_monoclinic_lattice_constants, base_centered
    )

    # --- Test convert_to_mC() for each permutation of the basis vectors

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
        monoclinic_lattice_constants = convert_to_mC(lattice_constants)
        monoclinic_unit_cell = UnitCell(monoclinic_lattice_constants, base_centered)

        # Check results
        @test is_equivalent_unit_cell(monoclinic_unit_cell, expected_monoclinic_unit_cell)
    end
end

@testset "convert_to_mC_case_1a(): {m_a, m_c} in aP basis, |m_a| < |m_c|" begin
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

        convert_to_mC_test_all_triclinic_basis_permutations(
            basis_a, basis_b, basis_c, expected_monoclinic_lattice_constants
        )
    end
end

@testset "convert_to_mC_case_1a(): {m_a, m_c} in aP basis, |m_a| > |m_c|" begin
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

        convert_to_mC_test_all_triclinic_basis_permutations(
            basis_a, basis_b, basis_c, expected_monoclinic_lattice_constants
        )
    end
end

@testset "convert_to_mC_case_1b(): {m_a, m_c} in aP basis, |m_a| < |m_c|" begin
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

        convert_to_mC_test_all_triclinic_basis_permutations(
            basis_a, basis_b, basis_c, expected_monoclinic_lattice_constants
        )
    end
end

@testset "convert_to_mC_case_1b(): {m_a, m_c} in aP basis, |m_a| > |m_c|" begin
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

        convert_to_mC_test_all_triclinic_basis_permutations(
            basis_a, basis_b, basis_c, expected_monoclinic_lattice_constants
        )
    end
end

@testset "convert_to_mC_case_2a(): only m_a in aP basis" begin
    # Case 2a:
    # * the sign of the coefficient of m_basis_a in basis_b is positive
    # * the signs of the coefficients of m_basis_b in basis_b and basis_c are the same

    # --- Tests: |m_a| < |m_c|

    # Construct basis for monoclinic unit cell
    m_a = 1.0
    m_b = 2.0
    m_c = 3.0
    m_β = 3π / 5

    expected_monoclinic_lattice_constants, _ = standardize(
        MonoclinicLatticeConstants(m_a, m_b, m_c, m_β), base_centered
    )

    m_basis_a, m_basis_b, m_basis_c = basis(expected_monoclinic_lattice_constants)

    basis_a = m_basis_a
    for signs in Iterators.product((-1, 1), (-1, 1), (-1, 1))
        basis_b = 0.5 * (m_basis_a + signs[2] * m_basis_b)
        basis_c = 0.5 * (signs[1] * m_basis_a + signs[2] * m_basis_b) + signs[3] * m_basis_c

        # Check that basis vectors are linearly independent
        if !is_basis(basis_a, basis_b, basis_c)
            continue
        end

        convert_to_mC_test_all_triclinic_basis_permutations(
            basis_a, basis_b, basis_c, expected_monoclinic_lattice_constants
        )
    end

    # --- Tests: |m_a| > |m_c|

    # Construct basis for monoclinic unit cell
    m_a = 3.0
    m_b = 2.0
    m_c = 1.0
    m_β = 3π / 5

    expected_monoclinic_lattice_constants, _ = standardize(
        MonoclinicLatticeConstants(m_a, m_b, m_c, m_β), base_centered
    )

    m_basis_a, m_basis_b, m_basis_c = basis(expected_monoclinic_lattice_constants)

    basis_a = m_basis_a
    for signs in Iterators.product((-1, 1), (-1, 1), (-1, 1))
        basis_b = 0.5 * (m_basis_a + signs[2] * m_basis_b)
        basis_c = 0.5 * (signs[1] * m_basis_a + signs[2] * m_basis_b) + signs[3] * m_basis_c

        # Check that basis vectors are linearly independent
        if !is_basis(basis_a, basis_b, basis_c)
            continue
        end

        convert_to_mC_test_all_triclinic_basis_permutations(
            basis_a, basis_b, basis_c, expected_monoclinic_lattice_constants
        )
    end
end

@testset "convert_to_mC_case_2b(): only m_a in aP basis" begin
    # Case 2b:
    # * the sign of the coefficient of m_basis_a in basis_b is positive
    # * the signs of the coefficients of m_basis_b in basis_b and basis_c are opposite

    # --- Tests: |m_a| < |m_c|

    # Construct basis for monoclinic unit cell
    m_a = 1.0
    m_b = 2.0
    m_c = 3.0
    m_β = 3π / 5

    expected_monoclinic_lattice_constants, _ = standardize(
        MonoclinicLatticeConstants(m_a, m_b, m_c, m_β), base_centered
    )

    m_basis_a, m_basis_b, m_basis_c = basis(expected_monoclinic_lattice_constants)

    basis_a = m_basis_a
    for signs in Iterators.product((-1, 1), (-1, 1), (-1, 1))
        basis_b = 0.5 * (m_basis_a + signs[2] * m_basis_b)
        basis_c = 0.5 * (signs[1] * m_basis_a - signs[2] * m_basis_b) + signs[3] * m_basis_c

        # Check that basis vectors are linearly independent
        if !is_basis(basis_a, basis_b, basis_c)
            continue
        end

        convert_to_mC_test_all_triclinic_basis_permutations(
            basis_a, basis_b, basis_c, expected_monoclinic_lattice_constants
        )
    end

    # --- Tests: |m_a| > |m_c|

    # Construct basis for monoclinic unit cell
    m_a = 3.0
    m_b = 2.0
    m_c = 1.0
    m_β = 3π / 5

    expected_monoclinic_lattice_constants, _ = standardize(
        MonoclinicLatticeConstants(m_a, m_b, m_c, m_β), base_centered
    )

    m_basis_a, m_basis_b, m_basis_c = basis(expected_monoclinic_lattice_constants)

    basis_a = m_basis_a
    for signs in Iterators.product((-1, 1), (-1, 1), (-1, 1))
        basis_b = 0.5 * (m_basis_a + signs[2] * m_basis_b)
        basis_c = 0.5 * (signs[1] * m_basis_a - signs[2] * m_basis_b) + signs[3] * m_basis_c

        # Check that basis vectors are linearly independent
        if !is_basis(basis_a, basis_b, basis_c)
            continue
        end

        convert_to_mC_test_all_triclinic_basis_permutations(
            basis_a, basis_b, basis_c, expected_monoclinic_lattice_constants
        )
    end
end

@testset "convert_to_mC_case_2c(): only m_a in aP basis" begin
    # Case 2c:
    # * the sign of the coefficient of m_basis_a in basis_b is negative
    # * the signs of the coefficients of m_basis_b in basis_b and basis_c are the same

    # --- Tests: |m_a| < |m_c|

    # Construct basis for monoclinic unit cell
    m_a = 1.0
    m_b = 2.0
    m_c = 3.0
    m_β = 3π / 5

    expected_monoclinic_lattice_constants, _ = standardize(
        MonoclinicLatticeConstants(m_a, m_b, m_c, m_β), base_centered
    )

    m_basis_a, m_basis_b, m_basis_c = basis(expected_monoclinic_lattice_constants)

    basis_a = m_basis_a
    for signs in Iterators.product((-1, 1), (-1, 1), (-1, 1))
        basis_b = 0.5 * (-m_basis_a + signs[2] * m_basis_b)
        basis_c = 0.5 * (signs[1] * m_basis_a + signs[2] * m_basis_b) + signs[3] * m_basis_c

        # Check that basis vectors are linearly independent
        if !is_basis(basis_a, basis_b, basis_c)
            continue
        end

        convert_to_mC_test_all_triclinic_basis_permutations(
            basis_a, basis_b, basis_c, expected_monoclinic_lattice_constants
        )
    end

    # --- Tests: |m_a| > |m_c|

    # Construct basis for monoclinic unit cell
    m_a = 3.0
    m_b = 2.0
    m_c = 1.0
    m_β = 3π / 5

    expected_monoclinic_lattice_constants, _ = standardize(
        MonoclinicLatticeConstants(m_a, m_b, m_c, m_β), base_centered
    )

    m_basis_a, m_basis_b, m_basis_c = basis(expected_monoclinic_lattice_constants)

    basis_a = m_basis_a
    for signs in Iterators.product((-1, 1), (-1, 1), (-1, 1))
        basis_b = 0.5 * (-m_basis_a + signs[2] * m_basis_b)
        basis_c = 0.5 * (signs[1] * m_basis_a + signs[2] * m_basis_b) + signs[3] * m_basis_c

        # Check that basis vectors are linearly independent
        if !is_basis(basis_a, basis_b, basis_c)
            continue
        end

        convert_to_mC_test_all_triclinic_basis_permutations(
            basis_a, basis_b, basis_c, expected_monoclinic_lattice_constants
        )
    end
end

@testset "convert_to_mC_case_2d(): only m_a in aP basis" begin
    # Case 2d:
    # * the sign of the coefficient of m_basis_a in basis_b is negative
    # * the signs of the coefficients of m_basis_b in basis_b and basis_c are opposite

    # --- Tests: |m_a| < |m_c|

    # Construct basis for monoclinic unit cell
    m_a = 1.0
    m_b = 2.0
    m_c = 3.0
    m_β = 3π / 5

    expected_monoclinic_lattice_constants, _ = standardize(
        MonoclinicLatticeConstants(m_a, m_b, m_c, m_β), base_centered
    )

    m_basis_a, m_basis_b, m_basis_c = basis(expected_monoclinic_lattice_constants)

    basis_a = m_basis_a
    for signs in Iterators.product((-1, 1), (-1, 1), (-1, 1))
        basis_b = 0.5 * (-m_basis_a + signs[2] * m_basis_b)
        basis_c = 0.5 * (signs[1] * m_basis_a - signs[2] * m_basis_b) + signs[3] * m_basis_c

        # Check that basis vectors are linearly independent
        if !is_basis(basis_a, basis_b, basis_c)
            continue
        end

        convert_to_mC_test_all_triclinic_basis_permutations(
            basis_a, basis_b, basis_c, expected_monoclinic_lattice_constants
        )
    end

    # --- Tests: |m_a| > |m_c|

    # Construct basis for monoclinic unit cell
    m_a = 3.0
    m_b = 2.0
    m_c = 1.0
    m_β = 3π / 5

    expected_monoclinic_lattice_constants, _ = standardize(
        MonoclinicLatticeConstants(m_a, m_b, m_c, m_β), base_centered
    )

    m_basis_a, m_basis_b, m_basis_c = basis(expected_monoclinic_lattice_constants)

    basis_a = m_basis_a
    for signs in Iterators.product((-1, 1), (-1, 1), (-1, 1))
        basis_b = 0.5 * (-m_basis_a + signs[2] * m_basis_b)
        basis_c = 0.5 * (signs[1] * m_basis_a - signs[2] * m_basis_b) + signs[3] * m_basis_c

        # Check that basis vectors are linearly independent
        if !is_basis(basis_a, basis_b, basis_c)
            continue
        end

        convert_to_mC_test_all_triclinic_basis_permutations(
            basis_a, basis_b, basis_c, expected_monoclinic_lattice_constants
        )
    end
end

@testset "convert_to_mC(): aP cell not equivalent to a mC cell" begin
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
        convert_to_mC(lattice_constants)

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
