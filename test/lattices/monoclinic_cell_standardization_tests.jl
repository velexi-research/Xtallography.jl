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
Tests for the unit cell standardization methods for monoclinic lattices
"""
# --- Imports

# Standard library
using Test
using LinearAlgebra: dot, norm

# XtallographyUtils package
using XtallographyUtils

# Notes
# =====
# These tests adopt the following variable conventions.
#
# - Unless otherwise noted, lattice constants and basis vectors refer to the monoclinic
#   (not orthorhombic or rhombohedral) unit cell.
#
# - Lattice constants and basis vectors for orthorhombic unit cells are indicated by
#   the "o_" prefix.
#
# - Lattice constants and basis vectors for rhombohedral unit cells are indicated by
#   the "r_" prefix.

# --- Tests

# ------ standardize()

@testset "standardize()" begin
    # --- Tests

    # ------ lattice constants already in standard form

    a = 6
    b = 10
    c = 8
    β = 1.1 * π / 2
    lattice_constants = MonoclinicLatticeConstants(a, b, c, β)

    # centering = PRIMITIVE
    standardized_lattice_constants, standardized_centering = standardize(
        lattice_constants, XtallographyUtils.PRIMITIVE
    )

    expected_lattice_constants = MonoclinicLatticeConstants(a, b, c, β)
    @test standardized_lattice_constants ≈ expected_lattice_constants

    @test standardized_centering == XtallographyUtils.PRIMITIVE

    # centering = BODY
    standardized_lattice_constants, standardized_centering = standardize(
        lattice_constants, XtallographyUtils.BODY
    )

    expected_lattice_constants = MonoclinicLatticeConstants(a, b, c, β)
    @test standardized_lattice_constants ≈ expected_lattice_constants

    @test standardized_centering == XtallographyUtils.BODY

    # ------ β ∉ [π/2, π]

    a = 6
    b = 10
    c = 8
    β = π - 1.1 * π / 2
    lattice_constants = MonoclinicLatticeConstants(a, b, c, β)

    # centering = PRIMITIVE
    standardized_lattice_constants, standardized_centering = standardize(
        lattice_constants, XtallographyUtils.PRIMITIVE
    )

    expected_lattice_constants = MonoclinicLatticeConstants(a, b, c, π - β)
    @test standardized_lattice_constants ≈ expected_lattice_constants

    @test standardized_centering == XtallographyUtils.PRIMITIVE

    # centering = BODY
    standardized_lattice_constants, standardized_centering = standardize(
        lattice_constants, XtallographyUtils.BODY
    )

    expected_lattice_constants = MonoclinicLatticeConstants(a, b, c, π - β)
    @test standardized_lattice_constants ≈ expected_lattice_constants

    @test standardized_centering == XtallographyUtils.BODY

    # ------ a > c

    a = 8
    b = 10
    c = 6
    β = 1.1 * π / 2
    lattice_constants = MonoclinicLatticeConstants(a, b, c, β)

    # centering = PRIMITIVE
    standardized_lattice_constants, standardized_centering = standardize(
        lattice_constants, XtallographyUtils.PRIMITIVE
    )

    expected_lattice_constants = MonoclinicLatticeConstants(c, b, a, β)
    @test standardized_lattice_constants ≈ expected_lattice_constants

    @test standardized_centering == XtallographyUtils.PRIMITIVE

    # centering = BODY
    standardized_lattice_constants, standardized_centering = standardize(
        lattice_constants, XtallographyUtils.BODY
    )

    expected_lattice_constants = MonoclinicLatticeConstants(c, b, a, β)
    @test standardized_lattice_constants ≈ expected_lattice_constants

    @test standardized_centering == XtallographyUtils.BODY

    # ------ base-centered unit cell, a < c

    a = 6
    b = 10
    c = 8
    β = 1.1 * π / 2
    lattice_constants = MonoclinicLatticeConstants(a, b, c, β)

    standardized_lattice_constants, standardized_centering = standardize(
        lattice_constants, XtallographyUtils.BASE
    )

    a_body = sqrt(a^2 + c^2 + 2 * a * c * cos(β))
    β_body = π - asin(sin(β) / a_body * a)
    expected_lattice_constants = MonoclinicLatticeConstants(c, b, a_body, β_body)
    @test standardized_lattice_constants ≈ expected_lattice_constants

    @test standardized_centering == XtallographyUtils.BODY

    # ------ base-centered unit cell, a > c

    a = 8
    b = 10
    c = 6
    β = 1.1 * π / 2
    lattice_constants = MonoclinicLatticeConstants(a, b, c, β)

    standardized_lattice_constants, standardized_centering = standardize(
        lattice_constants, XtallographyUtils.BASE
    )

    a_body = sqrt(a^2 + c^2 + 2 * a * c * cos(β))
    β_body = π - asin(sin(β) / a_body * a)
    expected_lattice_constants = MonoclinicLatticeConstants(c, b, a_body, β_body)
    @test standardized_lattice_constants ≈ expected_lattice_constants

    @test standardized_centering == XtallographyUtils.BODY

    # ------ unit cell requires single reduction in plane normal to b-axis

    a_ref = 6
    b_ref = 10
    c_ref = 8
    β_ref = 2π / 3

    # centering = PRIMITIVE
    a = a_ref
    b = b_ref
    c = c_ref
    β = β_ref
    lattice_constants = MonoclinicLatticeConstants(a, b, c, β)

    standardized_lattice_constants, standardized_centering = standardize(
        lattice_constants, XtallographyUtils.PRIMITIVE
    )

    expected_c = sqrt(a^2 + c^2 + 2 * a * c * cos(β))
    expected_β = π - asin(sin(β) / expected_c * c)
    expected_lattice_constants = MonoclinicLatticeConstants(a, b, expected_c, expected_β)
    @test standardized_lattice_constants ≈ expected_lattice_constants

    @test standardized_centering == XtallographyUtils.PRIMITIVE

    # centering = BODY
    a = a_ref
    b = b_ref
    c = sqrt((2 * a_ref)^2 + c_ref^2 + 2 * (2 * a_ref) * c_ref * cos(β_ref))
    β = π - asin(sin(β_ref) / c * c_ref)
    lattice_constants = MonoclinicLatticeConstants(a, b, c, β)

    standardized_lattice_constants, standardized_centering = standardize(
        lattice_constants, XtallographyUtils.BODY
    )

    expected_lattice_constants = MonoclinicLatticeConstants(a_ref, b_ref, c_ref, β_ref)
    @test standardized_lattice_constants ≈ expected_lattice_constants

    @test standardized_centering == XtallographyUtils.BODY

    # ------ unit cell requires multiple reductions in plane normal to b-axis

    a_ref = 6
    b_ref = 10
    c_ref = 8
    β_ref = 2π / 3

    # centering = PRIMITIVE
    a = a_ref
    b = b_ref
    c = sqrt((5 * a_ref)^2 + c_ref^2 + 2 * (5 * a_ref) * c_ref * cos(β_ref))
    β = asin(sin(β_ref) / c * c_ref)
    lattice_constants = MonoclinicLatticeConstants(a, b, c, β)

    standardized_lattice_constants, standardized_centering = standardize(
        lattice_constants, XtallographyUtils.PRIMITIVE
    )

    expected_c = sqrt(a_ref^2 + c_ref^2 + 2 * a_ref * c_ref * cos(β_ref))
    expected_β = π - asin(sin(β_ref) / expected_c * c_ref)
    expected_lattice_constants = MonoclinicLatticeConstants(a, b, expected_c, expected_β)
    @test standardized_lattice_constants ≈ expected_lattice_constants

    @test standardized_centering == XtallographyUtils.PRIMITIVE

    # centering = BODY
    a = a_ref
    b = b_ref
    c = sqrt((5 * a_ref)^2 + c_ref^2 + 2 * (5 * a_ref) * c_ref * cos(β_ref))
    β = asin(sin(β_ref) / c * c_ref)
    lattice_constants = MonoclinicLatticeConstants(a, b, c, β)

    standardized_lattice_constants, standardized_centering = standardize(
        lattice_constants, XtallographyUtils.BODY
    )

    expected_lattice_constants = MonoclinicLatticeConstants(a, b, expected_c, expected_β)
    @test standardized_lattice_constants ≈ expected_lattice_constants

    @test standardized_centering == XtallographyUtils.BODY

    # ------ Invalid centering

    # centering = FACE
    centering = XtallographyUtils.FACE
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
        "Invalid Bravais lattice: (lattice_system=Monoclinic, centering=$centering)"

    @test startswith(error_message, expected_error)
end

# ------ iucr_conventional_cell()

@testset "iucr_conventional_cell(): limiting cases, centering = PRIMITIVE" begin
    # --- Preparations

    # Construct basis for orthorhombic unit cell
    o_a = 1.0
    o_b = 2.0
    o_c = 3.0

    expected_orthorhombic_lattice_constants = OrthorhombicLatticeConstants(o_a, o_b, o_c)

    o_basis_a, o_basis_b, o_basis_c = basis(expected_orthorhombic_lattice_constants)

    # --- Tests

    # ------ β = π / 2

    # Construct lattice constants for monoclinic unit cell
    a = o_a
    b = o_b
    c = o_c
    β = π / 2
    lattice_constants = MonoclinicLatticeConstants(a, b, c, β)

    # Exercise functionality
    iucr_unit_cell = iucr_conventional_cell(
        UnitCell(lattice_constants, XtallographyUtils.PRIMITIVE)
    )

    # Check results
    @test iucr_unit_cell.lattice_constants ≈ OrthorhombicLatticeConstants(o_a, o_b, o_c)
    @test iucr_unit_cell.centering == XtallographyUtils.PRIMITIVE

    # ------ a = -2 * c * cos(β)

    # Construct lattice constants for monoclinic unit cell
    basis_a = o_basis_a
    basis_b = 0.5 * (o_basis_a + o_basis_b)
    basis_c = o_basis_c
    lattice_constants, _ = standardize(
        LatticeConstants(basis_a, basis_b, basis_c), XtallographyUtils.PRIMITIVE
    )

    # Check test conditions
    @test lattice_constants isa MonoclinicLatticeConstants
    @test lattice_constants.a ≈ -2 * lattice_constants.c * cos(lattice_constants.β)

    # Exercise functionality
    iucr_unit_cell = iucr_conventional_cell(
        UnitCell(lattice_constants, XtallographyUtils.PRIMITIVE)
    )

    # Check results
    standardized_lattice_constants, _ = standardize(
        expected_orthorhombic_lattice_constants, XtallographyUtils.BASE
    )
    @test iucr_unit_cell.lattice_constants ≈ standardized_lattice_constants
    @test iucr_unit_cell.centering == XtallographyUtils.BASE

    # ------ a = c, basis_b = o_basis_c

    # Construct lattice constants for monoclinic unit cell
    basis_a = 0.5 * (o_basis_a + o_basis_b)
    basis_b = o_basis_c
    basis_c = 0.5 * (o_basis_a - o_basis_b)
    lattice_constants = MonoclinicLatticeConstants(
        norm(basis_a),
        norm(basis_b),
        norm(basis_c),
        π - acos(dot(basis_a, basis_c) / norm(basis_a) / norm(basis_c)),
    )

    # Check test conditions
    @test lattice_constants isa MonoclinicLatticeConstants
    @test lattice_constants.a ≈ lattice_constants.c

    # Exercise functionality
    iucr_unit_cell = iucr_conventional_cell(
        UnitCell(lattice_constants, XtallographyUtils.PRIMITIVE)
    )

    # Check results
    standardized_lattice_constants, _ = standardize(
        expected_orthorhombic_lattice_constants, XtallographyUtils.BASE
    )
    @test iucr_unit_cell.lattice_constants ≈ standardized_lattice_constants
    @test iucr_unit_cell.centering == XtallographyUtils.BASE

    # ------ monoclinic unit cell is not equivalent to an orthorhombic unit cell

    # Construct lattice constants for monoclinic unit cell
    a = 1.0
    b = 2.0
    c = 3.0
    β = 4π / 7
    lattice_constants = MonoclinicLatticeConstants(a, b, c, β)

    # Exercise functionality
    iucr_unit_cell = iucr_conventional_cell(
        UnitCell(lattice_constants, XtallographyUtils.PRIMITIVE)
    )

    # Check results
    standardized_lattice_constants, _ = standardize(
        lattice_constants, XtallographyUtils.PRIMITIVE
    )
    @test iucr_unit_cell.lattice_constants ≈ standardized_lattice_constants
    @test iucr_unit_cell.centering == XtallographyUtils.PRIMITIVE
end

@testset "iucr_conventional_cell(): limiting cases, centering = BODY, orthorhombic limits" begin
    # --- Preparations

    # Construct basis for orthorhombic unit cell
    o_a = 1.0
    o_b = 2.0
    o_c = 3.0

    expected_orthorhombic_lattice_constants = OrthorhombicLatticeConstants(o_a, o_b, o_c)

    o_basis_a, o_basis_b, o_basis_c = basis(expected_orthorhombic_lattice_constants)

    # --- Tests

    # ------ β = π / 2

    # Construct lattice constants for monoclinic unit cell
    a = o_a
    b = o_b
    c = o_c
    β = π / 2
    lattice_constants = MonoclinicLatticeConstants(a, b, c, β)

    # Exercise functionality
    iucr_unit_cell = iucr_conventional_cell(
        UnitCell(lattice_constants, XtallographyUtils.BODY)
    )

    # Check results
    @test iucr_unit_cell.lattice_constants ≈ OrthorhombicLatticeConstants(a, b, c)
    @test iucr_unit_cell.centering == XtallographyUtils.BODY

    # ------ a = -c cos(β), basis_b = o_basis_a

    # Construct lattice constants for monoclinic unit cell
    basis_a = o_basis_c
    basis_b = o_basis_a
    basis_c = o_basis_c - o_basis_b
    lattice_constants = MonoclinicLatticeConstants(
        norm(basis_a),
        norm(basis_b),
        norm(basis_c),
        π - acos(dot(basis_a, basis_c) / norm(basis_a) / norm(basis_c)),
    )

    # Check test conditions
    @test lattice_constants isa MonoclinicLatticeConstants
    @test lattice_constants.a ≈ -lattice_constants.c * cos(lattice_constants.β)

    # Exercise functionality
    iucr_unit_cell = iucr_conventional_cell(
        UnitCell(lattice_constants, XtallographyUtils.BODY)
    )

    # Check results
    standardized_lattice_constants, _ = standardize(
        expected_orthorhombic_lattice_constants, XtallographyUtils.BASE
    )
    @test iucr_unit_cell.lattice_constants ≈ standardized_lattice_constants
    @test iucr_unit_cell.centering == XtallographyUtils.BASE

    # ------ a = -c cos(β), basis_b = o_basis_b

    # Construct lattice constants for monoclinic unit cell
    basis_a = o_basis_c
    basis_b = o_basis_b
    basis_c = o_basis_c - o_basis_a
    lattice_constants = MonoclinicLatticeConstants(
        norm(basis_a),
        norm(basis_b),
        norm(basis_c),
        π - acos(dot(basis_a, basis_c) / norm(basis_a) / norm(basis_c)),
    )

    # Check test conditions
    @test lattice_constants isa MonoclinicLatticeConstants
    @test lattice_constants.a ≈ -lattice_constants.c * cos(lattice_constants.β)

    # Exercise functionality
    iucr_unit_cell = iucr_conventional_cell(
        UnitCell(lattice_constants, XtallographyUtils.BODY)
    )

    # Check results
    standardized_lattice_constants, _ = standardize(
        expected_orthorhombic_lattice_constants, XtallographyUtils.BASE
    )
    @test iucr_unit_cell.lattice_constants ≈ standardized_lattice_constants
    @test iucr_unit_cell.centering == XtallographyUtils.BASE

    # ------ a = c

    # Construct lattice constants for monoclinic unit cell
    basis_a = 0.5 * (o_basis_a + o_basis_c)
    basis_b = o_basis_b
    basis_c = 0.5 * (o_basis_a - o_basis_c)
    lattice_constants = MonoclinicLatticeConstants(
        norm(basis_a),
        norm(basis_b),
        norm(basis_c),
        π - acos(dot(basis_a, basis_c) / norm(basis_a) / norm(basis_c)),
    )

    # Check test conditions
    @test lattice_constants isa MonoclinicLatticeConstants
    @test lattice_constants.a ≈ lattice_constants.c

    # Exercise functionality
    iucr_unit_cell = iucr_conventional_cell(
        UnitCell(lattice_constants, XtallographyUtils.BODY)
    )

    # Check results
    standardized_lattice_constants, _ = standardize(
        expected_orthorhombic_lattice_constants, XtallographyUtils.FACE
    )
    @test iucr_unit_cell.lattice_constants ≈ standardized_lattice_constants
    @test iucr_unit_cell.centering == XtallographyUtils.FACE
end

@testset "iucr_conventional_cell(): limiting cases, centering = BODY, rhombohedral limit cases" begin
    # --- Tests

    # ------ r_α < π/3: a^2 + b^2 = c^2, a^2 + a c cos(β) = b^2

    # Construct basis for rhombohedral unit cell
    r_a = 1.0
    r_α = π / 5

    expected_rhombohedral_lattice_constants = RhombohedralLatticeConstants(r_a, r_α)

    r_basis_a, r_basis_b, r_basis_c = basis(expected_rhombohedral_lattice_constants)

    # Construct lattice constants for monoclinic unit cell
    basis_a = r_basis_a
    basis_b = -r_basis_b + r_basis_c
    basis_c = r_basis_a - r_basis_b - r_basis_c
    lattice_constants = MonoclinicLatticeConstants(
        norm(basis_a),
        norm(basis_b),
        norm(basis_c),
        acos(dot(basis_a, basis_c) / norm(basis_a) / norm(basis_c)),
    )

    # Check test conditions
    @test lattice_constants isa MonoclinicLatticeConstants
    @test lattice_constants.a^2 + lattice_constants.b^2 ≈ lattice_constants.c^2
    @test lattice_constants.a^2 +
          lattice_constants.a * lattice_constants.c * cos(lattice_constants.β) ≈
        lattice_constants.b^2

    # Exercise functionality
    iucr_unit_cell = iucr_conventional_cell(
        UnitCell(lattice_constants, XtallographyUtils.BODY)
    )

    # Check results
    @test iucr_unit_cell.lattice_constants ≈
        standardize(expected_rhombohedral_lattice_constants)
    @test iucr_unit_cell.centering == XtallographyUtils.PRIMITIVE

    # ------ π/3 < r_α < π/2: a^2 + b^2 = c^2, b^2 + a c cos(β) = a^2

    # Construct basis for rhombohedral unit cell
    r_a = 1.0
    r_α = 2π / 5

    expected_rhombohedral_lattice_constants = RhombohedralLatticeConstants(r_a, r_α)

    r_basis_a, r_basis_b, r_basis_c = basis(expected_rhombohedral_lattice_constants)

    # Construct lattice constants for monoclinic unit cell
    basis_a = r_basis_a
    basis_b = r_basis_b - r_basis_c
    basis_c = -r_basis_a + r_basis_b + r_basis_c
    lattice_constants = MonoclinicLatticeConstants(
        norm(basis_a),
        norm(basis_b),
        norm(basis_c),
        acos(dot(basis_a, basis_c) / norm(basis_a) / norm(basis_c)),
    )

    # Check test conditions
    @test lattice_constants isa MonoclinicLatticeConstants
    @test lattice_constants.a^2 + lattice_constants.b^2 ≈ lattice_constants.c^2
    @test lattice_constants.b^2 +
          lattice_constants.a * lattice_constants.c * cos(lattice_constants.β) ≈
        lattice_constants.a^2

    # Exercise functionality
    iucr_unit_cell = iucr_conventional_cell(
        UnitCell(lattice_constants, XtallographyUtils.BODY)
    )

    # Check results
    @test iucr_unit_cell.lattice_constants ≈
        standardize(expected_rhombohedral_lattice_constants)
    @test iucr_unit_cell.centering == XtallographyUtils.PRIMITIVE

    # ------ π/2 < r_α < acos(-1/3):  a^2 + b^2 = c^2, b^2 + a c cos(β) = a^2

    # Construct basis for rhombohedral unit cell
    r_a = 1.0
    r_α = 3π / 5

    expected_rhombohedral_lattice_constants = RhombohedralLatticeConstants(r_a, r_α)

    r_basis_a, r_basis_b, r_basis_c = basis(expected_rhombohedral_lattice_constants)

    # Construct lattice constants for monoclinic unit cell
    basis_a = -r_basis_a
    basis_b = -r_basis_b + r_basis_c
    basis_c = r_basis_a + r_basis_b + r_basis_c
    lattice_constants = MonoclinicLatticeConstants(
        norm(basis_a),
        norm(basis_b),
        norm(basis_c),
        acos(dot(basis_a, basis_c) / norm(basis_a) / norm(basis_c)),
    )

    # Check test conditions
    @test lattice_constants isa MonoclinicLatticeConstants
    @test lattice_constants.c^2 + 3 * lattice_constants.b^2 ≈ 9 * lattice_constants.a^2
    @test lattice_constants.c ≈ -3 * lattice_constants.a * cos(lattice_constants.β)

    # Exercise functionality
    iucr_unit_cell = iucr_conventional_cell(
        UnitCell(lattice_constants, XtallographyUtils.BODY)
    )

    # Check results
    @test iucr_unit_cell.lattice_constants ≈
        standardize(expected_rhombohedral_lattice_constants)
    @test iucr_unit_cell.centering == XtallographyUtils.PRIMITIVE

    # ------ r_α > acos(-1/3): a^2 + 3 b^2 = 9 c^2, a = -3 c cos(β)

    # Construct basis for rhombohedral unit cell
    r_a = 1.0
    r_α = 1.05 * acos(-1 / 3)

    expected_rhombohedral_lattice_constants = RhombohedralLatticeConstants(r_a, r_α)

    r_basis_a, r_basis_b, r_basis_c = basis(expected_rhombohedral_lattice_constants)

    # Construct lattice constants for monoclinic unit cell
    basis_a = r_basis_a + r_basis_b + r_basis_c
    basis_b = r_basis_b - r_basis_c
    basis_c = -r_basis_a
    lattice_constants = MonoclinicLatticeConstants(
        norm(basis_a),
        norm(basis_b),
        norm(basis_c),
        acos(dot(basis_a, basis_c) / norm(basis_a) / norm(basis_c)),
    )

    # Check test conditions
    @test lattice_constants isa MonoclinicLatticeConstants
    @test lattice_constants.a^2 + 3 * lattice_constants.b^2 ≈ 9 * lattice_constants.c^2
    @test lattice_constants.a ≈ -3 * lattice_constants.c * cos(lattice_constants.β)

    # Exercise functionality
    iucr_unit_cell = iucr_conventional_cell(
        UnitCell(lattice_constants, XtallographyUtils.BODY)
    )

    # Check results
    @test iucr_unit_cell.lattice_constants ≈
        standardize(expected_rhombohedral_lattice_constants)
    @test iucr_unit_cell.centering == XtallographyUtils.PRIMITIVE
end

@testset "iucr_conventional_cell(): limiting cases, centering = BODY, non-limit cases" begin
    # --- Preparations

    # Construct basis for orthorhombic unit cell
    o_a = 1.0
    o_b = 2.0
    o_c = 3.0

    expected_orthorhombic_lattice_constants = OrthorhombicLatticeConstants(o_a, o_b, o_c)

    o_basis_a, o_basis_b, o_basis_c = basis(expected_orthorhombic_lattice_constants)

    # --- Tests

    # ------ monoclinic unit cell is not equivalent to an orthorhombic unit cell

    # Construct lattice constants for monoclinic unit cell
    a_ref = 1.0
    b_ref = 2.0
    c_ref = 3.0
    β_ref = 4π / 7

    a = a_ref
    b = b_ref
    c = sqrt((2 * a_ref)^2 + c_ref^2 + 2 * (2 * a_ref) * c_ref * cos(β_ref))
    β = asin(sin(β_ref) / c * c_ref)
    lattice_constants = MonoclinicLatticeConstants(a, b, c, β)

    # Check test conditions
    @test lattice_constants isa MonoclinicLatticeConstants

    # Exercise functionality
    iucr_unit_cell = iucr_conventional_cell(
        UnitCell(lattice_constants, XtallographyUtils.BODY)
    )

    # Check results
    expected_lattice_constants = MonoclinicLatticeConstants(a_ref, b_ref, c_ref, β_ref)

    @test iucr_unit_cell.lattice_constants isa MonoclinicLatticeConstants
    @test iucr_unit_cell.lattice_constants ≈ expected_lattice_constants
    @test iucr_unit_cell.centering == XtallographyUtils.BODY
end

@testset "iucr_conventional_cell(): invalid arguments" begin
    # --- Preparations

    # Construct lattice constants for monoclinic unit cell
    a = 1.0
    b = 2.0
    c = 3.0
    β = 3π / 5
    lattice_constants = MonoclinicLatticeConstants(a, b, c, β)

    # --- Tests

    # centering = FACE
    centering = XtallographyUtils.FACE
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
        "Invalid Bravais lattice: (lattice_system=Monoclinic, centering=$centering)"

    @test startswith(error_message, expected_error)
end

@testset "iucr_conventional_cell(): chain of limiting cases" begin
    # --- Preparations

    a = 5
    b = 10
    c = 7
    β = 5π / 8

    m_lattice_constants = MonoclinicLatticeConstants(a, b, c, β)
    basis_a, basis_b, basis_c = basis(m_lattice_constants)

    # --- Exercise functionality and check results

    # primitive unit cell: aP --> mP
    triclinic_unit_cell = UnitCell(
        LatticeConstants(basis_a, basis_b, basis_c; identify_lattice_system=false),
        XtallographyUtils.PRIMITIVE,
    )
    expected_unit_cell = standardize(
        UnitCell(m_lattice_constants, XtallographyUtils.PRIMITIVE)
    )
    @test triclinic_unit_cell.lattice_constants isa TriclinicLatticeConstants
    @test expected_unit_cell.lattice_constants isa MonoclinicLatticeConstants
    @test iucr_conventional_cell(triclinic_unit_cell) ≈ expected_unit_cell

    # body-centered unit cell: aP --> mI
    #
    # Case #1: m_basis_a and m_basis_b in triclinic basis
    triclinic_unit_cell = UnitCell(
        LatticeConstants(
            basis_a,
            basis_b,
            0.5 * (basis_a + basis_b + basis_c);
            identify_lattice_system=false,
        ),
        XtallographyUtils.PRIMITIVE,
    )
    expected_unit_cell = standardize(UnitCell(m_lattice_constants, XtallographyUtils.BODY))
    @test triclinic_unit_cell.lattice_constants isa TriclinicLatticeConstants
    @test expected_unit_cell.lattice_constants isa MonoclinicLatticeConstants
    @test iucr_conventional_cell(triclinic_unit_cell) ≈ expected_unit_cell

    # body-centered unit cell: aP --> mI
    #
    # Case #1: m_basis_a and m_basis_c in triclinic basis
    triclinic_unit_cell = UnitCell(
        LatticeConstants(
            basis_a,
            basis_c,
            0.5 * (basis_a + basis_b + basis_c);
            identify_lattice_system=false,
        ),
        XtallographyUtils.PRIMITIVE,
    )
    expected_unit_cell = standardize(UnitCell(m_lattice_constants, XtallographyUtils.BODY))
    @test triclinic_unit_cell.lattice_constants isa TriclinicLatticeConstants
    @test expected_unit_cell.lattice_constants isa MonoclinicLatticeConstants
    @test iucr_conventional_cell(triclinic_unit_cell) ≈ expected_unit_cell

    # base-centered unit cell: aP --> mS
    #
    # Case #1: m_basis_a and m_basis_c in triclinic basis
    triclinic_unit_cell = UnitCell(
        LatticeConstants(
            basis_a, 0.5 * (basis_a + basis_b), basis_c; identify_lattice_system=false
        ),
        XtallographyUtils.PRIMITIVE,
    )
    expected_unit_cell = standardize(UnitCell(m_lattice_constants, XtallographyUtils.BASE))
    @test triclinic_unit_cell.lattice_constants isa TriclinicLatticeConstants
    @test expected_unit_cell.lattice_constants isa MonoclinicLatticeConstants
    @test iucr_conventional_cell(triclinic_unit_cell) ≈ expected_unit_cell

    # base-centered unit cell: aP --> mS
    #
    # Case #2: m_basis_a and two base-centered lattice vectors in triclinic basis
    triclinic_unit_cell = UnitCell(
        LatticeConstants(
            basis_a,
            0.5 * (basis_a + basis_b),
            0.5 * (basis_a + basis_b) + basis_c;
            identify_lattice_system=false,
        ),
        XtallographyUtils.PRIMITIVE,
    )
    expected_unit_cell = standardize(UnitCell(m_lattice_constants, XtallographyUtils.BASE))
    @test triclinic_unit_cell.lattice_constants isa TriclinicLatticeConstants
    @test expected_unit_cell.lattice_constants isa MonoclinicLatticeConstants
    @test iucr_conventional_cell(triclinic_unit_cell) ≈ expected_unit_cell
end

# ------ convert_to_body_centering(), convert_to_base_centering()

@testset "convert_to_body_centering()" begin
    # --- Preparations

    a_base = 6
    b_base = 10
    c_base = 8
    β_base = 1.1 * π / 2
    lattice_constants_base_centered = MonoclinicLatticeConstants(
        a_base, b_base, c_base, β_base
    )

    # --- Exercise functionality and check results

    # Compute expected body-centered lattice constants
    a_body = sqrt(a_base^2 + c_base^2 - 2 * a_base * c_base * cos(β_base))
    b_body = b_base
    c_body = c_base
    β_body = asin(sin(β_base) / a_body * a_base)
    expected_lattice_constants, _ = standardize(
        MonoclinicLatticeConstants(a_body, b_body, c_body, β_body), XtallographyUtils.BODY
    )

    @test convert_to_body_centering(lattice_constants_base_centered) ≈
        expected_lattice_constants
end

@testset "convert_to_base_centering()" begin
    # --- Preparations

    a_body = 6
    b_body = 10
    c_body = 8
    β_body = 1.1 * π / 2
    lattice_constants_body_centered = MonoclinicLatticeConstants(
        a_body, b_body, c_body, β_body
    )

    # --- Exercise functionality and check results

    # Compute expected base-centered lattice constants
    a_base = sqrt(a_body^2 + c_body^2 + 2 * a_body * c_body * cos(β_body))
    b_base = b_body
    c_base = a_body
    β_base = π - asin(sin(β_body) / a_base * c_body)
    expected_lattice_constants = MonoclinicLatticeConstants(a_base, b_base, c_base, β_base)

    @test convert_to_base_centering(lattice_constants_body_centered) ≈
        expected_lattice_constants
end

# ------ reduced_cell()

@testset "reduced_cell()" begin
    # --- Preparations

    a = 2
    b = 3
    c = 5
    β = 4π / 7
    lattice_constants = MonoclinicLatticeConstants(a, b, c, β)
    basis_a, basis_b, basis_c = basis(lattice_constants)

    # --- Exercise functionality and check results

    # primitive unit cell defined by [basis_a, basis_b, basis_c]
    unit_cell = UnitCell(lattice_constants, XtallographyUtils.PRIMITIVE)

    expected_reduced_cell = reduced_cell(
        UnitCell(lattice_constants, XtallographyUtils.PRIMITIVE)
    )

    reduced_cell_ = reduced_cell(unit_cell)
    @test reduced_cell_.lattice_constants isa MonoclinicLatticeConstants
    @test volume(reduced_cell_) ≈ volume(unit_cell)
    @test reduced_cell_ ≈ expected_reduced_cell

    # primitive unit cell defined by linear combination of [basis_a, basis_b, basis_c]
    unit_cell = UnitCell(
        LatticeConstants(basis_a + basis_c, basis_b, basis_c), XtallographyUtils.PRIMITIVE
    )

    expected_reduced_cell = reduced_cell(
        UnitCell(lattice_constants, XtallographyUtils.PRIMITIVE)
    )

    reduced_cell_ = reduced_cell(unit_cell)
    @test reduced_cell_.lattice_constants isa MonoclinicLatticeConstants
    @test volume(reduced_cell_) ≈ volume(unit_cell)
    @test reduced_cell_ ≈ expected_reduced_cell

    # body-centered unit cell
    unit_cell = UnitCell(lattice_constants, XtallographyUtils.BODY)

    expected_reduced_cell = reduced_cell(
        UnitCell(
            LatticeConstants(basis_a, basis_b, 0.5 * (basis_a + basis_b + basis_c)),
            XtallographyUtils.PRIMITIVE,
        ),
    )

    reduced_cell_ = reduced_cell(unit_cell)
    @test reduced_cell_.lattice_constants isa TriclinicLatticeConstants
    @test volume(reduced_cell_) ≈ 0.5 * volume(unit_cell)
    @test reduced_cell_ ≈ expected_reduced_cell

    # base-centered unit cell
    unit_cell = UnitCell(lattice_constants, XtallographyUtils.BASE)

    expected_reduced_cell = reduced_cell(
        UnitCell(
            LatticeConstants(basis_a, 0.5 * (basis_a + basis_b), basis_c),
            XtallographyUtils.PRIMITIVE,
        ),
    )

    reduced_cell_ = reduced_cell(unit_cell)
    @test reduced_cell_.lattice_constants isa TriclinicLatticeConstants
    @test volume(reduced_cell_) ≈ 0.5 * volume(unit_cell)
    @test reduced_cell_ ≈ expected_reduced_cell

    # equivalent body-centered and base-centered monoclinic unit cells
    base_centered_unit_cell = UnitCell(lattice_constants, XtallographyUtils.BASE)
    body_centered_basis_a = basis_a + basis_c
    body_centered_basis_b = basis_b
    body_centered_basis_c = basis_c
    body_centered_unit_cell = UnitCell(
        MonoclinicLatticeConstants(
            norm(body_centered_basis_a),
            norm(body_centered_basis_b),
            norm(body_centered_basis_c),
            π - acos(
                dot(body_centered_basis_a, body_centered_basis_c) /
                norm(body_centered_basis_a) / norm(body_centered_basis_c),
            ),
        ),
        XtallographyUtils.BODY,
    )

    @test base_centered_unit_cell.lattice_constants isa MonoclinicLatticeConstants
    @test body_centered_unit_cell.lattice_constants isa MonoclinicLatticeConstants

    @test reduced_cell(base_centered_unit_cell) ≈ reduced_cell(body_centered_unit_cell)

    # face-centered unit cell (equivalent to smaller body-centered unit cell)
    face_centered_unit_cell = UnitCell(lattice_constants, XtallographyUtils.FACE)
    body_centered_basis_a = 0.5 * (basis_a + basis_c)
    body_centered_basis_b = basis_b
    body_centered_basis_c = 0.5 * (basis_a - basis_c)
    body_centered_unit_cell = UnitCell(
        MonoclinicLatticeConstants(
            norm(body_centered_basis_a),
            norm(body_centered_basis_b),
            norm(body_centered_basis_c),
            acos(
                dot(body_centered_basis_a, body_centered_basis_c) /
                norm(body_centered_basis_a) / norm(body_centered_basis_c),
            ),
        ),
        XtallographyUtils.BODY,
    )

    reduced_face_centered_unit_cell = reduced_cell(face_centered_unit_cell)
    reduced_body_centered_unit_cell = reduced_cell(body_centered_unit_cell)

    @test face_centered_unit_cell.lattice_constants isa MonoclinicLatticeConstants
    @test body_centered_unit_cell.lattice_constants isa MonoclinicLatticeConstants

    @test reduced_face_centered_unit_cell ≈ reduced_body_centered_unit_cell
    @test volume(reduced_body_centered_unit_cell) ≈ 0.25 * volume(face_centered_unit_cell)
    @test volume(reduced_face_centered_unit_cell) ≈ 0.25 * volume(face_centered_unit_cell)
end
