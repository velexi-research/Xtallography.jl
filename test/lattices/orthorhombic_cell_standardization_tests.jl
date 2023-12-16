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
Tests for the unit cell standardization methods for orthorhombic lattices
"""
# --- Imports

# Standard library
using Test
using LinearAlgebra: norm

# XtallographyUtils package
using XtallographyUtils

# Notes
# =====
# These tests adopt the following variable conventions.
#
# - Unless otherwise noted, lattice constants and basis vectors refer to the orthorhombic
#   (not tetragonal or hexagonal) unit cell.
#
# - Lattice constants and basis vectors for tetragonal unit cells are indicated by
#   the "t_" prefix.
#
# - Lattice constants and basis vectors for hexagonal unit cells are indicated by
#   the "h_" prefix.

# --- Tests

@testset "iucr_conventional_cell(): limiting cases, centering = PRIMITIVE" begin
    # --- Preparations

    # Construct lattice constants for tetragonal unit cell
    t_a = 1.0
    t_c = 3.0

    # --- Tests

    # ------ a = b

    # Construct lattice constants for orthorhombic unit cell
    a = t_a
    b = t_a
    c = t_c
    lattice_constants = OrthorhombicLatticeConstants(a, b, c)

    # Check test conditions
    @test lattice_constants isa OrthorhombicLatticeConstants
    @test lattice_constants.a ≈ lattice_constants.b

    # Exercise functionality
    iucr_unit_cell = iucr_conventional_cell(
        UnitCell(lattice_constants, XtallographyUtils.PRIMITIVE)
    )

    # Check results
    @test iucr_unit_cell.lattice_constants ≈ TetragonalLatticeConstants(t_a, t_c)
    @test iucr_unit_cell.centering == XtallographyUtils.PRIMITIVE

    # ------ b = c

    # Construct lattice constants for monoclinic unit cell
    a = t_a
    b = t_c
    c = t_c
    lattice_constants = OrthorhombicLatticeConstants(a, b, c)

    # Check test conditions
    @test lattice_constants isa OrthorhombicLatticeConstants
    @test lattice_constants.b ≈ lattice_constants.c

    # Exercise functionality
    iucr_unit_cell = iucr_conventional_cell(
        UnitCell(lattice_constants, XtallographyUtils.PRIMITIVE)
    )

    # Check results
    @test iucr_unit_cell.lattice_constants ≈ TetragonalLatticeConstants(t_c, t_a)
    @test iucr_unit_cell.centering == XtallographyUtils.PRIMITIVE

    # ------ orthorhombic unit cell is not equivalent to a tetragonal unit cell

    # Construct lattice constants for orthorhombic unit cell
    a = 2.0
    b = 1.0
    c = 3.0
    lattice_constants = OrthorhombicLatticeConstants(a, b, c)

    # Exercise functionality
    iucr_unit_cell = iucr_conventional_cell(
        UnitCell(lattice_constants, XtallographyUtils.PRIMITIVE)
    )

    # Check results
    @test iucr_unit_cell.lattice_constants ≈ standardize(lattice_constants)
    @test iucr_unit_cell.centering == XtallographyUtils.PRIMITIVE
end

@testset "iucr_conventional_cell(): limiting cases, centering = BODY" begin
    # --- Preparations

    # Construct lattice constants for tetragonal unit cell
    t_a = 1.0
    t_c = 3.0

    # --- Tests

    # ------ a = b

    # Construct lattice constants for orthorhombic unit cell
    a = t_a
    b = t_a
    c = t_c
    lattice_constants = OrthorhombicLatticeConstants(a, b, c)

    # Check test conditions
    @test lattice_constants isa OrthorhombicLatticeConstants
    @test lattice_constants.a ≈ lattice_constants.b

    # Exercise functionality
    iucr_unit_cell = iucr_conventional_cell(
        UnitCell(lattice_constants, XtallographyUtils.BODY)
    )

    # Check results
    @test iucr_unit_cell.lattice_constants ≈ TetragonalLatticeConstants(t_a, t_c)
    @test iucr_unit_cell.centering == XtallographyUtils.BODY

    # ------ b = c

    # Construct lattice constants for monoclinic unit cell
    a = t_a
    b = t_c
    c = t_c
    lattice_constants = OrthorhombicLatticeConstants(a, b, c)

    # Check test conditions
    @test lattice_constants isa OrthorhombicLatticeConstants
    @test lattice_constants.b ≈ lattice_constants.c

    # Exercise functionality
    iucr_unit_cell = iucr_conventional_cell(
        UnitCell(lattice_constants, XtallographyUtils.BODY)
    )

    # Check results
    @test iucr_unit_cell.lattice_constants ≈ TetragonalLatticeConstants(t_c, t_a)
    @test iucr_unit_cell.centering == XtallographyUtils.BODY

    # ------ orthorhombic unit cell is not equivalent to a tetragonal unit cell

    # Construct lattice constants for orthorhombic unit cell
    a = 2.0
    b = 1.0
    c = 3.0
    lattice_constants = OrthorhombicLatticeConstants(a, b, c)

    # Exercise functionality
    iucr_unit_cell = iucr_conventional_cell(
        UnitCell(lattice_constants, XtallographyUtils.BODY)
    )

    # Check results
    expected_lattice_constants, _ = standardize(lattice_constants, XtallographyUtils.BODY)
    @test iucr_unit_cell.lattice_constants ≈ expected_lattice_constants
    @test iucr_unit_cell.centering == XtallographyUtils.BODY
end

@testset "iucr_conventional_cell(): limiting cases, centering = FACE" begin
    # --- Tests

    # ------ a = b

    # Construct basis for tetragonal unit cell
    t_a = 1.0
    t_c = 3.0
    t_basis_a, t_basis_b, t_basis_c = basis(TetragonalLatticeConstants(t_a, t_c))

    # Construct lattice constants for orthorhombic unit cell
    basis_a = t_basis_a - t_basis_b
    basis_b = t_basis_a + t_basis_b
    basis_c = t_basis_c
    lattice_constants = OrthorhombicLatticeConstants(
        norm(basis_a), norm(basis_b), norm(basis_c)
    )

    # Check test conditions
    @test lattice_constants isa OrthorhombicLatticeConstants
    @test lattice_constants.a ≈ lattice_constants.b

    # Exercise functionality
    iucr_unit_cell = iucr_conventional_cell(
        UnitCell(lattice_constants, XtallographyUtils.FACE)
    )

    # Check results
    @test iucr_unit_cell.lattice_constants ≈ TetragonalLatticeConstants(t_a, t_c)
    @test iucr_unit_cell.centering == XtallographyUtils.BODY

    # ------ b = c

    # Construct basis for tetragonal unit cell
    t_a = 3.0
    t_c = 1.0
    t_basis_a, t_basis_b, t_basis_c = basis(TetragonalLatticeConstants(t_a, t_c))

    # Construct lattice constants for monoclinic unit cell
    basis_a = t_basis_c
    basis_b = t_basis_a - t_basis_b
    basis_c = t_basis_a + t_basis_b
    lattice_constants = OrthorhombicLatticeConstants(
        norm(basis_a), norm(basis_b), norm(basis_c)
    )

    # Check test conditions
    @test lattice_constants isa OrthorhombicLatticeConstants
    @test lattice_constants.b ≈ lattice_constants.c

    # Exercise functionality
    iucr_unit_cell = iucr_conventional_cell(
        UnitCell(lattice_constants, XtallographyUtils.FACE)
    )

    # Check results
    @test iucr_unit_cell.lattice_constants ≈ TetragonalLatticeConstants(t_a, t_c)
    @test iucr_unit_cell.centering == XtallographyUtils.BODY

    # ------ orthorhombic unit cell is not equivalent to a tetragonal unit cell

    # Construct lattice constants for orthorhombic unit cell
    a = 2.0
    b = 1.0
    c = 3.0
    lattice_constants = OrthorhombicLatticeConstants(a, b, c)

    # Exercise functionality
    iucr_unit_cell = iucr_conventional_cell(
        UnitCell(lattice_constants, XtallographyUtils.FACE)
    )

    # Check results
    expected_lattice_constants, _ = standardize(lattice_constants, XtallographyUtils.FACE)
    @test iucr_unit_cell.lattice_constants ≈ expected_lattice_constants
    @test iucr_unit_cell.centering == XtallographyUtils.FACE
end

@testset "iucr_conventional_cell(): limiting cases, centering = BASE" begin
    # --- Preparations

    # Construct basis for tetragonal unit cell
    t_a = 1.0
    t_c = 3.0
    t_basis_a, t_basis_b, t_basis_c = basis(TetragonalLatticeConstants(t_a, t_c))

    # --- Tests

    # ------ a = b < c

    # Construct lattice constants for orthorhombic unit cell
    a = 1.0
    b = a
    c = 2 * a
    lattice_constants, _ = standardize(
        OrthorhombicLatticeConstants(a, b, c), XtallographyUtils.BASE
    )

    # Check test conditions
    @test lattice_constants isa OrthorhombicLatticeConstants
    @test lattice_constants.a ≈ lattice_constants.b

    # Exercise functionality
    iucr_unit_cell = iucr_conventional_cell(
        UnitCell(lattice_constants, XtallographyUtils.BASE)
    )

    # Check results
    @test iucr_unit_cell.lattice_constants ≈ TetragonalLatticeConstants(a / sqrt(2), c)
    @test iucr_unit_cell.centering == XtallographyUtils.PRIMITIVE

    # ------ a = b > c

    # Construct lattice constants for orthorhombic unit cell
    a = 1.0
    b = a
    c = 0.5 * a
    lattice_constants, _ = standardize(
        OrthorhombicLatticeConstants(a, b, c), XtallographyUtils.BASE
    )

    # Check test conditions
    @test lattice_constants isa OrthorhombicLatticeConstants
    @test lattice_constants.a ≈ lattice_constants.b

    # Exercise functionality
    iucr_unit_cell = iucr_conventional_cell(
        UnitCell(lattice_constants, XtallographyUtils.BASE)
    )

    # Check results
    @test iucr_unit_cell.lattice_constants ≈ TetragonalLatticeConstants(a / sqrt(2), c)
    @test iucr_unit_cell.centering == XtallographyUtils.PRIMITIVE

    # ------ b = a * sqrt(3)

    # Construct basis for hexagonal unit cell
    h_a = 3.0
    h_c = 1.0
    h_basis_a, h_basis_b, h_basis_c = basis(HexagonalLatticeConstants(h_a, h_c))

    # Construct lattice constants for monoclinic unit cell
    basis_a = h_basis_a
    basis_b = h_basis_a + 2 * h_basis_b
    basis_c = h_basis_c
    lattice_constants = OrthorhombicLatticeConstants(
        norm(basis_a), norm(basis_b), norm(basis_c)
    )

    # Check test conditions
    @test lattice_constants isa OrthorhombicLatticeConstants
    @test lattice_constants.b ≈ lattice_constants.a * sqrt(3)

    # Exercise functionality
    iucr_unit_cell = iucr_conventional_cell(
        UnitCell(lattice_constants, XtallographyUtils.BASE)
    )

    # Check results
    @test iucr_unit_cell.lattice_constants ≈ HexagonalLatticeConstants(h_a, h_c)
    @test iucr_unit_cell.centering == XtallographyUtils.PRIMITIVE

    # ------ orthorhombic unit cell is not equivalent to a tetragonal unit cell

    # Construct lattice constants for orthorhombic unit cell
    a = 2.0
    b = 1.0
    c = 3.0
    lattice_constants = OrthorhombicLatticeConstants(a, b, c)

    # Exercise functionality
    iucr_unit_cell = iucr_conventional_cell(
        UnitCell(lattice_constants, XtallographyUtils.BASE)
    )

    # Check results
    expected_lattice_constants, _ = standardize(lattice_constants, XtallographyUtils.BASE)
    @test iucr_unit_cell.lattice_constants ≈ expected_lattice_constants
    @test iucr_unit_cell.centering == XtallographyUtils.BASE
end

@testset "iucr_conventional_cell(): chain of limiting cases" begin
    # --- Preparations

    a = 5
    b = 7
    c = 9

    lattice_constants = OrthorhombicLatticeConstants(a, b, c)
    basis_a, basis_b, basis_c = basis(lattice_constants)

    # --- Exercise functionality and check results

    # primitive unit cell: aP --> mP --> oP
    triclinic_unit_cell = UnitCell(
        LatticeConstants(basis_a, basis_b, basis_c; identify_lattice_system=false),
        XtallographyUtils.PRIMITIVE,
    )
    expected_unit_cell = standardize(
        UnitCell(lattice_constants, XtallographyUtils.PRIMITIVE)
    )
    @test triclinic_unit_cell.lattice_constants isa TriclinicLatticeConstants
    @test expected_unit_cell.lattice_constants isa OrthorhombicLatticeConstants
    @test iucr_conventional_cell(triclinic_unit_cell) ≈ expected_unit_cell

    # body-centered unit cell: aP --> mI --> oI
    triclinic_unit_cell = UnitCell(
        LatticeConstants(
            basis_a,
            basis_b,
            0.5 * (basis_a + basis_b + basis_c);
            identify_lattice_system=false,
        ),
        XtallographyUtils.PRIMITIVE,
    )
    expected_unit_cell = standardize(UnitCell(lattice_constants, XtallographyUtils.BODY))
    @test triclinic_unit_cell.lattice_constants isa TriclinicLatticeConstants
    @test expected_unit_cell.lattice_constants isa OrthorhombicLatticeConstants
    @test iucr_conventional_cell(triclinic_unit_cell) ≈ expected_unit_cell

    # face-centered unit cell: aP --> mI --> oF
    triclinic_unit_cell = UnitCell(
        LatticeConstants(
            0.5 * (basis_a + basis_b),
            0.5 * (basis_a - basis_b),
            0.5 * (basis_a - basis_c);
            identify_lattice_system=false,
        ),
        XtallographyUtils.PRIMITIVE,
    )
    expected_unit_cell = standardize(UnitCell(lattice_constants, XtallographyUtils.FACE))
    @test triclinic_unit_cell.lattice_constants isa TriclinicLatticeConstants
    @test expected_unit_cell.lattice_constants isa OrthorhombicLatticeConstants
    @test iucr_conventional_cell(triclinic_unit_cell) ≈ expected_unit_cell

    # base-centered unit cell: aP --> mP --> oC
    triclinic_unit_cell = UnitCell(
        LatticeConstants(
            basis_a, 0.5 * (basis_a + basis_b), basis_c; identify_lattice_system=false
        ),
        XtallographyUtils.PRIMITIVE,
    )
    expected_unit_cell = standardize(UnitCell(lattice_constants, XtallographyUtils.BASE))
    @test triclinic_unit_cell.lattice_constants isa TriclinicLatticeConstants
    @test expected_unit_cell.lattice_constants isa OrthorhombicLatticeConstants
    @test iucr_conventional_cell(triclinic_unit_cell) ≈ expected_unit_cell

    # base-centered unit cell: aP --> mI --> oC
    triclinic_unit_cell = UnitCell(
        LatticeConstants(
            basis_a,
            0.5 * (basis_a + basis_b),
            0.5 * (basis_a + basis_b) + basis_c;
            identify_lattice_system=false,
        ),
        XtallographyUtils.PRIMITIVE,
    )
    expected_unit_cell = standardize(UnitCell(lattice_constants, XtallographyUtils.BASE))
    @test triclinic_unit_cell.lattice_constants isa TriclinicLatticeConstants
    @test expected_unit_cell.lattice_constants isa OrthorhombicLatticeConstants
    @test iucr_conventional_cell(triclinic_unit_cell) ≈ expected_unit_cell
end
