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
using LinearAlgebra: norm
using Logging
using Test

# Xtallography package
using Xtallography

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

@testset "conventional_cell(::OrthorhombicUnitCell): non-limiting cases" begin
    # orthorhombic unit cells that are not equivalent to any higher symmetry unit cell

    # --- centering = primitive

    # Construct orthorhombic unit cell
    a = 2.0
    b = 1.0
    c = 3.0
    unit_cell = OrthorhombicUnitCell(a, b, c; centering=primitive_centering)

    # Exercise functionality
    iucr_unit_cell = conventional_cell(unit_cell)

    # Check results
    @test iucr_unit_cell ≈ standardize(unit_cell)

    # --- centering = body-centering

    # Construct orthorhombic unit cell
    a = 2.0
    b = 1.0
    c = 3.0
    unit_cell = OrthorhombicUnitCell(a, b, c; centering=body_centering)

    # Exercise functionality
    iucr_unit_cell = conventional_cell(unit_cell)

    # Check results
    @test iucr_unit_cell ≈ standardize(unit_cell)

    # --- centering = face_centering

    # Construct orthorhombic unit cell
    a = 2.0
    b = 1.0
    c = 3.0
    unit_cell = OrthorhombicUnitCell(a, b, c; centering=face_centering)

    # Exercise functionality
    iucr_unit_cell = conventional_cell(unit_cell)

    # Check results
    @test iucr_unit_cell ≈ standardize(unit_cell)

    # --- centering = base_centering

    # Construct orthorhombic unit cell
    a = 2.0
    b = 1.0
    c = 3.0
    unit_cell = OrthorhombicUnitCell(a, b, c; centering=base_centering)

    # Exercise functionality
    iucr_unit_cell = conventional_cell(unit_cell)

    # Check results
    @test iucr_unit_cell ≈ standardize(unit_cell)
end

@testset "conventional_cell(::OrthorhombicUnitCell): limiting cases - oP --> tP" begin
    # --- Preparations

    # Construct tetragonal unit cell
    t_a = 1.0
    t_c = 3.0

    # --- Tests

    # ------ a = b

    # Construct orthorhombic unit cell
    a = t_a
    b = t_a
    c = t_c
    unit_cell = OrthorhombicUnitCell(a, b, c; centering=primitive_centering)

    # Check test conditions
    @test unit_cell isa OrthorhombicUnitCell
    @test centering(unit_cell) === primitive_centering
    @test lattice_constants(unit_cell).a ≈ lattice_constants(unit_cell).b

    # Exercise functionality
    iucr_unit_cell = conventional_cell(unit_cell)

    # Check results
    @test iucr_unit_cell ≈ TetragonalUnitCell(t_a, t_c; centering=primitive_centering)

    # ------ b = c

    # Construct orthorhombic unit cell
    a = t_a
    b = t_c
    c = t_c
    unit_cell = OrthorhombicUnitCell(a, b, c; centering=primitive_centering)

    # Check test conditions
    @test unit_cell isa OrthorhombicUnitCell
    @test centering(unit_cell) === primitive_centering
    @test lattice_constants(unit_cell).b ≈ lattice_constants(unit_cell).c

    # Exercise functionality
    iucr_unit_cell = conventional_cell(unit_cell)

    # Check results
    @test iucr_unit_cell ≈ TetragonalUnitCell(t_c, t_a; centering=primitive_centering)
end

@testset "conventional_cell(::OrthorhombicUnitCell): limiting cases - oI --> tI" begin
    # --- Preparations

    # Construct tetragonal unit cell
    t_a = 1.0
    t_c = 3.0

    # --- Tests

    # ------ a = b

    # Construct orthorhombic unit cell
    a = t_a
    b = t_a
    c = t_c
    unit_cell = OrthorhombicUnitCell(a, b, c; centering=body_centering)

    # Check test conditions
    @test unit_cell isa OrthorhombicUnitCell
    @test centering(unit_cell) === body_centering
    @test lattice_constants(unit_cell).a ≈ lattice_constants(unit_cell).b

    # Exercise functionality
    iucr_unit_cell = conventional_cell(unit_cell)

    # Check results
    expected_unit_cell = TetragonalUnitCell(t_a, t_c; centering=body_centering)
    @test iucr_unit_cell ≈ expected_unit_cell

    # ------ b = c

    # Construct orthorhombic unit cell
    a = t_a
    b = t_c
    c = t_c
    unit_cell = OrthorhombicUnitCell(a, b, c; centering=body_centering)

    # Check test conditions
    @test unit_cell isa OrthorhombicUnitCell
    @test centering(unit_cell) === body_centering
    @test lattice_constants(unit_cell).b ≈ lattice_constants(unit_cell).c

    # Exercise functionality
    iucr_unit_cell = conventional_cell(unit_cell)

    # Check results
    expected_unit_cell = TetragonalUnitCell(t_c, t_a; centering=body_centering)
    @test iucr_unit_cell ≈ expected_unit_cell
end

@testset "conventional_cell(::OrthorhombicUnitCell): limiting cases - oF --> tI" begin
    # --- Tests

    # ------ a = b

    # Construct basis for tetragonal unit cell
    t_a = 1.0
    t_c = 3.0
    t_basis_a, t_basis_b, t_basis_c = basis(TetragonalUnitCell(t_a, t_c))

    # Construct orthorhombic unit cell
    basis_a = t_basis_a - t_basis_b
    basis_b = t_basis_a + t_basis_b
    basis_c = t_basis_c
    unit_cell = OrthorhombicUnitCell(
        norm(basis_a), norm(basis_b), norm(basis_c); centering=face_centering
    )

    # Check test conditions
    @test unit_cell isa OrthorhombicUnitCell
    @test centering(unit_cell) === face_centering
    @test lattice_constants(unit_cell).a ≈ lattice_constants(unit_cell).b

    # Exercise functionality
    iucr_unit_cell = conventional_cell(unit_cell)

    # Check results
    expected_unit_cell = TetragonalUnitCell(t_a, t_c; centering=body_centering)
    @test iucr_unit_cell ≈ expected_unit_cell

    # ------ b = c

    # Construct basis for tetragonal unit cell
    t_a = 3.0
    t_c = 1.0
    t_basis_a, t_basis_b, t_basis_c = basis(TetragonalUnitCell(t_a, t_c))

    # Construct orthorhombic unit cell
    basis_a = t_basis_c
    basis_b = t_basis_a - t_basis_b
    basis_c = t_basis_a + t_basis_b
    unit_cell = OrthorhombicUnitCell(
        norm(basis_a), norm(basis_b), norm(basis_c); centering=face_centering
    )

    # Check test conditions
    @test unit_cell isa OrthorhombicUnitCell
    @test centering(unit_cell) === face_centering
    @test lattice_constants(unit_cell).b ≈ lattice_constants(unit_cell).c

    # Exercise functionality
    iucr_unit_cell = conventional_cell(unit_cell)

    # Check results
    expected_unit_cell = TetragonalUnitCell(t_a, t_c; centering=body_centering)
    @test iucr_unit_cell ≈ expected_unit_cell
end

@testset "conventional_cell(::OrthorhombicUnitCell): limiting cases - oC --> tP,hP" begin
    # --- Preparations

    # Construct basis for tetragonal unit cell
    t_a = 1.0
    t_c = 3.0
    t_basis_a, t_basis_b, t_basis_c = basis(TetragonalUnitCell(t_a, t_c))

    # --- Tests

    # ------ a = b < c

    # Construct orthorhombic unit cell
    a = 1.0
    b = a
    c = 2 * a
    unit_cell = standardize(OrthorhombicUnitCell(a, b, c; centering=base_centering))

    # Check test conditions
    @test unit_cell isa OrthorhombicUnitCell
    @test centering(unit_cell) === base_centering
    @test lattice_constants(unit_cell).a ≈ lattice_constants(unit_cell).b

    # Exercise functionality
    iucr_unit_cell = conventional_cell(unit_cell)

    # Check results
    @test iucr_unit_cell ≈ TetragonalUnitCell(a / sqrt(2), c; centering=primitive_centering)

    # ------ a = b > c

    # Construct orthorhombic unit cell
    a = 1.0
    b = a
    c = 0.5 * a
    unit_cell = standardize(OrthorhombicUnitCell(a, b, c; centering=base_centering))

    # Check test conditions
    @test unit_cell isa OrthorhombicUnitCell
    @test centering(unit_cell) === base_centering
    @test lattice_constants(unit_cell).a ≈ lattice_constants(unit_cell).b

    # Exercise functionality
    iucr_unit_cell = conventional_cell(unit_cell)

    # Check results
    @test iucr_unit_cell ≈ TetragonalUnitCell(a / sqrt(2), c; centering=primitive_centering)

    # ------ b = a * sqrt(3)

    # Construct basis for hexagonal unit cell
    h_a = 3.0
    h_c = 1.0
    h_basis_a, h_basis_b, h_basis_c = basis(HexagonalUnitCell(h_a, h_c))

    # Construct orthorhombic unit cell
    basis_a = h_basis_a
    basis_b = h_basis_a + 2 * h_basis_b
    basis_c = h_basis_c
    unit_cell = OrthorhombicUnitCell(
        norm(basis_a), norm(basis_b), norm(basis_c); centering=base_centering
    )

    # Check test conditions
    @test unit_cell isa OrthorhombicUnitCell
    @test centering(unit_cell) === base_centering
    @test lattice_constants(unit_cell).b ≈ lattice_constants(unit_cell).a * sqrt(3)

    # Exercise functionality
    iucr_unit_cell = conventional_cell(unit_cell)

    # Check results
    @test iucr_unit_cell ≈ HexagonalUnitCell(h_a, h_c)
end

@testset "conventional_cell()::OrthorhombicUnitCell: chain of limiting cases" begin
    # --- Preparations

    a = 5
    b = 7
    c = 9

    orthorhombic_unit_cell = OrthorhombicUnitCell(a, b, c)
    basis_a, basis_b, basis_c = basis(orthorhombic_unit_cell)

    # --- Exercise functionality and check results

    # primitive unit cell: aP --> mP --> oP
    triclinic_unit_cell = UnitCell(
        basis_a,
        basis_b,
        basis_c;
        identify_lattice_system=false,
        centering=primitive_centering,
    )
    expected_unit_cell = standardize(orthorhombic_unit_cell)
    @test triclinic_unit_cell isa TriclinicUnitCell
    @test expected_unit_cell isa OrthorhombicUnitCell
    @test centering(expected_unit_cell) === primitive_centering
    @debug "chain of limiting cases: aP --> mP --> oP"
    @test conventional_cell(triclinic_unit_cell) ≈ expected_unit_cell

    # body-centered unit cell: aP --> mI --> oI
    triclinic_unit_cell = UnitCell(
        basis_a,
        basis_b,
        0.5 * (basis_a + basis_b + basis_c);
        identify_lattice_system=false,
        centering=primitive_centering,
    )
    expected_unit_cell = standardize(
        OrthorhombicUnitCell(a, b, c; centering=body_centering)
    )
    @test triclinic_unit_cell isa TriclinicUnitCell
    @test expected_unit_cell isa OrthorhombicUnitCell
    @test centering(expected_unit_cell) === body_centering
    @debug "chain of limiting cases: aP --> mI --> oI"
    @test conventional_cell(triclinic_unit_cell) ≈ expected_unit_cell

    # face-centered unit cell: aP --> mI --> oF
    triclinic_unit_cell = UnitCell(
        0.5 * (basis_a + basis_b),
        0.5 * (basis_a - basis_b),
        0.5 * (basis_a - basis_c);
        identify_lattice_system=false,
        centering=primitive_centering,
    )
    expected_unit_cell = standardize(
        OrthorhombicUnitCell(a, b, c; centering=face_centering)
    )
    @test triclinic_unit_cell isa TriclinicUnitCell
    @test expected_unit_cell isa OrthorhombicUnitCell
    @test centering(expected_unit_cell) === face_centering
    @debug "chain of limiting cases: aP --> mI --> oF"
    @test conventional_cell(triclinic_unit_cell) ≈ expected_unit_cell

    # base-centered unit cell: aP --> mP --> oC
    triclinic_unit_cell = UnitCell(
        basis_a,
        0.5 * (basis_a + basis_b),
        basis_c;
        identify_lattice_system=false,
        centering=primitive_centering,
    )
    expected_unit_cell = standardize(
        OrthorhombicUnitCell(a, b, c; centering=base_centering)
    )
    @test triclinic_unit_cell isa TriclinicUnitCell
    @test expected_unit_cell isa OrthorhombicUnitCell
    @test centering(expected_unit_cell) === base_centering
    @debug "chain of limiting cases: aP --> mP --> oC"
    @test conventional_cell(triclinic_unit_cell) ≈ expected_unit_cell

    # base-centered unit cell: aP --> mI --> oC
    triclinic_unit_cell = UnitCell(
        basis_a,
        0.5 * (basis_a + basis_b),
        0.5 * (basis_a + basis_b) + basis_c;
        identify_lattice_system=false,
        centering=primitive_centering,
    )
    expected_unit_cell = standardize(
        OrthorhombicUnitCell(a, b, c; centering=base_centering)
    )
    @test triclinic_unit_cell isa TriclinicUnitCell
    @test expected_unit_cell isa OrthorhombicUnitCell
    @test centering(expected_unit_cell) === base_centering
    @debug "chain of limiting cases: aP --> mI --> oC"
    @test conventional_cell(triclinic_unit_cell) ≈ expected_unit_cell
end
