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
Tests for the unit cell standardization methods for rhombohedral lattices
"""
# --- Imports

# Standard library
using Logging
using Test

# Xtallography package
using Xtallography

# Notes
# =====
# These tests adopt the following variable conventions.
#
# - Unless otherwise noted, lattice constants and basis vectors refer to the rhombohedral
#   (not cubic) unit cell.
#
# - Lattice constants and basis vectors for cubic unit cells are indicated by the "c_"
#   prefix.

# --- Tests

@testset "conventional_cell(::RhombohedralUnitCell): invalid arguments" begin
    # --- Preparations

    # Construct rhombohedral unit cell
    a = 1.0
    α = 2π / 5

    # --- Tests

    # ------ Invalid centering

    for centering_ in (body_centering, face_centering, base_centering)
        expected_message =
            "Invalid Bravais lattice: " *
            "(lattice_system=Rhombohedral, centering=$(nameof(typeof(centering_))))"

        @test_throws ArgumentError(expected_message) conventional_cell(
            RhombohedralUnitCell(a, α; centering=centering_)
        )
    end
end

@testset "conventional_cell(::RhombohedralUnitCell): non-limiting cases" begin
    # rhombohedral unit cells that are not equivalent to any higher symmetry unit cell

    # ------ rhombohedral unit cell is not equivalent to an cubic unit cell

    # Construct rhombohedral unit cell
    a = 1.0
    α = 2π / 5
    unit_cell = RhombohedralUnitCell(a, α)

    # Exercise functionality
    iucr_unit_cell = conventional_cell(unit_cell)

    # Check results
    @test iucr_unit_cell == standardize(unit_cell)
end

@testset "conventional_cell(::RhombohedralUnitCell): limiting cases" begin
    # --- Preparations

    # Construct basis for cubic unit cell
    c_a = 1.0

    primitive_cubic_unit_cell = CubicUnitCell(c_a)
    c_basis_a, c_basis_b, c_basis_c = basis(primitive_cubic_unit_cell)

    # --- Tests

    # ------ α = π/3

    # Construct rhombohedral unit cell
    basis_a = 0.5 * (c_basis_b + c_basis_c)
    basis_b = 0.5 * (c_basis_a + c_basis_c)
    basis_c = 0.5 * (c_basis_a + c_basis_b)
    unit_cell = standardize(
        UnitCell(basis_a, basis_b, basis_c; centering=primitive_centering)
    )

    # Check test conditions
    @test unit_cell isa RhombohedralUnitCell
    @test lattice_constants(unit_cell).α ≈ π / 3

    # Exercise functionality
    iucr_unit_cell = conventional_cell(unit_cell)

    # Check results
    @test iucr_unit_cell ≈ CubicUnitCell(c_a; centering=face_centering)

    # ------ α = π/2

    # Construct rhombohedral unit cell
    unit_cell = RhombohedralUnitCell(c_a, π / 2)

    # Check test conditions
    @test unit_cell isa RhombohedralUnitCell
    @test lattice_constants(unit_cell).α ≈ π / 2

    # Exercise functionality
    iucr_unit_cell = conventional_cell(unit_cell)

    # Check results
    @test iucr_unit_cell ≈ CubicUnitCell(c_a; centering=primitive_centering)

    # ------ α = acos(-1/3)

    # Construct rhombohedral unit cell
    basis_a = 0.5 * (-c_basis_a + c_basis_b + c_basis_c)
    basis_b = 0.5 * (c_basis_a - c_basis_b + c_basis_c)
    basis_c = 0.5 * (c_basis_a + c_basis_b - c_basis_c)
    unit_cell = standardize(
        UnitCell(basis_a, basis_b, basis_c; centering=primitive_centering)
    )

    # Check test conditions
    @test unit_cell isa RhombohedralUnitCell
    @test lattice_constants(unit_cell).α ≈ acos(-1 / 3)

    # Exercise functionality
    iucr_unit_cell = conventional_cell(unit_cell)

    # Check results
    @test iucr_unit_cell ≈ CubicUnitCell(c_a; centering=body_centering)
end

@testset "conventional_cell()::RhombohedralUnitCell: chain of limiting cases" begin
    # --- Exercise functionality and check results

    # primitive unit cell: aP --> mI --> hR
    a = 5
    α = 4π / 7
    expected_unit_cell = RhombohedralUnitCell(a, α)
    basis_a, basis_b, basis_c = basis(expected_unit_cell)
    triclinic_unit_cell = UnitCell(
        basis_a,
        basis_b,
        basis_c;
        identify_lattice_system=false,
        centering=primitive_centering,
    )

    @test triclinic_unit_cell isa TriclinicUnitCell
    @test expected_unit_cell isa RhombohedralUnitCell
    @debug "chain of limiting cases: aP --> mI --> hR"
    @test conventional_cell(triclinic_unit_cell) ≈ expected_unit_cell
end
