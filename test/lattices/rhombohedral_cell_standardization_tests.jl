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
using Test

# XtallographyUtils package
using XtallographyUtils

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

@testset "conventional_cell(): limiting cases" begin
    # --- Preparations

    # Construct basis for cubic unit cell
    c_a = 1.0

    expected_cubic_lattice_constants = CubicLatticeConstants(c_a)

    c_basis_a, c_basis_b, c_basis_c = basis(expected_cubic_lattice_constants)

    # --- Tests

    # ------ α = π/3

    # Construct lattice constants for rhombohedral unit cell
    basis_a = 0.5 * (c_basis_b + c_basis_c)
    basis_b = 0.5 * (c_basis_a + c_basis_c)
    basis_c = 0.5 * (c_basis_a + c_basis_b)
    lattice_constants, _ = standardize(
        LatticeConstants(basis_a, basis_b, basis_c), Primitive()
    )

    # Check test conditions
    @test lattice_constants isa RhombohedralLatticeConstants
    @test lattice_constants.α ≈ π / 3

    # Exercise functionality
    iucr_unit_cell = conventional_cell(UnitCell(lattice_constants, Primitive()))

    # Check results
    @test iucr_unit_cell.lattice_constants ≈ CubicLatticeConstants(c_a)
    @test iucr_unit_cell.centering == FaceCentered()

    # ------ α = π/2

    # Construct lattice constants for rhombohedral unit cell
    lattice_constants = RhombohedralLatticeConstants(c_a, π / 2)

    # Check test conditions
    @test lattice_constants isa RhombohedralLatticeConstants
    @test lattice_constants.α ≈ π / 2

    # Exercise functionality
    iucr_unit_cell = conventional_cell(UnitCell(lattice_constants, Primitive()))

    # Check results
    @test iucr_unit_cell.lattice_constants ≈ expected_cubic_lattice_constants
    @test iucr_unit_cell.centering == Primitive()

    # ------ α = acos(-1/3)

    # Construct lattice constants for rhombohedral unit cell
    basis_a = 0.5 * (-c_basis_a + c_basis_b + c_basis_c)
    basis_b = 0.5 * (c_basis_a - c_basis_b + c_basis_c)
    basis_c = 0.5 * (c_basis_a + c_basis_b - c_basis_c)
    lattice_constants, _ = standardize(
        LatticeConstants(basis_a, basis_b, basis_c), Primitive()
    )

    # Check test conditions
    @test lattice_constants isa RhombohedralLatticeConstants
    @test lattice_constants.α ≈ acos(-1 / 3)

    # Exercise functionality
    iucr_unit_cell = conventional_cell(UnitCell(lattice_constants, Primitive()))

    # Check results
    @test iucr_unit_cell.lattice_constants ≈ expected_cubic_lattice_constants
    @test iucr_unit_cell.centering == BodyCentered()

    # ------ rhombohedral unit cell is not equivalent to an cubic unit cell

    # Construct lattice constants for rhombohedral unit cell
    a = 1.0
    α = 2π / 5
    lattice_constants = RhombohedralLatticeConstants(a, α)

    # Exercise functionality
    iucr_unit_cell = conventional_cell(UnitCell(lattice_constants, Primitive()))

    # Check results
    @test iucr_unit_cell.lattice_constants ≈ lattice_constants
    @test iucr_unit_cell.centering == Primitive()
end

@testset "conventional_cell(): invalid arguments" begin
    # --- Preparations

    # Construct lattice constants for rhombohedral unit cell
    a = 1.0
    α = 2π / 5
    lattice_constants = RhombohedralLatticeConstants(a, α)

    # --- Tests

    for centering in (BodyCentered, FaceCentered, BaseCentered)
        local error = nothing
        local error_message = ""
        try
            conventional_cell(UnitCell(lattice_constants, centering()))
        catch error
            bt = catch_backtrace()
            error_message = sprint(showerror, error, bt)
        end

        @test error isa ArgumentError

        expected_error =
            "ArgumentError: " *
            "Invalid Bravais lattice: " *
            "(lattice_system=Rhombohedral, centering=$(nameof(centering)))"

        @test startswith(error_message, expected_error)
    end
end

@testset "conventional_cell(): chain of limiting cases" begin
    # --- Exercise functionality and check results

    # primitive unit cell
    a = 5
    α = 4π / 7
    lattice_constants = RhombohedralLatticeConstants(a, α)
    basis_a, basis_b, basis_c = basis(lattice_constants)
    triclinic_unit_cell = UnitCell(
        LatticeConstants(basis_a, basis_b, basis_c; identify_lattice_system=false),
        Primitive(),
    )
    expected_unit_cell = UnitCell(lattice_constants, Primitive())
    @test triclinic_unit_cell.lattice_constants isa TriclinicLatticeConstants
    @test expected_unit_cell.lattice_constants isa RhombohedralLatticeConstants
    @test conventional_cell(triclinic_unit_cell) ≈ expected_unit_cell
end
