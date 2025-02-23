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
Tests for unit cell standardization methods for tetragonal lattices
"""
# --- Imports

# Standard library
using Logging
using Test

# Xtallography package
using Xtallography

# --- Tests

@testset "conventional_cell():tetragonal: limiting cases, centering = primitive" begin
    # --- Tests

    # ------ a = b = c

    a = 5.0
    c = a
    lattice_constants = TetragonalLatticeConstants(a, c)
    iucr_unit_cell = conventional_cell(UnitCell(lattice_constants, primitive))

    @test iucr_unit_cell.lattice_constants ≈ CubicLatticeConstants(a)
    @test iucr_unit_cell.centering === primitive

    # ------ tetragonal unit cell is not equivalent to a cubic unit cell

    a = 5.0
    c = 10.0
    lattice_constants = TetragonalLatticeConstants(a, c)
    iucr_unit_cell = conventional_cell(UnitCell(lattice_constants, primitive))

    @test iucr_unit_cell.lattice_constants ≈ lattice_constants
    @test iucr_unit_cell.centering === primitive
end

@testset "conventional_cell():tetragonal: limiting cases, centering = body" begin
    # --- Tests

    # ------ a = b = c

    a = 5.0
    c = a
    lattice_constants = TetragonalLatticeConstants(a, c)
    iucr_unit_cell = conventional_cell(UnitCell(lattice_constants, body_centered))

    @test iucr_unit_cell.lattice_constants ≈ CubicLatticeConstants(a)
    @test iucr_unit_cell.centering === body_centered

    # ------ c = a √2

    a = 5.0
    c = a * sqrt(2)
    lattice_constants = TetragonalLatticeConstants(a, c)
    iucr_unit_cell = conventional_cell(UnitCell(lattice_constants, body_centered))

    @test iucr_unit_cell.lattice_constants ≈ CubicLatticeConstants(c)
    @test iucr_unit_cell.centering === face_centered

    # ------ tetragonal unit cell is not equivalent to a cubic unit cell

    a = 5.0
    c = 10.0
    lattice_constants = TetragonalLatticeConstants(a, c)
    iucr_unit_cell = conventional_cell(UnitCell(lattice_constants, body_centered))

    @test iucr_unit_cell.lattice_constants ≈ lattice_constants
    @test iucr_unit_cell.centering === body_centered
end

@testset "conventional_cell():tetragonal: chain of limiting cases" begin
    # --- Preparations

    a = 5
    c = 9
    lattice_constants = TetragonalLatticeConstants(a, c)
    basis_a, basis_b, basis_c = basis(lattice_constants)

    # --- Exercise functionality and check results

    # ------ primitive unit cell: aP --> mP --> oP --> tP

    triclinic_unit_cell = UnitCell(
        LatticeConstants(basis_a, basis_b, basis_c; identify_lattice_system=false),
        primitive,
    )
    expected_unit_cell = standardize(UnitCell(lattice_constants, primitive))
    @test triclinic_unit_cell.lattice_constants isa TriclinicLatticeConstants
    @test expected_unit_cell.lattice_constants isa TetragonalLatticeConstants
    @debug "chain of limiting cases: aP --> mP --> oP --> tP"
    @test conventional_cell(triclinic_unit_cell) ≈ expected_unit_cell

    # ------ primitive unit cell: aP --> mP --> oC --> tP

    triclinic_unit_cell = UnitCell(
        LatticeConstants(
            basis_a, basis_a + basis_b, basis_c; identify_lattice_system=false
        ),
        primitive,
    )
    expected_unit_cell = standardize(UnitCell(lattice_constants, primitive))
    @test triclinic_unit_cell.lattice_constants isa TriclinicLatticeConstants
    @test expected_unit_cell.lattice_constants isa TetragonalLatticeConstants
    @debug "chain of limiting cases: aP --> mP --> oC --> tP"
    @test conventional_cell(triclinic_unit_cell) ≈ expected_unit_cell

    # ------ primitive unit cell: aP --> mI --> oC --> tP

    triclinic_unit_cell = UnitCell(
        LatticeConstants(
            basis_a, basis_b + basis_a, basis_c + basis_b; identify_lattice_system=false
        ),
        primitive,
    )
    expected_unit_cell = standardize(UnitCell(lattice_constants, primitive))
    @test triclinic_unit_cell.lattice_constants isa TriclinicLatticeConstants
    @test expected_unit_cell.lattice_constants isa TetragonalLatticeConstants
    @debug "chain of limiting cases: aP --> mI --> oC --> tP"
    @test conventional_cell(triclinic_unit_cell) ≈ expected_unit_cell

    # ------ body-centered unit cell: aP --> mI --> oI --> tI

    triclinic_unit_cell = UnitCell(
        LatticeConstants(
            basis_a,
            basis_b,
            0.5 * (basis_a + basis_b + basis_c);
            identify_lattice_system=false,
        ),
        primitive,
    )
    expected_unit_cell = standardize(UnitCell(lattice_constants, body_centered))
    @test triclinic_unit_cell.lattice_constants isa TriclinicLatticeConstants
    @test expected_unit_cell.lattice_constants isa TetragonalLatticeConstants
    @debug "chain of limiting cases: aP --> mI --> oI --> tI"
    @test conventional_cell(triclinic_unit_cell) ≈ expected_unit_cell

    # ------ body-centered unit cell: aP --> mI --> oF --> tI

    triclinic_unit_cell = UnitCell(
        LatticeConstants(
            basis_a,
            0.5 * (basis_a + basis_b + basis_c),
            0.5 * (basis_a + basis_b - basis_c);
            identify_lattice_system=false,
        ),
        primitive,
    )
    expected_unit_cell = standardize(UnitCell(lattice_constants, body_centered))
    @test triclinic_unit_cell.lattice_constants isa TriclinicLatticeConstants
    @test expected_unit_cell.lattice_constants isa TetragonalLatticeConstants
    @debug "chain of limiting cases: aP --> mI --> oF --> tI"
    @test conventional_cell(triclinic_unit_cell) ≈ expected_unit_cell
end

@testset "conventional_cell():tetragonal: invalid arguments" begin
    # --- Preparations

    # Construct lattice constants for tetragonal unit cell
    a = 1.0
    c = 3.0
    lattice_constants = TetragonalLatticeConstants(a, c)

    # --- Tests

    for centering in (face_centered, base_centered)
        local error = nothing
        local error_message = ""
        try
            conventional_cell(UnitCell(lattice_constants, centering))
        catch error
            bt = catch_backtrace()
            error_message = sprint(showerror, error, bt)
        end

        @test error isa ArgumentError

        expected_error =
            "ArgumentError: " *
            "Invalid Bravais lattice: " *
            "(lattice_system=Tetragonal, centering=$(nameof(typeof(centering))))"

        @test startswith(error_message, expected_error)
    end
end
