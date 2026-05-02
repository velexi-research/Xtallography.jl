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

@testset "conventional_cell(::TetragonalUnitCell): invalid arguments" begin
    # --- Preparations

    # Construct lattice constants for tetragonal unit cell
    a = 1.0
    c = 3.0

    # --- Tests

    # ------ Invalid centering

    for centering_ in (face_centering, base_centering)
        expected_message =
            "Invalid Bravais lattice: " *
            "(lattice_system=Tetragonal, centering=$(nameof(typeof(centering_))))"

        @test_throws ArgumentError(expected_message) conventional_cell(
            TetragonalUnitCell(a, c; centering=centering_)
        )
    end
end

@testset "conventional_cell(::TetragonalUnitCell): non-limiting cases" begin
    # tetragonal unit cells that are not equivalent to any higher symmetry unit cell

    # --- centering = primitive_centering

    a = 5.0
    c = 10.0
    unit_cell = TetragonalUnitCell(a, c)

    # Check test conditions
    @test unit_cell isa TetragonalUnitCell
    @test centering(unit_cell) === primitive_centering

    # Exercise functionality
    iucr_unit_cell = conventional_cell(unit_cell)

    # Check results
    expected_unit_cell = standardize(unit_cell)
    @test iucr_unit_cell == expected_unit_cell

    # --- centering = body_centering

    a = 5.0
    c = 10.0
    unit_cell = TetragonalUnitCell(a, c; centering=body_centering)

    # Check test conditions
    @test unit_cell isa TetragonalUnitCell
    @test centering(unit_cell) === body_centering

    # Exercise functionality
    iucr_unit_cell = conventional_cell(unit_cell)

    # Check results
    expected_unit_cell = standardize(unit_cell)
    @test iucr_unit_cell == expected_unit_cell
end

@testset "conventional_cell(::TetragonalUnitCell): limiting cases - tP --> cP" begin
    # --- Tests

    # ------ a = b = c

    a = 5.0
    c = a
    unit_cell = TetragonalUnitCell(a, c; centering=primitive_centering)

    iucr_unit_cell = conventional_cell(unit_cell)

    @test iucr_unit_cell ≈ CubicUnitCell(a; centering=primitive_centering)
end

@testset "conventional_cell(::TetragonalUnitCell): limiting cases - tI --> cI" begin
    # --- Tests

    # ------ a = b = c

    a = 5.0
    c = a
    unit_cell = TetragonalUnitCell(a, c; centering=body_centering)

    iucr_unit_cell = conventional_cell(unit_cell)

    @test iucr_unit_cell ≈ CubicUnitCell(a; centering=body_centering)

    # ------ c = a √2

    a = 5.0
    c = a * sqrt(2)
    unit_cell = TetragonalUnitCell(a, c; centering=body_centering)

    iucr_unit_cell = conventional_cell(unit_cell)

    @test iucr_unit_cell ≈ CubicUnitCell(c; centering=face_centering)
end

@testset "conventional_cell()::TetragonalUnitCell: chain of limiting cases" begin
    # --- Preparations

    a = 5
    c = 9
    primitive_tetragonal_unit_cell = TetragonalUnitCell(a, c)
    basis_a, basis_b, basis_c = basis(primitive_tetragonal_unit_cell)

    # --- Exercise functionality and check results

    # ------ primitive unit cell: aP --> mP --> oP --> tP

    triclinic_unit_cell = UnitCell(
        basis_a,
        basis_b,
        basis_c;
        identify_lattice_system=false,
        centering=primitive_centering,
    )
    expected_unit_cell = standardize(primitive_tetragonal_unit_cell)
    @test triclinic_unit_cell isa TriclinicUnitCell
    @test expected_unit_cell isa TetragonalUnitCell
    @test centering(expected_unit_cell) === primitive_centering
    @debug "chain of limiting cases: aP --> mP --> oP --> tP"
    @test conventional_cell(triclinic_unit_cell) ≈ expected_unit_cell

    # ------ primitive unit cell: aP --> mP --> oC --> tP

    triclinic_unit_cell = UnitCell(
        basis_a,
        basis_a + basis_b,
        basis_c;
        identify_lattice_system=false,
        centering=primitive_centering,
    )
    expected_unit_cell = standardize(primitive_tetragonal_unit_cell)
    @test triclinic_unit_cell isa TriclinicUnitCell
    @test expected_unit_cell isa TetragonalUnitCell
    @test centering(expected_unit_cell) === primitive_centering
    @debug "chain of limiting cases: aP --> mP --> oC --> tP"
    @test conventional_cell(triclinic_unit_cell) ≈ expected_unit_cell

    # ------ primitive unit cell: aP --> mI --> oC --> tP

    triclinic_unit_cell = UnitCell(
        basis_a,
        basis_b + basis_a,
        basis_c + basis_b;
        identify_lattice_system=false,
        centering=primitive_centering,
    )
    expected_unit_cell = standardize(primitive_tetragonal_unit_cell)
    @test triclinic_unit_cell isa TriclinicUnitCell
    @test expected_unit_cell isa TetragonalUnitCell
    @test centering(expected_unit_cell) === primitive_centering
    @debug "chain of limiting cases: aP --> mI --> oC --> tP"
    @test conventional_cell(triclinic_unit_cell) ≈ expected_unit_cell

    # ------ body-centered unit cell: aP --> mI --> oI --> tI

    triclinic_unit_cell = UnitCell(
        basis_a,
        basis_b,
        0.5 * (basis_a + basis_b + basis_c);
        identify_lattice_system=false,
        centering=primitive_centering,
    )
    expected_unit_cell = standardize(TetragonalUnitCell(a, c; centering=body_centering))
    @test triclinic_unit_cell isa TriclinicUnitCell
    @test expected_unit_cell isa TetragonalUnitCell
    @test centering(expected_unit_cell) === body_centering
    @debug "chain of limiting cases: aP --> mI --> oI --> tI"
    @test conventional_cell(triclinic_unit_cell) ≈ expected_unit_cell

    # ------ body-centered unit cell: aP --> mI --> oF --> tI

    triclinic_unit_cell = UnitCell(
        basis_a,
        0.5 * (basis_a + basis_b + basis_c),
        0.5 * (basis_a + basis_b - basis_c);
        identify_lattice_system=false,
        centering=primitive_centering,
    )
    expected_unit_cell = standardize(TetragonalUnitCell(a, c; centering=body_centering))
    @test triclinic_unit_cell isa TriclinicUnitCell
    @test expected_unit_cell isa TetragonalUnitCell
    @test centering(expected_unit_cell) === body_centering
    @debug "chain of limiting cases: aP --> mI --> oF --> tI"
    @test conventional_cell(triclinic_unit_cell) ≈ expected_unit_cell
end
