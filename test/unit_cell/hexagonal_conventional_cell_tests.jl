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
Tests for unit cell standardization methods for hexagonal lattices
"""
# --- Imports

# Standard library
using Logging
using Test

# Xtallography package
using Xtallography

# --- Tests

@testset "conventional_cell(::MonoclinicUnitCell): invalid arguments" begin
    # ------ Invalid centering

    for centering_ in (body_centering, face_centering, base_centering)
        expected_message =
            "Invalid Bravais lattice: " *
            "(lattice_system=Hexagonal, centering=$(nameof(typeof(centering_))))"

        unit_cell = HexagonalUnitCell(1.0, 2.0; centering=centering_)
        @test_throws ArgumentError(expected_message) conventional_cell(unit_cell)
    end
end

@testset "conventional_cell(::HexagonalUnitCell): non-limiting cases" begin
    # --- Tests

    # ------ Hexagonal lattices have no limiting cases for primitive centering

    unit_cell = HexagonalUnitCell(1.0, 2.0)

    iucr_unit_cell = conventional_cell(unit_cell)

    @test iucr_unit_cell == standardize(unit_cell)
end

@testset "conventional_cell()::HexagonalUnitCell: chain of limiting cases" begin
    # --- Exercise functionality and check results

    # ------ primitive unit cell: aP --> mP --> oC --> hP
    #
    # Case #1: a < c

    a = 5
    c = 7
    expected_unit_cell = HexagonalUnitCell(a, c)
    basis_a, basis_b, basis_c = basis(expected_unit_cell)
    triclinic_unit_cell = UnitCell(
        basis_a,
        basis_b,
        basis_c;
        identify_lattice_system=false,
        centering=primitive_centering,
    )
    @test triclinic_unit_cell isa TriclinicUnitCell
    @test expected_unit_cell isa HexagonalUnitCell
    @debug "chain of limiting cases: aP --> mP --> oC --> hP"
    @test conventional_cell(triclinic_unit_cell) ≈ expected_unit_cell

    # ------ primitive unit cell: aP --> mP --> oC --> hP
    #
    # Case #2: a > c

    a = 7
    c = 5
    expected_unit_cell = HexagonalUnitCell(a, c)
    basis_a, basis_b, basis_c = basis(expected_unit_cell)
    triclinic_unit_cell = UnitCell(
        basis_a,
        basis_b,
        basis_c;
        identify_lattice_system=false,
        centering=primitive_centering,
    )
    @test triclinic_unit_cell isa TriclinicUnitCell
    @test expected_unit_cell isa HexagonalUnitCell
    @debug "chain of limiting cases: aP --> mP --> oC --> hP"
    @test conventional_cell(triclinic_unit_cell) ≈ expected_unit_cell

    # ------ primitive unit cell: aP --> mI --> oC --> hP

    a = 7
    c = 5
    expected_unit_cell = HexagonalUnitCell(a, c)
    basis_a, basis_b, basis_c = basis(expected_unit_cell)
    triclinic_unit_cell = UnitCell(
        basis_a,
        basis_b,
        basis_c + basis_a;
        identify_lattice_system=false,
        centering=primitive_centering,
    )
    @test triclinic_unit_cell isa TriclinicUnitCell
    @test expected_unit_cell isa HexagonalUnitCell
    @debug "chain of limiting cases: aP --> mI --> oC --> hP"
    @test conventional_cell(triclinic_unit_cell) ≈ expected_unit_cell
end
