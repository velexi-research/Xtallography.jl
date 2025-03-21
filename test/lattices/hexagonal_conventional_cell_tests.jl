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

@testset "conventional_cell():hexagonal: limiting cases" begin
    # --- Tests

    # ------ Hexagonal lattices have no limiting cases for primitive centering

    lattice_constants = HexagonalLatticeConstants(1.0, 2.0)

    iucr_unit_cell = conventional_cell(UnitCell(lattice_constants, primitive))

    @test iucr_unit_cell.lattice_constants == lattice_constants
    @test iucr_unit_cell.centering === primitive

    # ------ Invalid centering

    for centering in (body_centered, face_centered, base_centered)
        expected_message =
            "Invalid Bravais lattice: " *
            "(lattice_system=Hexagonal, centering=$(nameof(typeof(centering))))"

        @test_throws ArgumentError(expected_message) conventional_cell(
            UnitCell(lattice_constants, centering)
        )
    end
end

@testset "conventional_cell():hexagonal: chain of limiting cases" begin
    # --- Exercise functionality and check results

    # ------ primitive unit cell: aP --> mP --> oC --> hP
    #
    # Case #1: a < c

    a = 5
    c = 7
    lattice_constants = HexagonalLatticeConstants(a, c)
    basis_a, basis_b, basis_c = basis(lattice_constants)
    triclinic_unit_cell = UnitCell(
        LatticeConstants(basis_a, basis_b, basis_c; identify_lattice_system=false),
        primitive,
    )
    expected_unit_cell = UnitCell(lattice_constants, primitive)
    @test triclinic_unit_cell.lattice_constants isa TriclinicLatticeConstants
    @test expected_unit_cell.lattice_constants isa HexagonalLatticeConstants
    @debug "chain of limiting cases: aP --> mP --> oC --> hP"
    @test conventional_cell(triclinic_unit_cell) ≈ expected_unit_cell

    # ------ primitive unit cell: aP --> mP --> oC --> hP
    #
    # Case #2: a > c

    a = 7
    c = 5
    lattice_constants = HexagonalLatticeConstants(a, c)
    basis_a, basis_b, basis_c = basis(lattice_constants)
    triclinic_unit_cell = UnitCell(
        LatticeConstants(basis_a, basis_b, basis_c; identify_lattice_system=false),
        primitive,
    )
    expected_unit_cell = UnitCell(lattice_constants, primitive)
    @test triclinic_unit_cell.lattice_constants isa TriclinicLatticeConstants
    @test expected_unit_cell.lattice_constants isa HexagonalLatticeConstants
    @debug "chain of limiting cases: aP --> mP --> oC --> hP"
    @test conventional_cell(triclinic_unit_cell) ≈ expected_unit_cell

    # ------ primitive unit cell: aP --> mI --> oC --> hP

    a = 7
    c = 5
    lattice_constants = HexagonalLatticeConstants(a, c)
    basis_a, basis_b, basis_c = basis(lattice_constants)
    triclinic_unit_cell = UnitCell(
        LatticeConstants(
            basis_a, basis_b, basis_c + basis_a; identify_lattice_system=false
        ),
        primitive,
    )
    expected_unit_cell = UnitCell(lattice_constants, primitive)
    @test triclinic_unit_cell.lattice_constants isa TriclinicLatticeConstants
    @test expected_unit_cell.lattice_constants isa HexagonalLatticeConstants
    @debug "chain of limiting cases: aP --> mI --> oC --> hP"
    @test conventional_cell(triclinic_unit_cell) ≈ expected_unit_cell
end
