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
using Test
using Logging

# XtallographyUtils package
using XtallographyUtils

# --- Tests

@testset "conventional_cell(): limiting cases" begin
    # --- Tests

    # ------ Hexagonal lattices have no limiting cases for primitive centering

    lattice_constants = HexagonalLatticeConstants(1.0, 2.0)

    iucr_unit_cell = conventional_cell(UnitCell(lattice_constants, Primitive()))

    @test iucr_unit_cell.lattice_constants == lattice_constants
    @test iucr_unit_cell.centering == Primitive()

    # ------ Invalid centerings

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
            "(lattice_system=Hexagonal, centering=$(nameof(centering)))"

        @test startswith(error_message, expected_error)
    end
end

@testset "conventional_cell(): chain of limiting cases" begin
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
        Primitive(),
    )
    expected_unit_cell = UnitCell(lattice_constants, Primitive())
    @test triclinic_unit_cell.lattice_constants isa TriclinicLatticeConstants
    @test expected_unit_cell.lattice_constants isa HexagonalLatticeConstants
    @info "aP --> mP --> oS --> hP"
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
        Primitive(),
    )
    expected_unit_cell = UnitCell(lattice_constants, Primitive())
    @test triclinic_unit_cell.lattice_constants isa TriclinicLatticeConstants
    @test expected_unit_cell.lattice_constants isa HexagonalLatticeConstants
    @info "aP --> mP --> oS --> hP"
    @test conventional_cell(triclinic_unit_cell) ≈ expected_unit_cell

    # ------ primitive unit cell: aP --> mI --> oC --> hP

    # TODO
end
