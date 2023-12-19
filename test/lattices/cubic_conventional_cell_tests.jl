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
Tests for unit cell standardization methods for cubic lattices
"""
# --- Imports

# Standard library
using Logging
using Test

# XtallographyUtils package
using XtallographyUtils

# --- Tests

@testset "conventional_cell(): limiting cases" begin
    # --- Tests

    # ------ Cubic lattices have no limiting cases for primitive, body, and face centerings

    lattice_constants = CubicLatticeConstants(1.0)

    for centering in (Primitive(), BodyCentered(), FaceCentered())
        unit_cell = UnitCell(lattice_constants, centering)

        iucr_unit_cell = conventional_cell(unit_cell)

        @test iucr_unit_cell == unit_cell
    end

    # ------ Invalid centering

    local error = nothing
    local error_message = ""
    try
        conventional_cell(UnitCell(lattice_constants, BaseCentered()))
    catch error
        bt = catch_backtrace()
        error_message = sprint(showerror, error, bt)
    end

    @test error isa ArgumentError

    expected_error =
        "ArgumentError: " *
        "Invalid Bravais lattice: (lattice_system=Cubic, centering=BaseCentered)"

    @test startswith(error_message, expected_error)
end

@testset "conventional_cell(): chain of limiting cases" begin
    # --- Preparations

    a = 5
    lattice_constants = CubicLatticeConstants(a)
    basis_a, basis_b, basis_c = basis(lattice_constants)

    # --- Exercise functionality and check results

    # ------ primitive unit cell: aP --> mP --> oP --> tP --> cP

    triclinic_unit_cell = UnitCell(
        LatticeConstants(basis_a, basis_b, basis_c; identify_lattice_system=false),
        Primitive(),
    )
    expected_unit_cell = UnitCell(lattice_constants, Primitive())
    @test triclinic_unit_cell.lattice_constants isa TriclinicLatticeConstants
    @test expected_unit_cell.lattice_constants isa CubicLatticeConstants
    @debug "chain of limiting cases: aP --> mP --> oP --> tP --> cP"
    @test conventional_cell(triclinic_unit_cell) ≈ expected_unit_cell

    # ------ body-centered unit cell: aP --> mI --> oI --> tI --> cI

    triclinic_unit_cell = UnitCell(
        LatticeConstants(
            basis_a,
            basis_b,
            0.5 * (basis_a + basis_b + basis_c);
            identify_lattice_system=false,
        ),
        Primitive(),
    )
    expected_unit_cell = UnitCell(lattice_constants, BodyCentered())
    @test triclinic_unit_cell.lattice_constants isa TriclinicLatticeConstants
    @test expected_unit_cell.lattice_constants isa CubicLatticeConstants
    @debug "chain of limiting cases: aP --> mI --> oI --> tI --> cI"
    @test conventional_cell(triclinic_unit_cell) ≈ expected_unit_cell

    # ------ face-centered unit cell: aP --> mI --> oI --> tI --> cF

    triclinic_unit_cell = UnitCell(
        LatticeConstants(
            0.5 * (basis_a + basis_b),
            0.5 * (basis_a - basis_b),
            0.5 * (basis_a + basis_c);
            identify_lattice_system=false,
        ),
        Primitive(),
    )
    expected_unit_cell = UnitCell(lattice_constants, FaceCentered())
    @test triclinic_unit_cell.lattice_constants isa TriclinicLatticeConstants
    @test expected_unit_cell.lattice_constants isa CubicLatticeConstants
    @debug "chain of limiting cases: aP --> mI --> oI --> tI --> cF"
    @test conventional_cell(triclinic_unit_cell) ≈ expected_unit_cell
end
