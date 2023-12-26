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

    for centering in (primitive, body_centered, face_centered)
        unit_cell = UnitCell(lattice_constants, centering)

        iucr_unit_cell = conventional_cell(unit_cell)

        @test iucr_unit_cell == unit_cell
    end

    # ------ Invalid centering

    local error = nothing
    local error_message = ""
    try
        conventional_cell(UnitCell(lattice_constants, base_centered))
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
    #
    # Note: presents as the case aP --> mP --> oC --> tP --> cP because the mP --> oP and
    #       mP --> oC limiting cases are equivalent when β = π / 2 and a = c.

    triclinic_unit_cell = UnitCell(
        LatticeConstants(basis_a, basis_b, basis_c; identify_lattice_system=false),
        primitive,
    )
    expected_unit_cell = UnitCell(lattice_constants, primitive)
    @test triclinic_unit_cell.lattice_constants isa TriclinicLatticeConstants
    @test expected_unit_cell.lattice_constants isa CubicLatticeConstants
    @debug "chain of limiting cases: aP --> mP --> oP --> tP --> cP " *
        "(presents as equivalent chain aP --> mP --> oC --> tP --> cP)"
    @test conventional_cell(triclinic_unit_cell) ≈ expected_unit_cell

    # ------ primitive unit cell: aP --> mP --> oC --> tP --> cP

    triclinic_unit_cell = UnitCell(
        LatticeConstants(
            basis_a, basis_a + basis_b, basis_c; identify_lattice_system=false
        ),
        primitive,
    )
    expected_unit_cell = UnitCell(lattice_constants, primitive)
    @test triclinic_unit_cell.lattice_constants isa TriclinicLatticeConstants
    @test expected_unit_cell.lattice_constants isa CubicLatticeConstants
    @debug "chain of limiting cases: aP --> mP --> oC --> tP --> cP"
    @test conventional_cell(triclinic_unit_cell) ≈ expected_unit_cell

    # ------ primitive unit cell: aP --> mI --> oC --> tP --> cP
    #
    # Note: covers the case aP --> mI --> hR --> cP because the mI --> hR and mI --> oC
    #       limiting cases are equivalent when
    #
    #       - a^2 + b^2 = c^2 and b^2 + a * c * cos(β) = a^2
    #         (hR: π/3 < α < π/2)
    #
    #         or
    #       - c^2 + 3 * b^2 = 9 * a^2 and c = -3 * a * cos(β)
    #         (hR: π/2 < α < acos(-1/3)

    triclinic_unit_cell = UnitCell(
        LatticeConstants(
            basis_a, basis_b + basis_c, basis_a + basis_c; identify_lattice_system=false
        ),
        primitive,
    )
    expected_unit_cell = UnitCell(lattice_constants, primitive)
    @test triclinic_unit_cell.lattice_constants isa TriclinicLatticeConstants
    @test expected_unit_cell.lattice_constants isa CubicLatticeConstants
    @debug "chain of limiting cases: aP --> mI --> oC --> tP --> cP " *
        "(equivalent to the case aP --> mI --> hR --> cP)"
    @test conventional_cell(triclinic_unit_cell) ≈ expected_unit_cell

    # ------ body-centered unit cell: aP --> mI --> oI --> tI --> cI

    triclinic_unit_cell = UnitCell(
        LatticeConstants(
            basis_a,
            basis_b,
            0.5 * (basis_a + basis_b + basis_c);
            identify_lattice_system=false,
        ),
        primitive,
    )
    expected_unit_cell = UnitCell(lattice_constants, body_centered)
    @test triclinic_unit_cell.lattice_constants isa TriclinicLatticeConstants
    @test expected_unit_cell.lattice_constants isa CubicLatticeConstants
    @debug "chain of limiting cases: aP --> mI --> oI --> tI --> cI"
    @test conventional_cell(triclinic_unit_cell) ≈ expected_unit_cell

    # ------ body-centered unit cell: aP --> mI --> oF --> tI --> cI
    #
    # Note: covers the case aP --> mI --> hR --> cI because the mI --> hR and mI --> oF
    #       limiting cases are equivalent when
    #
    #           c^2 + 3 * b^2 = 9 * a^2 and c = -3 * a * cos(β)

    triclinic_unit_cell = UnitCell(
        LatticeConstants(
            0.5 * (basis_a + basis_b - basis_c),
            0.5 * (basis_a - basis_b + basis_c),
            0.5 * (-basis_a + basis_b + basis_c);
            identify_lattice_system=false,
        ),
        primitive,
    )
    expected_unit_cell = UnitCell(lattice_constants, body_centered)
    @test triclinic_unit_cell.lattice_constants isa TriclinicLatticeConstants
    @test expected_unit_cell.lattice_constants isa CubicLatticeConstants
    @debug "chain of limiting cases: aP --> mI --> oF --> tI --> cI " *
        "(equivalent to the case aP --> mI --> hR --> cI)"
    @test conventional_cell(triclinic_unit_cell) ≈ expected_unit_cell

    # ------ face-centered unit cell: aP --> mI --> oI --> tI --> cF
    #
    # Note: covers the case aP --> mI --> hR --> cI because the mI --> hR and mI --> oI
    #       limiting cases are equivalent when
    #
    #         a^2 + b^2 = c^2 && a^2 + a * c * cos(β) = b^2

    triclinic_unit_cell = UnitCell(
        LatticeConstants(
            0.5 * (basis_a + basis_b),
            0.5 * (basis_b + basis_c),
            0.5 * (basis_c + basis_a);
            identify_lattice_system=false,
        ),
        primitive,
    )
    expected_unit_cell = UnitCell(lattice_constants, face_centered)
    @test triclinic_unit_cell.lattice_constants isa TriclinicLatticeConstants
    @test expected_unit_cell.lattice_constants isa CubicLatticeConstants
    @debug "chain of limiting cases: aP --> mI --> oI --> tI --> cF " *
        "(equivalent to the case aP --> mI --> hR --> cF)"
    @test conventional_cell(triclinic_unit_cell) ≈ expected_unit_cell

    # ------ face-centered unit cell: aP --> mI --> oF --> tI --> cF
    #
    # Note: presents as the case aP --> mI --> oI --> tI --> cF because the mI --> oF and
    #       mP --> oI limiting cases are equivalent when β = π / 2 and a = c.

    triclinic_unit_cell = UnitCell(
        LatticeConstants(
            basis_a,
            0.5 * (basis_a + basis_b),
            0.5 * (basis_b + basis_c);
            identify_lattice_system=false,
        ),
        primitive,
    )
    expected_unit_cell = UnitCell(lattice_constants, face_centered)
    @test triclinic_unit_cell.lattice_constants isa TriclinicLatticeConstants
    @test expected_unit_cell.lattice_constants isa CubicLatticeConstants
    @debug "chain of limiting cases: aP --> mI --> oF --> tI --> cF " *
        "(presents as equivalent chain aP --> mI --> oI --> tI --> cF)"
    @test conventional_cell(triclinic_unit_cell) ≈ expected_unit_cell
end
