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
using Test
using Logging  # DEBUG

# XtallographyUtils package
using XtallographyUtils

# --- Tests

# ------ standardize()

@testset "standardize()" begin
    # --- Tests

    # ------ Cubic lattices have no lattice constants conventions for PRIMITIVE, BODY, and
    #        FACE centerings

    lattice_constants = CubicLatticeConstants(1.0)

    @test standardize(lattice_constants, XtallographyUtils.PRIMITIVE) ==
        (lattice_constants, XtallographyUtils.PRIMITIVE)
    @test standardize(lattice_constants, XtallographyUtils.BODY) ==
        (lattice_constants, XtallographyUtils.BODY)
    @test standardize(lattice_constants, XtallographyUtils.FACE) ==
        (lattice_constants, XtallographyUtils.FACE)

    # ------ Invalid centering

    local error = nothing
    local error_message = ""
    try
        standardize(lattice_constants, XtallographyUtils.BASE)
    catch error
        bt = catch_backtrace()
        error_message = sprint(showerror, error, bt)
    end

    @test error isa ArgumentError

    expected_error =
        "ArgumentError: " *
        "Invalid Bravais lattice: (lattice_system=Cubic, centering=BASE)"

    @test startswith(error_message, expected_error)
end

# ------ iucr_conventional_cell()

@testset "iucr_conventional_cell(): limiting cases" begin
    # --- Tests

    # ------ Cubic lattices have no limiting cases for PRIMITIVE, BODY, and FACE centerings

    lattice_constants = CubicLatticeConstants(1.0)

    for centering in
        (XtallographyUtils.PRIMITIVE, XtallographyUtils.BODY, XtallographyUtils.FACE)
        unit_cell = UnitCell(lattice_constants, centering)

        iucr_unit_cell = iucr_conventional_cell(unit_cell)

        @test iucr_unit_cell == unit_cell
    end

    # ------ Invalid centering

    local error = nothing
    local error_message = ""
    try
        iucr_conventional_cell(UnitCell(lattice_constants, XtallographyUtils.BASE))
    catch error
        bt = catch_backtrace()
        error_message = sprint(showerror, error, bt)
    end

    @test error isa ArgumentError

    expected_error =
        "ArgumentError: " *
        "Invalid Bravais lattice: (lattice_system=Cubic, centering=BASE)"

    @test startswith(error_message, expected_error)
end

@testset "iucr_conventional_cell(): chain of limiting cases" begin
    # --- Preparations

    a = 5
    lattice_constants = CubicLatticeConstants(a)
    basis_a, basis_b, basis_c = basis(lattice_constants)

    # --- Exercise functionality and check results

    # ------ primitive unit cell: aP --> mP --> oP --> tP --> cP

    triclinic_unit_cell = UnitCell(
        LatticeConstants(basis_a, basis_b, basis_c; identify_lattice_system=false),
        XtallographyUtils.PRIMITIVE,
    )
    expected_unit_cell = UnitCell(lattice_constants, XtallographyUtils.PRIMITIVE)
    @test triclinic_unit_cell.lattice_constants isa TriclinicLatticeConstants
    @test expected_unit_cell.lattice_constants isa CubicLatticeConstants
    @info "aP --> mP --> oP --> tP --> cP"
    @test iucr_conventional_cell(triclinic_unit_cell) ≈ expected_unit_cell

    # ------ body-centered unit cell: aP --> mI --> oI --> tI --> cI

    triclinic_unit_cell = UnitCell(
        LatticeConstants(
            basis_a,
            basis_b,
            0.5 * (basis_a + basis_b + basis_c);
            identify_lattice_system=false,
        ),
        XtallographyUtils.PRIMITIVE,
    )
    expected_unit_cell = UnitCell(lattice_constants, XtallographyUtils.BODY)
    @test triclinic_unit_cell.lattice_constants isa TriclinicLatticeConstants
    @test expected_unit_cell.lattice_constants isa CubicLatticeConstants
    @info "aP --> mI --> oI --> tI --> cI"
    @test iucr_conventional_cell(triclinic_unit_cell) ≈ expected_unit_cell

    # ------ face-centered unit cell: aP --> mI --> oI --> tI --> cF

    triclinic_unit_cell = UnitCell(
        LatticeConstants(
            0.5 * (basis_a + basis_b),
            0.5 * (basis_a - basis_b),
            0.5 * (basis_a + basis_c);
            identify_lattice_system=false,
        ),
        XtallographyUtils.PRIMITIVE,
    )
    expected_unit_cell = UnitCell(lattice_constants, XtallographyUtils.FACE)
    @test triclinic_unit_cell.lattice_constants isa TriclinicLatticeConstants
    @test expected_unit_cell.lattice_constants isa CubicLatticeConstants
    @info "aP --> mI --> oI --> tI --> cF"
    @test iucr_conventional_cell(triclinic_unit_cell) ≈ expected_unit_cell
end

# ------ reduced_cell()

@testset "reduced_cell()" begin
    # --- Preparations

    a = 5
    lattice_constants = CubicLatticeConstants(a)
    basis_a, basis_b, basis_c = basis(lattice_constants)

    # --- Exercise functionality and check results

    # primitive unit cell
    unit_cell = UnitCell(lattice_constants, XtallographyUtils.PRIMITIVE)

    expected_reduced_cell = reduced_cell(
        UnitCell(
            LatticeConstants(basis_a, basis_b, basis_c; identify_lattice_system=false),
            XtallographyUtils.PRIMITIVE,
        ),
    )

    reduced_cell_ = reduced_cell(unit_cell)
    @test reduced_cell_.lattice_constants isa CubicLatticeConstants
    @test volume(reduced_cell_) ≈ volume(unit_cell)
    @test reduced_cell_ ≈ expected_reduced_cell

    # body-centered unit cell
    unit_cell = UnitCell(lattice_constants, XtallographyUtils.BODY)

    expected_reduced_cell = reduced_cell(
        UnitCell(
            LatticeConstants(basis_a, basis_b, 0.5 * (basis_a + basis_b + basis_c)),
            XtallographyUtils.PRIMITIVE,
        ),
    )

    reduced_cell_ = reduced_cell(unit_cell)
    @test reduced_cell_.lattice_constants isa RhombohedralLatticeConstants
    @test volume(reduced_cell_) ≈ 0.5 * volume(unit_cell)
    @test reduced_cell_ ≈ expected_reduced_cell

    # face-centered unit cell
    unit_cell = UnitCell(lattice_constants, XtallographyUtils.FACE)

    expected_reduced_cell = reduced_cell(
        UnitCell(
            LatticeConstants(
                0.5 * (basis_a + basis_b),
                0.5 * (basis_b - basis_c),
                0.5 * (basis_b + basis_c);
                identify_lattice_system=false,
            ),
            XtallographyUtils.PRIMITIVE,
        ),
    )

    reduced_cell_ = reduced_cell(unit_cell)
    @test reduced_cell_.lattice_constants isa TriclinicLatticeConstants
    @test volume(reduced_cell_) ≈ 0.25 * volume(unit_cell)
    @test reduced_cell_ ≈ expected_reduced_cell
end
