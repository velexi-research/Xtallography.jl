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
using Test
using Logging  # DEBUG

# XtallographyUtils package
using XtallographyUtils

# --- Tests

# ------ standardize()

@testset "standardize()" begin
    # --- Tests

    # ------ Tetragonal lattices have no lattice constants conventions for PRIMITIVE and
    #        BODY centerings

    a = 1.0
    c = 10.0
    lattice_constants = TetragonalLatticeConstants(a, c)

    # centering == PRIMITIVE
    standardized_lattice_constants, standardized_centering = standardize(
        lattice_constants, XtallographyUtils.PRIMITIVE
    )

    expected_lattice_constants = lattice_constants
    @test standardized_lattice_constants ≈ expected_lattice_constants

    @test standardized_centering == XtallographyUtils.PRIMITIVE

    # centering == BODY
    standardized_lattice_constants, standardized_centering = standardize(
        lattice_constants, XtallographyUtils.BODY
    )

    expected_lattice_constants = lattice_constants
    @test standardized_lattice_constants ≈ expected_lattice_constants

    @test standardized_centering == XtallographyUtils.BODY

    # ------ Invalid centerings

    for centering in (XtallographyUtils.FACE, XtallographyUtils.BASE)
        local error = nothing
        local error_message = ""
        try
            standardize(lattice_constants, centering)
        catch error
            bt = catch_backtrace()
            error_message = sprint(showerror, error, bt)
        end

        @test error isa ArgumentError

        expected_error =
            "ArgumentError: " *
            "Invalid Bravais lattice: (lattice_system=Tetragonal, centering=$centering)"

        @test startswith(error_message, expected_error)
    end
end

# ------ iucr_conventional_cell()

@testset "iucr_conventional_cell(): limiting cases, centering = PRIMITIVE" begin
    # --- Tests

    # ------ a = b = c

    a = 5.0
    c = a
    lattice_constants = TetragonalLatticeConstants(a, c)
    iucr_unit_cell = iucr_conventional_cell(
        UnitCell(lattice_constants, XtallographyUtils.PRIMITIVE)
    )

    @test iucr_unit_cell.lattice_constants ≈ CubicLatticeConstants(a)
    @test iucr_unit_cell.centering == XtallographyUtils.PRIMITIVE

    # ------ tetragonal unit cell is not equivalent to a cubic unit cell

    a = 5.0
    c = 10.0
    lattice_constants = TetragonalLatticeConstants(a, c)
    iucr_unit_cell = iucr_conventional_cell(
        UnitCell(lattice_constants, XtallographyUtils.PRIMITIVE)
    )

    @test iucr_unit_cell.lattice_constants ≈ lattice_constants
    @test iucr_unit_cell.centering == XtallographyUtils.PRIMITIVE
end

@testset "iucr_conventional_cell(): limiting cases, centering = BODY" begin
    # --- Tests

    # ------ a = b = c

    a = 5.0
    c = a
    lattice_constants = TetragonalLatticeConstants(a, c)
    iucr_unit_cell = iucr_conventional_cell(
        UnitCell(lattice_constants, XtallographyUtils.BODY)
    )

    @test iucr_unit_cell.lattice_constants ≈ CubicLatticeConstants(a)
    @test iucr_unit_cell.centering == XtallographyUtils.BODY

    # ------ c = a √2

    a = 5.0
    c = a * sqrt(2)
    lattice_constants = TetragonalLatticeConstants(a, c)
    iucr_unit_cell = iucr_conventional_cell(
        UnitCell(lattice_constants, XtallographyUtils.BODY)
    )

    @test iucr_unit_cell.lattice_constants ≈ CubicLatticeConstants(c)
    @test iucr_unit_cell.centering == XtallographyUtils.FACE

    # ------ tetragonal unit cell is not equivalent to a cubic unit cell

    a = 5.0
    c = 10.0
    lattice_constants = TetragonalLatticeConstants(a, c)
    iucr_unit_cell = iucr_conventional_cell(
        UnitCell(lattice_constants, XtallographyUtils.BODY)
    )

    @test iucr_unit_cell.lattice_constants ≈ lattice_constants
    @test iucr_unit_cell.centering == XtallographyUtils.BODY
end

@testset "iucr_conventional_cell(): chain of limiting cases" begin
    # --- Preparations

    a = 5
    c = 9
    lattice_constants = TetragonalLatticeConstants(a, c)
    basis_a, basis_b, basis_c = basis(lattice_constants)

    # --- Exercise functionality and check results

    # ------ primitive unit cell: aP --> mP --> oP --> tP

    triclinic_unit_cell = UnitCell(
        LatticeConstants(basis_a, basis_b, basis_c; identify_lattice_system=false),
        XtallographyUtils.PRIMITIVE,
    )
    expected_unit_cell = standardize(
        UnitCell(lattice_constants, XtallographyUtils.PRIMITIVE)
    )
    @test triclinic_unit_cell.lattice_constants isa TriclinicLatticeConstants
    @test expected_unit_cell.lattice_constants isa TetragonalLatticeConstants
    @info "chain of limiting cases: aP --> mP --> oP --> tP"
    @test iucr_conventional_cell(triclinic_unit_cell) ≈ expected_unit_cell

    # ------ primitive unit cell: aP --> mP --> oC --> tI

    # TODO

    # ------ body-centered unit cell: aP --> mI --> oI --> tI

    triclinic_unit_cell = UnitCell(
        LatticeConstants(
            basis_a,
            basis_b,
            0.5 * (basis_a + basis_b + basis_c);
            identify_lattice_system=false,
        ),
        XtallographyUtils.PRIMITIVE,
    )
    expected_unit_cell = standardize(UnitCell(lattice_constants, XtallographyUtils.BODY))
    @test triclinic_unit_cell.lattice_constants isa TriclinicLatticeConstants
    @test expected_unit_cell.lattice_constants isa TetragonalLatticeConstants
    @info "chain of limiting cases: aP --> mI --> oI --> tI"
    @test iucr_conventional_cell(triclinic_unit_cell) ≈ expected_unit_cell

    # ------ body-centered unit cell: aP --> mI --> oF --> tI"

    triclinic_unit_cell = UnitCell(
        LatticeConstants(
            basis_a,
            0.5 * (basis_a + basis_b + basis_c),
            0.5 * (basis_a + basis_b - basis_c);
            identify_lattice_system=false,
        ),
        XtallographyUtils.PRIMITIVE,
    )
    expected_unit_cell = standardize(UnitCell(lattice_constants, XtallographyUtils.BODY))
    @test triclinic_unit_cell.lattice_constants isa TriclinicLatticeConstants
    @test expected_unit_cell.lattice_constants isa TetragonalLatticeConstants
    @info "chain of limiting cases: aP --> mI --> oF --> tI"
    @test iucr_conventional_cell(triclinic_unit_cell) ≈ expected_unit_cell
end

@testset "iucr_conventional_cell(): invalid arguments" begin
    # --- Preparations

    # Construct lattice constants for tetragonal unit cell
    a = 1.0
    c = 3.0
    lattice_constants = TetragonalLatticeConstants(a, c)

    # --- Tests

    for centering in (XtallographyUtils.FACE, XtallographyUtils.BASE)
        local error = nothing
        local error_message = ""
        try
            iucr_conventional_cell(UnitCell(lattice_constants, centering))
        catch error
            bt = catch_backtrace()
            error_message = sprint(showerror, error, bt)
        end

        @test error isa ArgumentError

        expected_error =
            "ArgumentError: " *
            "Invalid Bravais lattice: (lattice_system=Tetragonal, centering=$centering)"

        @test startswith(error_message, expected_error)
    end
end

# ------ reduced_cell()

@testset "reduced_cell()" begin
    # --- Preparations

    a = 5
    c = 7
    lattice_constants = TetragonalLatticeConstants(a, c)
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
    @test reduced_cell_.lattice_constants isa TetragonalLatticeConstants
    @test volume(reduced_cell_) ≈ volume(unit_cell)
    @test reduced_cell_ ≈ expected_reduced_cell

    # body-centered unit cell
    unit_cell = UnitCell(lattice_constants, XtallographyUtils.BODY)

    expected_reduced_cell = reduced_cell(
        UnitCell(
            LatticeConstants(
                basis_a,
                basis_b,
                0.5 * (basis_a + basis_b + basis_c);
                identify_lattice_system=false,
            ),
            XtallographyUtils.PRIMITIVE,
        ),
    )

    reduced_cell_ = reduced_cell(unit_cell)
    @test reduced_cell_.lattice_constants isa TriclinicLatticeConstants
    @test volume(reduced_cell_) ≈ 0.5 * volume(unit_cell)
    @test reduced_cell_ ≈ expected_reduced_cell

    # face-centered unit cell
    unit_cell = UnitCell(lattice_constants, XtallographyUtils.FACE)

    expected_reduced_cell = reduced_cell(
        UnitCell(TetragonalLatticeConstants(a / sqrt(2), c), XtallographyUtils.BODY)
    )

    reduced_cell_ = reduced_cell(unit_cell)
    @test reduced_cell_.lattice_constants isa TriclinicLatticeConstants
    @test volume(reduced_cell_) ≈ 0.25 * volume(unit_cell)
    @test reduced_cell_ ≈ expected_reduced_cell
end
