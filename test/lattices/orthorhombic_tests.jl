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
Tests for methods in lattice/orthorhombic.jl (except for cell standardization methods)
"""
# --- Imports

# Standard library
using Test

# XtallographyUtils package
using XtallographyUtils

# --- Tests

# ------ Types

@testset "OrthorhombicLatticeConstants constructor: valid arguments" begin
    # --- Preparations

    a = 1
    b = 2
    c = 3

    # --- Exercise functionality

    lattice_constants = OrthorhombicLatticeConstants(a, b, c)

    # --- Check results

    @test lattice_constants.a == a
    @test lattice_constants.b == b
    @test lattice_constants.c == c
end

@testset "OrthorhombicLatticeConstants constructor: invalid arguments" begin
    # --- Preparations

    # Valid arguments
    a = 1
    b = 2
    c = 3

    # --- Tests

    # ------ a

    # a = 0
    local error = nothing
    local error_message = ""
    try
        lattice_constants = OrthorhombicLatticeConstants(0, b, c)
    catch error
        bt = catch_backtrace()
        error_message = sprint(showerror, error, bt)
    end

    @test error isa ArgumentError

    expected_error = "ArgumentError: `a` must be positive"
    @test startswith(error_message, expected_error)

    # a < 0
    local error = nothing
    local error_message = ""
    try
        lattice_constants = OrthorhombicLatticeConstants(-1.0, b, c)
    catch error
        bt = catch_backtrace()
        error_message = sprint(showerror, error, bt)
    end

    @test error isa ArgumentError

    expected_error = "ArgumentError: `a` must be positive"
    @test startswith(error_message, expected_error)

    # ------ b

    # b = 0
    local error = nothing
    local error_message = ""
    try
        lattice_constants = OrthorhombicLatticeConstants(a, 0, c)
    catch error
        bt = catch_backtrace()
        error_message = sprint(showerror, error, bt)
    end

    @test error isa ArgumentError

    expected_error = "ArgumentError: `b` must be positive"
    @test startswith(error_message, expected_error)

    # b < 0
    local error = nothing
    local error_message = ""
    try
        lattice_constants = OrthorhombicLatticeConstants(a, -1.0, c)
    catch error
        bt = catch_backtrace()
        error_message = sprint(showerror, error, bt)
    end

    @test error isa ArgumentError

    expected_error = "ArgumentError: `b` must be positive"
    @test startswith(error_message, expected_error)

    # ------ c

    # c = 0
    local error = nothing
    local error_message = ""
    try
        lattice_constants = OrthorhombicLatticeConstants(a, b, 0)
    catch error
        bt = catch_backtrace()
        error_message = sprint(showerror, error, bt)
    end

    @test error isa ArgumentError

    expected_error = "ArgumentError: `c` must be positive"
    @test startswith(error_message, expected_error)

    # c < 0
    local error = nothing
    local error_message = ""
    try
        lattice_constants = OrthorhombicLatticeConstants(a, b, -1.0)
    catch error
        bt = catch_backtrace()
        error_message = sprint(showerror, error, bt)
    end

    @test error isa ArgumentError

    expected_error = "ArgumentError: `c` must be positive"
    @test startswith(error_message, expected_error)
end

# ------ LatticeConstants functions

@testset "isapprox(::LatticeConstants)" begin
    # --- Preparations

    x = OrthorhombicLatticeConstants(1.0, 2.0, 3.0)
    y = OrthorhombicLatticeConstants(1.5, 2.5, 3.5)

    # --- Exercise functionality and check results

    # x ≈ (x + delta)
    @test x ≈ OrthorhombicLatticeConstants(1.0 + 1e-9, 2.0, 3.0)
    @test x ≈ OrthorhombicLatticeConstants(1.0, 2.0 + 1e-9, 3.0)
    @test x ≈ OrthorhombicLatticeConstants(1.0, 2.0, 3.0 - 1e-9)

    # x !≈ y
    @test !(x ≈ y)

    # x ≈ y: atol = 1
    @test isapprox(x, y; atol=1)

    # x ≈ y: rtol = 1
    @test isapprox(x, y; rtol=1)

    # x ≈ y: atol = 0.01, rtol = 1
    @test isapprox(x, y; atol=0.01, rtol=1)

    # x ≈ y: atol = 1, rtol = 0.01
    @test isapprox(x, y; atol=1, rtol=0.01)

    # x !≈ y: atol = 0.01, rtol = 0.01
    @test !isapprox(x, y; atol=0.01, rtol=0.01)
end

@testset "lattice_system()" begin
    # --- Tests

    lattice_constants = OrthorhombicLatticeConstants(1, 2, 3)
    @test lattice_system(lattice_constants) == Orthorhombic
end

@testset "standardize(): primitive, body-centered, face-centered" begin
    # --- Tests

    # ------ lattice constants already in standard form

    a = 1.0
    b = 5.0
    c = 10.0
    lattice_constants = OrthorhombicLatticeConstants(a, b, c)

    # centering == PRIMITIVE
    standardized_lattice_constants, standardized_centering = standardize(
        lattice_constants, XtallographyUtils.PRIMITIVE
    )

    expected_lattice_constants = OrthorhombicLatticeConstants(a, b, c)
    @test standardized_lattice_constants ≈ expected_lattice_constants

    @test standardized_centering == XtallographyUtils.PRIMITIVE

    # centering == BODY
    standardized_lattice_constants, standardized_centering = standardize(
        lattice_constants, XtallographyUtils.BODY
    )

    expected_lattice_constants = OrthorhombicLatticeConstants(a, b, c)
    @test standardized_lattice_constants ≈ expected_lattice_constants

    @test standardized_centering == XtallographyUtils.BODY

    # centering == FACE
    standardized_lattice_constants, standardized_centering = standardize(
        lattice_constants, XtallographyUtils.FACE
    )

    expected_lattice_constants = OrthorhombicLatticeConstants(a, b, c)
    @test standardized_lattice_constants ≈ expected_lattice_constants

    @test standardized_centering == XtallographyUtils.FACE

    # ------ lattice constants not sorted

    a = 1.0
    b = 10.0
    c = 5.0
    lattice_constants = OrthorhombicLatticeConstants(a, b, c)

    # centering == PRIMITIVE
    standardized_lattice_constants, standardized_centering = standardize(
        lattice_constants, XtallographyUtils.PRIMITIVE
    )

    expected_lattice_constants = OrthorhombicLatticeConstants(a, c, b)
    @test standardized_lattice_constants ≈ expected_lattice_constants

    @test standardized_centering == XtallographyUtils.PRIMITIVE

    # centering == BODY
    standardized_lattice_constants, standardized_centering = standardize(
        lattice_constants, XtallographyUtils.BODY
    )

    expected_lattice_constants = OrthorhombicLatticeConstants(a, c, b)
    @test standardized_lattice_constants ≈ expected_lattice_constants

    @test standardized_centering == XtallographyUtils.BODY

    # centering == FACE
    standardized_lattice_constants, standardized_centering = standardize(
        lattice_constants, XtallographyUtils.FACE
    )

    expected_lattice_constants = OrthorhombicLatticeConstants(a, c, b)
    @test standardized_lattice_constants ≈ expected_lattice_constants

    @test standardized_centering == XtallographyUtils.FACE
end

@testset "standardize(): base-centered" begin
    # --- Tests

    # ------ lattice constants already in standard form

    a = 1.0
    b = 5.0
    c = 10.0
    lattice_constants = OrthorhombicLatticeConstants(a, b, c)

    standardized_lattice_constants, standardized_centering = standardize(
        lattice_constants, XtallographyUtils.BASE
    )

    expected_lattice_constants = OrthorhombicLatticeConstants(a, b, c)
    @test standardized_lattice_constants ≈ expected_lattice_constants

    @test standardized_centering == XtallographyUtils.BASE

    # ------ lattice constants not sorted

    a = 5.0
    b = 1.0
    c = 10.0
    lattice_constants = OrthorhombicLatticeConstants(a, b, c)

    standardized_lattice_constants, standardized_centering = standardize(
        lattice_constants, XtallographyUtils.BASE
    )

    expected_lattice_constants = OrthorhombicLatticeConstants(b, a, c)
    @test standardized_lattice_constants ≈ expected_lattice_constants

    @test standardized_centering == XtallographyUtils.BASE
end

# ------ Unit cell computations

@testset "basis()" begin
    # --- Preparations

    a = 2
    b = 3
    c = 5
    lattice_constants = OrthorhombicLatticeConstants(a, b, c)

    # --- Exercise functionality

    basis_a, basis_b, basis_c = basis(lattice_constants)

    # --- Check results

    @test basis_a ≈ [a, 0, 0]
    @test basis_b ≈ [0, b, 0]
    @test basis_c ≈ [0, 0, c]
end

@testset "volume()" begin
    # --- Preparations

    a = 6
    b = 10
    c = 8
    lattice_constants = OrthorhombicLatticeConstants(a, b, c)

    # --- Exercise functionality and check results

    @test volume(lattice_constants) ≈ a * b * c
end

@testset "surface_area()" begin
    # --- Preparations

    a = 6
    b = 10
    c = 8
    lattice_constants = OrthorhombicLatticeConstants(a, b, c)

    # --- Exercise functionality and check results

    @test surface_area(lattice_constants) ≈ 2 * a * b + 2 * b * c + 2 * c * a
end

@testset "reduced_cell()" begin
    # --- Preparations

    a = 2
    b = 3
    c = 5
    lattice_constants = OrthorhombicLatticeConstants(a, b, c)
    basis_a, basis_b, basis_c = basis(lattice_constants)

    # --- Exercise functionality and check results

    # primitive unit cell
    unit_cell = UnitCell(lattice_constants, XtallographyUtils.PRIMITIVE)

    expected_reduced_cell = reduced_cell(
        UnitCell(lattice_constants, XtallographyUtils.PRIMITIVE)
    )

    reduced_cell_ = reduced_cell(unit_cell)
    @test reduced_cell_.lattice_constants isa OrthorhombicLatticeConstants
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
    @test reduced_cell_.lattice_constants isa TriclinicLatticeConstants
    @test volume(reduced_cell_) ≈ 0.5 * volume(unit_cell)
    @test reduced_cell_ ≈ expected_reduced_cell

    # face-centered unit cell
    unit_cell = UnitCell(lattice_constants, XtallographyUtils.FACE)

    expected_reduced_cell = reduced_cell(
        UnitCell(
            LatticeConstants(
                0.5 * (basis_a + basis_b),
                0.5 * (basis_a - basis_b),
                0.5 * (basis_b + basis_c),
            ),
            XtallographyUtils.PRIMITIVE,
        ),
    )

    reduced_cell_ = reduced_cell(unit_cell)
    @test reduced_cell_.lattice_constants isa TriclinicLatticeConstants
    @test volume(reduced_cell_) ≈ 0.25 * volume(unit_cell)
    @test reduced_cell_ ≈ expected_reduced_cell

    # base-centered unit cell
    unit_cell = UnitCell(lattice_constants, XtallographyUtils.BASE)

    expected_reduced_cell = reduced_cell(
        UnitCell(
            LatticeConstants(basis_a, 0.5 * (basis_a + basis_b), basis_c),
            XtallographyUtils.PRIMITIVE,
        ),
    )

    reduced_cell_ = reduced_cell(unit_cell)
    @test reduced_cell_.lattice_constants isa MonoclinicLatticeConstants
    @test volume(reduced_cell_) ≈ 0.5 * volume(unit_cell)
    @test reduced_cell_ ≈ expected_reduced_cell
end

@testset "is_equivalent_unit_cell(::UnitCell, ::UnitCell)" begin
    # --- Preparations

    a = 2
    b = 3
    c = 5
    lattice_constants = OrthorhombicLatticeConstants(a, b, c)
    basis_a, basis_b, basis_c = basis(lattice_constants)

    # --- Exercise functionality and check results

    # equivalent orthorhombic and triclinic unit cells
    orthorhombic_unit_cell = UnitCell(lattice_constants, XtallographyUtils.PRIMITIVE)
    triclinic_unit_cell = UnitCell(
        LatticeConstants(basis_a, basis_b, basis_c; identify_lattice_system=false),
        XtallographyUtils.PRIMITIVE,
    )
    @test is_equivalent_unit_cell(orthorhombic_unit_cell, triclinic_unit_cell)

    # body-centered unit cell
    body_centered_unit_cell = UnitCell(lattice_constants, XtallographyUtils.BODY)
    primitive_unit_cell = UnitCell(
        LatticeConstants(
            basis_a,
            basis_b,
            0.5 * (basis_a + basis_b + basis_c);
            identify_lattice_system=false,
        ),
        XtallographyUtils.PRIMITIVE,
    )
    @test is_equivalent_unit_cell(body_centered_unit_cell, primitive_unit_cell)

    # face-centered unit cell
    face_centered_unit_cell = UnitCell(lattice_constants, XtallographyUtils.FACE)
    primitive_unit_cell = UnitCell(
        LatticeConstants(
            0.5 * (basis_a + basis_b), 0.5 * (basis_a - basis_b), 0.5 * (basis_b + basis_c)
        ),
        XtallographyUtils.PRIMITIVE,
    )
    @test is_equivalent_unit_cell(face_centered_unit_cell, primitive_unit_cell)

    # base-centered unit cell
    base_centered_unit_cell = UnitCell(lattice_constants, XtallographyUtils.BASE)
    primitive_unit_cell = UnitCell(
        LatticeConstants(basis_a, 0.5 * (basis_a + basis_b), basis_c),
        XtallographyUtils.PRIMITIVE,
    )
    @test is_equivalent_unit_cell(base_centered_unit_cell, primitive_unit_cell)
end

@testset "is_equivalent_unit_cell(::LatticeConstants, ::LatticeConstants)" begin
    # --- Preparations

    a_ref = 2
    b_ref = 7
    c_ref = 5
    lattice_constants_ref = OrthorhombicLatticeConstants(a_ref, b_ref, c_ref)

    # --- Exercise functionality and check results

    # unit cells are equivalent
    lattice_constants_test = OrthorhombicLatticeConstants(
        a_ref + 1e-9, b_ref + 3e-8, c_ref - 1e-9
    )
    @test is_equivalent_unit_cell(lattice_constants_test, lattice_constants_ref)

    # test unit cell and reference unit cell are unrelated
    lattice_constants_test = OrthorhombicLatticeConstants(2 * a_ref, b_ref / 2, c_ref)
    @test !is_equivalent_unit_cell(lattice_constants_test, lattice_constants_ref)

    # test unit cell and reference unit cell are for different lattice systems
    lattice_constants_test = CubicLatticeConstants(1)
    @test !is_equivalent_unit_cell(lattice_constants_test, lattice_constants_ref)
end
