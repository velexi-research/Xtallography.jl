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
Tests for methods in lattice/tetragonal.jl (except for cell standardization methods)
"""
# --- Imports

# Standard library
using Test

# XtallographyUtils package
using XtallographyUtils

# --- Tests

# ------ Types

@testset "TetragonalLatticeConstants constructor: valid arguments" begin
    # --- Preparations

    a = 1
    c = 3

    # --- Exercise functionality

    lattice_constants = TetragonalLatticeConstants(a, c)

    # --- Check results

    @test lattice_constants.a == a
    @test lattice_constants.c == c
end

@testset "TetragonalLatticeConstants constructor: invalid arguments" begin
    # --- Preparations

    # Valid arguments
    a = 1
    c = 3

    # --- Tests

    # ------ a

    # a = 0
    local error = nothing
    local error_message = ""
    try
        lattice_constants = TetragonalLatticeConstants(0, c)
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
        lattice_constants = TetragonalLatticeConstants(-1.0, c)
    catch error
        bt = catch_backtrace()
        error_message = sprint(showerror, error, bt)
    end

    @test error isa ArgumentError

    expected_error = "ArgumentError: `a` must be positive"
    @test startswith(error_message, expected_error)

    # ------ c

    # c = 0
    local error = nothing
    local error_message = ""
    try
        lattice_constants = TetragonalLatticeConstants(a, 0)
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
        lattice_constants = TetragonalLatticeConstants(a, -1.0)
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

    x = TetragonalLatticeConstants(1.0, 2.0)
    y = TetragonalLatticeConstants(1.5, 2.5)

    # --- Exercise functionality and check results

    # x ≈ (x + delta)
    @test x ≈ TetragonalLatticeConstants(1.0 + 1e-9, 2.0)
    @test x ≈ TetragonalLatticeConstants(1.0, 2.0 + 1e-9)

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

    lattice_constants = TetragonalLatticeConstants(1, 2)
    @test lattice_system(lattice_constants) == Tetragonal
end

@testset "standardize()" begin
    # --- Tests

    # ------ Tetragonal lattices have no lattice constants conventions for primitive and
    #        body centerings

    a = 1.0
    c = 10.0
    lattice_constants = TetragonalLatticeConstants(a, c)

    # centering = primitive
    standardized_lattice_constants, standardized_centering = standardize(
        lattice_constants, Primitive()
    )

    expected_lattice_constants = lattice_constants
    @test standardized_lattice_constants ≈ expected_lattice_constants

    @test standardized_centering == Primitive()

    # centering = body-centered
    standardized_lattice_constants, standardized_centering = standardize(
        lattice_constants, BodyCentered()
    )

    expected_lattice_constants = lattice_constants
    @test standardized_lattice_constants ≈ expected_lattice_constants

    @test standardized_centering == BodyCentered()

    # ------ Invalid centerings

    for centering in (FaceCentered(), BaseCentered())
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

# ------ Unit cell computations

@testset "basis()" begin
    # --- Preparations

    a = 5
    c = 7
    lattice_constants = TetragonalLatticeConstants(a, c)

    # --- Exercise functionality

    basis_a, basis_b, basis_c = basis(lattice_constants)

    # --- Check results

    @test basis_a ≈ [a, 0, 0]
    @test basis_b ≈ [0, a, 0]
    @test basis_c ≈ [0, 0, c]
end

@testset "volume()" begin
    # --- Preparations

    lattice_constants = TetragonalLatticeConstants(5, 7)

    # --- Exercise functionality and check results

    @test volume(lattice_constants) ≈ lattice_constants.a^2 * lattice_constants.c
end

@testset "surface_area()" begin
    # --- Preparations

    lattice_constants = TetragonalLatticeConstants(5, 7)

    # --- Exercise functionality and check results

    @test surface_area(lattice_constants) ≈
        2 * lattice_constants.a^2 + 4 * lattice_constants.a * lattice_constants.c
end

@testset "reduced_cell()" begin
    # --- Preparations

    a = 5
    c = 7
    lattice_constants = TetragonalLatticeConstants(a, c)
    basis_a, basis_b, basis_c = basis(lattice_constants)

    # --- Exercise functionality and check results

    # primitive unit cell
    unit_cell = UnitCell(lattice_constants, Primitive())

    expected_reduced_cell = reduced_cell(
        UnitCell(
            LatticeConstants(basis_a, basis_b, basis_c; identify_lattice_system=false),
            Primitive(),
        ),
    )

    reduced_cell_ = reduced_cell(unit_cell)
    @test reduced_cell_.lattice_constants isa TetragonalLatticeConstants
    @test volume(reduced_cell_) ≈ volume(unit_cell)
    @test reduced_cell_ ≈ expected_reduced_cell

    # body-centered unit cell
    unit_cell = UnitCell(lattice_constants, BodyCentered())

    expected_reduced_cell = reduced_cell(
        UnitCell(
            LatticeConstants(
                basis_a,
                basis_b,
                0.5 * (basis_a + basis_b + basis_c);
                identify_lattice_system=false,
            ),
            Primitive(),
        ),
    )

    reduced_cell_ = reduced_cell(unit_cell)
    @test reduced_cell_.lattice_constants isa TriclinicLatticeConstants
    @test volume(reduced_cell_) ≈ 0.5 * volume(unit_cell)
    @test reduced_cell_ ≈ expected_reduced_cell

    # face-centered unit cell
    unit_cell = UnitCell(lattice_constants, FaceCentered())

    expected_reduced_cell = reduced_cell(
        UnitCell(TetragonalLatticeConstants(a / sqrt(2), c), BodyCentered())
    )

    reduced_cell_ = reduced_cell(unit_cell)
    @test reduced_cell_.lattice_constants isa TriclinicLatticeConstants
    @test volume(reduced_cell_) ≈ 0.25 * volume(unit_cell)
    @test reduced_cell_ ≈ expected_reduced_cell
end

@testset "is_equivalent_unit_cell(::UnitCell, ::UnitCell)" begin
    # --- Preparations

    a = 2
    c = 5
    lattice_constants = TetragonalLatticeConstants(a, c)
    basis_a, basis_b, basis_c = basis(lattice_constants)

    # --- Exercise functionality and check results

    # equivalent tetragonal and triclinic unit cells
    tetragonal_unit_cell = UnitCell(lattice_constants, Primitive())
    triclinic_unit_cell = UnitCell(
        LatticeConstants(basis_a, basis_b, basis_c; identify_lattice_system=false),
        Primitive(),
    )
    @test is_equivalent_unit_cell(tetragonal_unit_cell, triclinic_unit_cell)

    # body-centered unit cell
    body_centered_unit_cell = UnitCell(lattice_constants, BodyCentered())
    primitive_unit_cell = UnitCell(
        LatticeConstants(
            basis_a,
            basis_b,
            0.5 * (basis_a + basis_b + basis_c);
            identify_lattice_system=false,
        ),
        Primitive(),
    )
    @test is_equivalent_unit_cell(body_centered_unit_cell, primitive_unit_cell)

    # face-centered unit cell
    face_centered_unit_cell = UnitCell(lattice_constants, FaceCentered())
    primitive_unit_cell = UnitCell(
        TetragonalLatticeConstants(a / sqrt(2), c), BodyCentered()
    )
    @test is_equivalent_unit_cell(face_centered_unit_cell, primitive_unit_cell)
end

@testset "is_equivalent_unit_cell(::LatticeConstants, ::LatticeConstants)" begin
    # --- Preparations

    a_ref = 2
    c_ref = 5
    lattice_constants_ref = TetragonalLatticeConstants(a_ref, c_ref)

    # --- Exercise functionality and check results

    # unit cells are equivalent
    lattice_constants_test = TetragonalLatticeConstants(a_ref + 1e-9, c_ref - 1e-9)
    @test is_equivalent_unit_cell(lattice_constants_test, lattice_constants_ref)

    # test unit cell and reference unit cell are unrelated
    lattice_constants_test = TetragonalLatticeConstants(2 * a_ref, c_ref)
    @test !is_equivalent_unit_cell(lattice_constants_test, lattice_constants_ref)

    # test unit cell and reference unit cell are for different lattice systems
    lattice_constants_test = CubicLatticeConstants(1)
    @test !is_equivalent_unit_cell(lattice_constants_test, lattice_constants_ref)
end
