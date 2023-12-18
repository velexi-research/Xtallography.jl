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
Tests for methods in lattice/cubic.jl (except for cell standardization methods)
"""
# --- Imports

# Standard library
using Test

# XtallographyUtils package
using XtallographyUtils

# --- Tests

# ------ Types

@testset "CubicLatticeConstants constructor: valid arguments" begin
    # --- Preparations

    a = 1

    # --- Exercise functionality

    lattice_constants = CubicLatticeConstants(a)

    # --- Check results

    @test lattice_constants.a == a
end

@testset "CubicLatticeConstants constructor: invalid arguments" begin
    # --- Tests

    # ------ a

    # a = 0
    local error = nothing
    local error_message = ""
    try
        lattice_constants = CubicLatticeConstants(0)
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
        lattice_constants = CubicLatticeConstants(-1.0)
    catch error
        bt = catch_backtrace()
        error_message = sprint(showerror, error, bt)
    end

    @test error isa ArgumentError

    expected_error = "ArgumentError: `a` must be positive"
    @test startswith(error_message, expected_error)
end

# ------ LatticeConstants functions

@testset "isapprox(::LatticeConstants)" begin
    # --- Preparations

    x = CubicLatticeConstants(1.0)
    y = CubicLatticeConstants(2.0)

    # --- Exercise functionality and check results

    # x ≈ (x + delta)
    @test x ≈ CubicLatticeConstants(1 + 1e-9)

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

    lattice_constants = CubicLatticeConstants(1)
    @test lattice_system(lattice_constants) === Cubic()
end

@testset "standardize()" begin
    # --- Tests

    # ------ Cubic lattices have no lattice constants conventions for primitive, body, and
    #        face centerings

    lattice_constants = CubicLatticeConstants(1.0)

    @test standardize(lattice_constants, Primitive()) == (lattice_constants, Primitive())
    @test standardize(lattice_constants, BodyCentered()) ==
        (lattice_constants, BodyCentered())
    @test standardize(lattice_constants, FaceCentered()) ==
        (lattice_constants, FaceCentered())

    # ------ Invalid centering

    local error = nothing
    local error_message = ""
    try
        standardize(lattice_constants, BaseCentered())
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

# ------ Unit cell computations

@testset "basis()" begin
    # --- Preparations

    a = 5
    lattice_constants = CubicLatticeConstants(a)

    # --- Exercise functionality

    basis_a, basis_b, basis_c = basis(lattice_constants)

    # --- Check results

    @test basis_a ≈ [a, 0, 0]
    @test basis_b ≈ [0, a, 0]
    @test basis_c ≈ [0, 0, a]
end

@testset "volume()" begin
    # --- Preparations

    lattice_constants = CubicLatticeConstants(5)

    # --- Exercise functionality and check results

    @test volume(lattice_constants) ≈ lattice_constants.a^3
end

@testset "surface_area()" begin
    # --- Preparations

    lattice_constants = CubicLatticeConstants(5)

    # --- Exercise functionality and check results

    @test surface_area(lattice_constants) ≈ 6 * lattice_constants.a^2
end

@testset "reduced_cell()" begin
    # --- Preparations

    a = 5
    lattice_constants = CubicLatticeConstants(a)
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
    @test reduced_cell_.lattice_constants isa CubicLatticeConstants
    @test volume(reduced_cell_) ≈ volume(unit_cell)
    @test reduced_cell_ ≈ expected_reduced_cell

    # body-centered unit cell
    unit_cell = UnitCell(lattice_constants, BodyCentered())

    expected_reduced_cell = reduced_cell(
        UnitCell(
            LatticeConstants(basis_a, basis_b, 0.5 * (basis_a + basis_b + basis_c)),
            Primitive(),
        ),
    )

    reduced_cell_ = reduced_cell(unit_cell)
    @test reduced_cell_.lattice_constants isa RhombohedralLatticeConstants
    @test volume(reduced_cell_) ≈ 0.5 * volume(unit_cell)
    @test reduced_cell_ ≈ expected_reduced_cell

    # face-centered unit cell
    unit_cell = UnitCell(lattice_constants, FaceCentered())

    expected_reduced_cell = reduced_cell(
        UnitCell(
            LatticeConstants(
                0.5 * (basis_a + basis_b),
                0.5 * (basis_b - basis_c),
                0.5 * (basis_b + basis_c);
                identify_lattice_system=false,
            ),
            Primitive(),
        ),
    )

    reduced_cell_ = reduced_cell(unit_cell)
    @test reduced_cell_.lattice_constants isa TriclinicLatticeConstants
    @test volume(reduced_cell_) ≈ 0.25 * volume(unit_cell)
    @test reduced_cell_ ≈ expected_reduced_cell
end

@testset "is_equivalent_unit_cell(::UnitCell, ::UnitCell)" begin
    # --- Preparations

    a = 5
    lattice_constants = CubicLatticeConstants(a)
    basis_a, basis_b, basis_c = basis(lattice_constants)

    # --- Tests

    # equivalent cubic and triclinic unit cells
    cubic_unit_cell = UnitCell(lattice_constants, Primitive())
    triclinic_unit_cell = UnitCell(
        LatticeConstants(basis_a, basis_b, basis_c; identify_lattice_system=false),
        Primitive(),
    )
    @test is_equivalent_unit_cell(cubic_unit_cell, triclinic_unit_cell)

    # equivalent body-centered unit cell and primitive rhombohedral unit cell
    body_centered_unit_cell = UnitCell(lattice_constants, BodyCentered())
    primitive_unit_cell = UnitCell(
        LatticeConstants(basis_a, basis_b, 0.5 * (basis_a + basis_b + basis_c)), Primitive()
    )
    @test is_equivalent_unit_cell(body_centered_unit_cell, primitive_unit_cell)

    # equivalent face-centered unit cell and primitive triclinic unit cell
    face_centered_unit_cell = UnitCell(lattice_constants, FaceCentered())
    primitive_unit_cell = UnitCell(
        LatticeConstants(
            0.5 * (basis_a + basis_b), 0.5 * (basis_b - basis_c), 0.5 * (basis_b + basis_c)
        ),
        Primitive(),
    )
    @test is_equivalent_unit_cell(face_centered_unit_cell, primitive_unit_cell)
end

@testset "is_equivalent_unit_cell(::LatticeConstants, ::LatticeConstants)" begin
    # --- Preparations

    lattice_constants_ref = CubicLatticeConstants(2.0)

    # --- Exercise functionality and check results

    # unit cells are equivalent
    lattice_constants_test = CubicLatticeConstants(2.0 + 1e-9)
    @test is_equivalent_unit_cell(lattice_constants_test, lattice_constants_ref)

    # test unit cell is a supercell of the reference unit cell
    lattice_constants_test = CubicLatticeConstants(10.0)
    @test !is_equivalent_unit_cell(lattice_constants_test, lattice_constants_ref)

    # test unit cell and reference unit cell are unrelated
    lattice_constants_test = CubicLatticeConstants(3)
    @test !is_equivalent_unit_cell(lattice_constants_test, lattice_constants_ref)

    # test unit cell and reference unit cell are for different lattice systems
    lattice_constants_test = OrthorhombicLatticeConstants(1, 2, 3)
    @test !is_equivalent_unit_cell(lattice_constants_test, lattice_constants_ref)
end

@testset "is_supercell(): valid arguments" begin
    # --- Preparations

    lattice_constants_ref = CubicLatticeConstants(2.5)

    # --- Exercise functionality and check results

    # ------ default keyword arguments

    # test unit cell is a supercell of referernce unit cell
    lattice_constants_test = CubicLatticeConstants(5.0 + 1e-9)
    @test is_supercell(lattice_constants_test, lattice_constants_ref)

    # test unit cell is not a supercell of referernce unit cell
    lattice_constants_test = CubicLatticeConstants(6.0)
    @test !is_supercell(lattice_constants_test, lattice_constants_ref)

    # ------ keyword arguments tests

    # `tol` large enough for multiplier deviation to pass check
    lattice_constants_test = CubicLatticeConstants(6.0)
    @test is_supercell(lattice_constants_test, lattice_constants_ref; tol=0.45)

    # `tol` large enough for multiplier deviation to pass check but where multiplier
    # is close to 1 so that the test unit cell is not a proper supercell
    lattice_constants_test = CubicLatticeConstants(3.0)
    @test !is_supercell(lattice_constants_test, lattice_constants_ref; tol=0.2)
end

@testset "is_supercell(): invalid arguments" begin
    # --- Preparations

    lattice_constants_ref = CubicLatticeConstants(2.5)
    lattice_constants_test = CubicLatticeConstants(6.0)

    # --- Exercise functionality and check results

    # ------ `tol`

    # tol = 0
    local error, error_message
    try
        is_supercell(lattice_constants_test, lattice_constants_ref; tol=0)
    catch error
        bt = catch_backtrace()
        error_message = sprint(showerror, error, bt)
    end

    @test error isa ArgumentError

    expected_error = "ArgumentError: `tol` must be positive"
    @test startswith(error_message, expected_error)

    # tol < 0
    local error, error_message
    try
        is_supercell(lattice_constants_test, lattice_constants_ref; tol=-0.1)
    catch error
        bt = catch_backtrace()
        error_message = sprint(showerror, error, bt)
    end

    @test error isa ArgumentError

    expected_error = "ArgumentError: `tol` must be positive"
    @test startswith(error_message, expected_error)
end
