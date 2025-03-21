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

# Xtallography package
using Xtallography

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
    expected_message = "`a` must be positive"
    @test_throws DomainError(0, expected_message) CubicLatticeConstants(0)

    # a < 0
    expected_message = "`a` must be positive"
    @test_throws DomainError(-1, expected_message) CubicLatticeConstants(-1.0)
end

@testset "CubicLatticeConstantDeltas constructor" begin
    # --- Tests

    Δa = 1
    Δlattice_constants = CubicLatticeConstantDeltas(Δa)

    @test Δlattice_constants.Δa == Δa
end

# ------ LatticeConstants functions

@testset "isapprox(::CubicLatticeConstants)" begin
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

@testset "convert(::CubicLatticeConstants)" begin
    # --- Tests

    x = CubicLatticeConstants(rand())

    x_vector = convert(Vector, x)
    @test x_vector == [x.a]
end

@testset "-(::CubicLatticeConstants)" begin
    # --- Tests

    x = CubicLatticeConstants(1)
    y = CubicLatticeConstants(2)
    @test x - y == CubicLatticeConstantDeltas(x.a - y.a)
end

@testset "lattice_system(::CubicLatticeConstants)" begin
    # --- Tests

    lattice_constants = CubicLatticeConstants(1)
    @test lattice_system(lattice_constants) === cubic
end

@testset "standardize(): cubic" begin
    # --- Tests

    # ------ Cubic lattices have no lattice constants conventions for primitive, body, and
    #        face centerings

    lattice_constants = CubicLatticeConstants(1.0)

    @test standardize(lattice_constants, primitive) == (lattice_constants, primitive)
    @test standardize(lattice_constants, body_centered) ==
        (lattice_constants, body_centered)
    @test standardize(lattice_constants, face_centered) ==
        (lattice_constants, face_centered)

    # ------ Invalid centering

    expected_message = (
        "Invalid Bravais lattice: (lattice_system=Cubic, centering=BaseCentered)"
    )
    @test_throws ArgumentError(expected_message) standardize(
        lattice_constants, base_centered
    )
end

# ------ LatticeConstantDeltas functions

@testset "isapprox(::CubicLatticeConstantDeltas)" begin
    # --- Preparations

    x = CubicLatticeConstantDeltas(1.0)
    y = CubicLatticeConstantDeltas(2.0)

    # --- Exercise functionality and check results

    # x ≈ (x + delta)
    @test x ≈ CubicLatticeConstantDeltas(1 + 1e-9)

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

@testset "convert(::CubicLatticeConstantDeltas)" begin
    # --- Tests

    x = CubicLatticeConstantDeltas(rand())

    x_vector = convert(Vector, x)
    @test x_vector == [x.Δa]
end

@testset "lattice_system(::CubicLatticeConstantDeltas)" begin
    # --- Tests

    Δlattice_constants = CubicLatticeConstantDeltas(1)
    @test lattice_system(Δlattice_constants) === cubic
end

# ------ Unit cell computations

@testset "basis(::CubicLatticeConstants)" begin
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

@testset "volume(::CubicLatticeConstants)" begin
    # --- Preparations

    lattice_constants = CubicLatticeConstants(5)

    # --- Exercise functionality and check results

    @test volume(lattice_constants) ≈ lattice_constants.a^3
end

@testset "surface_area(::CubicLatticeConstants)" begin
    # --- Preparations

    lattice_constants = CubicLatticeConstants(5)

    # --- Exercise functionality and check results

    @test surface_area(lattice_constants) ≈ 6 * lattice_constants.a^2
end

@testset "reduced_cell(): cubic" begin
    # --- Preparations

    a = 5
    lattice_constants = CubicLatticeConstants(a)
    basis_a, basis_b, basis_c = basis(lattice_constants)

    # --- Exercise functionality and check results

    # primitive unit cell
    unit_cell = UnitCell(lattice_constants, primitive)

    expected_reduced_cell = reduced_cell(
        UnitCell(
            LatticeConstants(basis_a, basis_b, basis_c; identify_lattice_system=false),
            primitive,
        ),
    )

    reduced_cell_ = reduced_cell(unit_cell)
    @test reduced_cell_.lattice_constants isa CubicLatticeConstants
    @test volume(reduced_cell_) ≈ volume(unit_cell)
    @test reduced_cell_ ≈ expected_reduced_cell

    # body-centered unit cell
    unit_cell = UnitCell(lattice_constants, body_centered)

    expected_reduced_cell = reduced_cell(
        UnitCell(
            LatticeConstants(basis_a, basis_b, 0.5 * (basis_a + basis_b + basis_c)),
            primitive,
        ),
    )

    reduced_cell_ = reduced_cell(unit_cell)
    @test reduced_cell_.lattice_constants isa RhombohedralLatticeConstants
    @test volume(reduced_cell_) ≈ 0.5 * volume(unit_cell)
    @test reduced_cell_ ≈ expected_reduced_cell

    # face-centered unit cell
    unit_cell = UnitCell(lattice_constants, face_centered)

    expected_reduced_cell = reduced_cell(
        UnitCell(
            LatticeConstants(
                0.5 * (basis_a + basis_b),
                0.5 * (basis_b - basis_c),
                0.5 * (basis_b + basis_c);
                identify_lattice_system=false,
            ),
            primitive,
        ),
    )

    reduced_cell_ = reduced_cell(unit_cell)
    @test reduced_cell_.lattice_constants isa TriclinicLatticeConstants
    @test volume(reduced_cell_) ≈ 0.25 * volume(unit_cell)
    @test reduced_cell_ ≈ expected_reduced_cell
end

@testset "is_equivalent_unit_cell(::UnitCell): cubic" begin
    # --- Preparations

    a = 5
    lattice_constants = CubicLatticeConstants(a)
    basis_a, basis_b, basis_c = basis(lattice_constants)

    # --- Tests

    # equivalent cubic and triclinic unit cells
    cubic_unit_cell = UnitCell(lattice_constants, primitive)
    triclinic_unit_cell = UnitCell(
        LatticeConstants(basis_a, basis_b, basis_c; identify_lattice_system=false),
        primitive,
    )
    @test is_equivalent_unit_cell(cubic_unit_cell, triclinic_unit_cell)

    # equivalent body-centered unit cell and primitive rhombohedral unit cell
    body_centered_unit_cell = UnitCell(lattice_constants, body_centered)
    primitive_unit_cell = UnitCell(
        LatticeConstants(basis_a, basis_b, 0.5 * (basis_a + basis_b + basis_c)), primitive
    )
    @test is_equivalent_unit_cell(body_centered_unit_cell, primitive_unit_cell)

    # equivalent face-centered unit cell and primitive triclinic unit cell
    face_centered_unit_cell = UnitCell(lattice_constants, face_centered)
    primitive_unit_cell = UnitCell(
        LatticeConstants(
            0.5 * (basis_a + basis_b), 0.5 * (basis_b - basis_c), 0.5 * (basis_b + basis_c)
        ),
        primitive,
    )
    @test is_equivalent_unit_cell(face_centered_unit_cell, primitive_unit_cell)
end

@testset "is_equivalent_unit_cell(::CubicLatticeConstants)" begin
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

@testset "is_supercell(::CubicLatticeConstants): valid arguments" begin
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

@testset "is_supercell(::CubicLatticeConstants): invalid arguments" begin
    # --- Preparations

    lattice_constants_ref = CubicLatticeConstants(2.5)
    lattice_constants_test = CubicLatticeConstants(6.0)

    # --- Exercise functionality and check results

    # ------ `tol`

    # tol = 0
    expected_message = "`tol` must be positive"
    @test_throws DomainError(0, expected_message) is_supercell(
        lattice_constants_test, lattice_constants_ref; tol=0
    )

    # tol < 0
    expected_message = "`tol` must be positive"
    @test_throws DomainError(-0.1, expected_message) is_supercell(
        lattice_constants_test, lattice_constants_ref; tol=-0.1
    )
end
