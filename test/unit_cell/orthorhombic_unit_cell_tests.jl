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
Tests for OrthorhombicUnitCell types and methods in unit_cell/orthorhombic.jl (except for
conventional_cell() tests)
"""
# --- Imports

# Standard library
using Test
using LinearAlgebra: qr, I

# Xtallography package
using Xtallography

# Testing utilities
include("testing_utilities.jl")

# --- Tests

# ------ Types

@testset "OrthorhombicUnitCell(::NamedTuple,::UnitCellSymmetry) inner constructor" begin
    # --- Tests

    # Valid arguments
    lattice_constants_ = (a=1, b=3, c=5)
    symmetry_ = primitive_unit_cell_symmetry
    unit_cell = OrthorhombicUnitCell(lattice_constants_, symmetry_)

    @test lattice_constants(unit_cell) == lattice_constants_
    @test symmetry(unit_cell) == symmetry_

    # Invalid arguments
    lattice_constants_ = (a=1,)
    expected_message = (
        "Invalid lattice_constants argument passed to UnitCell{Orthorhombic} " *
        "constructor. Expected keys: (:a, :b, :c). " *
        "Provided keys: $(keys(lattice_constants_))."
    )
    @test_throws ArgumentError(expected_message) OrthorhombicUnitCell(
        lattice_constants_, symmetry_
    )
end

@testset "OrthorhombicUnitCell(::NamedTuple;::Centering,::Set) outer constructor" begin
    # --- Preparations

    lattice_constants_ = (a=1, b=3, c=5)

    # --- Tests

    # Default keyword arguments
    unit_cell = OrthorhombicUnitCell(lattice_constants_)

    @test lattice_constants(unit_cell) == lattice_constants_
    @test symmetry(unit_cell) == primitive_unit_cell_symmetry

    # Non-default centering keyword argument
    for centering_ in CENTERINGS
        unit_cell = OrthorhombicUnitCell(lattice_constants_; centering=centering_)

        @test lattice_constants(unit_cell) == lattice_constants_
        @test symmetry(unit_cell) == UnitCellSymmetry(centering_, Set{SymmetryElement}())
    end

    # Non-default symmetry_elements keyword argument
    symmetry_elements_ = [a_4_2]
    unit_cell = OrthorhombicUnitCell(
        lattice_constants_; symmetry_elements=symmetry_elements_
    )

    @test lattice_constants(unit_cell) == lattice_constants_
    @test symmetry(unit_cell) ==
        UnitCellSymmetry(primitive_centering, Set{SymmetryElement}(symmetry_elements_))
end

@testset "OrthorhombicUnitCell(::Real,::Real,::Real;::Centering,::Set) outer constructor: valid arguments" begin

    # --- Preparations

    a = 1
    b = 2
    c = 3

    # --- Exercise functionality and check results

    # Default keyword arguments
    unit_cell = OrthorhombicUnitCell(a, b, c)

    @test lattice_constants(unit_cell).a == a
    @test lattice_constants(unit_cell).b == b
    @test lattice_constants(unit_cell).c == c

    @test centering(unit_cell) === primitive_centering
    @test symmetry_elements(unit_cell) isa Set
    @test isempty(symmetry_elements(unit_cell))

    # Non-default centering keyword argument
    for centering_ in CENTERINGS
        unit_cell = OrthorhombicUnitCell(a, b, c; centering=centering_)

        @test lattice_constants(unit_cell).a == a
        @test lattice_constants(unit_cell).b == b
        @test lattice_constants(unit_cell).c == c
        @test centering(unit_cell) === centering_
        @test symmetry_elements(unit_cell) isa Set
        @test isempty(symmetry_elements(unit_cell))
    end

    # Non-default symmetry_elements keyword argument
    symmetry_elements_ = [a_4_2]
    unit_cell = OrthorhombicUnitCell(a, b, c; symmetry_elements=symmetry_elements_)

    @test lattice_constants(unit_cell).a == a
    @test lattice_constants(unit_cell).b == b
    @test lattice_constants(unit_cell).c == c
    @test centering(unit_cell) === primitive_centering
    @test symmetry_elements(unit_cell) isa Set
    @test symmetry_elements(unit_cell) == Set(symmetry_elements_)
end

@testset "OrthorhombicUnitCell(::Real,::Real,::Real;::Centering,::Set) outer constructor: invalid arguments" begin

    # --- Preparations

    # Valid arguments
    a = 1
    b = 2
    c = 3

    # --- Tests

    # ------ a

    # a = 0
    expected_message = "`a` must be positive"
    @test_throws DomainError(0, expected_message) OrthorhombicUnitCell(0, b, c)

    # a < 0
    expected_message = "`a` must be positive"
    @test_throws DomainError(-1, expected_message) OrthorhombicUnitCell(-1, b, c)

    # ------ b

    # b = 0
    expected_message = "`b` must be positive"
    @test_throws DomainError(0, expected_message) OrthorhombicUnitCell(a, 0, c)

    # b < 0
    expected_message = "`b` must be positive"
    @test_throws DomainError(-2.0, expected_message) OrthorhombicUnitCell(a, -2.0, c)

    # ------ c

    # c = 0
    expected_message = "`c` must be positive"
    @test_throws DomainError(0, expected_message) OrthorhombicUnitCell(a, b, 0.0)

    # c < 0
    expected_message = "`c` must be positive"
    @test_throws DomainError(-3, expected_message) OrthorhombicUnitCell(a, b, -3.0)

    # ------ symmetry elements

    # TODO
end

@testset "UnitCell(::UnitCell) copy constructor: orthorhombic" begin
    # --- Preparations

    lattice_constants_ = (a=1, b=2, c=3)
    centering = primitive_centering
    symmetry_elements_ = [a_4_2]
    reference_unit_cell = OrthorhombicUnitCell(
        lattice_constants_,
        UnitCellSymmetry(centering; symmetry_elements=symmetry_elements_),
    )

    # --- Tests

    unit_cell = UnitCell(reference_unit_cell)

    @test unit_cell === reference_unit_cell
end

@testset "UnitCell(::Vector, ::Vector, ::Vector;::Bool,::Centering) constructor: orthorhombic" begin

    # --- Preparations

    # Generate random rotation matrix
    rotations = [I, qr(rand(3)).Q, qr(rand(3)).Q]

    # --- Tests

    a = 5
    b = 8
    c = 10

    basis_a = [a, 0, 0]
    basis_b = [0, b, 0]
    basis_c = [0, 0, c]

    # default centering keyword argument
    expected_unit_cell = OrthorhombicUnitCell(a, b, c)
    test_basis_rotations_and_permutations(
        rotations, expected_unit_cell, basis_a, basis_b, basis_c
    )

    # centering = primitive-centering, body-centering, or face-centering
    for centering_ in (primitive_centering, body_centering, face_centering)
        expected_unit_cell = OrthorhombicUnitCell(a, b, c; centering=centering_)
        test_basis_rotations_and_permutations(
            rotations, expected_unit_cell, basis_a, basis_b, basis_c; centering=centering_
        )
    end

    # centering = base-centering; a, b < c
    a = 5
    b = 8
    c = 10

    basis_a = [a, 0, 0]
    basis_b = [0, b, 0]
    basis_c = [0, 0, c]

    expected_unit_cell = OrthorhombicUnitCell(a, b, c; centering=base_centering)

    test_basis_rotations_and_permutations(
        rotations, expected_unit_cell, basis_a, basis_b, basis_c; centering=base_centering
    )

    # centering = base-centering; a < c < b
    a = 5
    b = 10
    c = 8

    basis_a = [a, 0, 0]
    basis_b = [0, b, 0]
    basis_c = [0, 0, c]

    expected_unit_cell = OrthorhombicUnitCell(a, b, c; centering=base_centering)
    test_basis_rotations_and_permutations(
        rotations, expected_unit_cell, basis_a, basis_b, basis_c; centering=base_centering
    )

    # centering = base-centering; b < c < a
    a = 10
    b = 5
    c = 8

    basis_a = [a, 0, 0]
    basis_b = [0, b, 0]
    basis_c = [0, 0, c]

    expected_unit_cell = OrthorhombicUnitCell(b, a, c; centering=base_centering)
    test_basis_rotations_and_permutations(
        rotations, expected_unit_cell, basis_a, basis_b, basis_c; centering=base_centering
    )
end

# ------ Methods

@testset "lattice_system(::OrthorhombicUnitCell)" begin
    lattice_constants_ = (a=1, b=2, c=3)
    unit_cell = OrthorhombicUnitCell(lattice_constants_)
    @test lattice_system(unit_cell) === orthorhombic
end

@testset "lattice_constants(::OrthorhombicUnitCell)" begin
    lattice_constants_ = (a=1, b=2, c=3)
    unit_cell = OrthorhombicUnitCell(lattice_constants_)
    @test lattice_constants(unit_cell) == lattice_constants_
end

@testset "symmetry(::OrthorhombicUnitCell)" begin
    # --- Preparations

    lattice_constants_ = (a=1, b=2, c=3)
    symmetry_elements_ = Set([a_2_1, d_perp_110])

    # --- Tests

    for centering_ in CENTERINGS
        symmetry_ = UnitCellSymmetry(;
            centering=centering_, symmetry_elements=symmetry_elements_
        )
        unit_cell = OrthorhombicUnitCell(lattice_constants_, symmetry_)

        @test symmetry(unit_cell) == symmetry_
    end
end

@testset "centering(::OrthorhombicUnitCell)" begin
    # --- Preparations

    lattice_constants_ = (a=1, b=2, c=3)

    # --- Tests

    # Default centering
    unit_cell = OrthorhombicUnitCell(lattice_constants_)
    @test centering(unit_cell) == primitive_centering

    # Non-default centering
    for centering_ in CENTERINGS
        unit_cell = OrthorhombicUnitCell(lattice_constants_; centering=centering_)
        @test centering(unit_cell) === centering_
    end
end

@testset "symmetry_elements(::OrthorhombicUnitCell)" begin
    # --- Preparations

    lattice_constants_ = (a=1, b=2, c=3)

    # --- Tests

    # Default symmetry elements
    unit_cell = OrthorhombicUnitCell(lattice_constants_)
    @test symmetry_elements(unit_cell) isa Set{SymmetryElement}
    @test isempty(symmetry_elements(unit_cell))

    # Non-default symmetry elements
    symmetry_elements_ = Set([a_2_1, d_perp_110])
    unit_cell = OrthorhombicUnitCell(
        lattice_constants_; symmetry_elements=symmetry_elements_
    )
    @test symmetry_elements(unit_cell) == symmetry_elements_
end

@testset "is_bravais_lattice(::OrthorhombicUnitCell)" begin
    # --- Preparations

    a = 1
    b = 2
    c = 3

    # --- Tests

    # Valid Bravais lattices
    for centering_ in CENTERINGS
        unit_cell = OrthorhombicUnitCell(a, b, c; centering=centering_)
        @test is_bravais_lattice(unit_cell)
    end
end

@testset "basis(::OrthorhombicUnitCell)" begin
    # --- Preparations

    a = 2
    b = 3
    c = 5
    lattice_constants = OrthorhombicUnitCell(a, b, c)

    # --- Exercise functionality

    basis_a, basis_b, basis_c = basis(lattice_constants)

    # --- Check results

    @test basis_a ≈ [a, 0, 0]
    @test basis_b ≈ [0, b, 0]
    @test basis_c ≈ [0, 0, c]
end

@testset "volume(::OrthorhombicUnitCell)" begin
    # --- Preparations

    a = 6
    b = 10
    c = 8
    unit_cell = OrthorhombicUnitCell(a, b, c)

    # --- Exercise functionality and check results

    @test volume(unit_cell) ≈ a * b * c
end

@testset "surface_area(::OrthorhombicUnitCell)" begin
    # --- Preparations

    a = 6
    b = 10
    c = 8
    unit_cell = OrthorhombicUnitCell(a, b, c)

    # --- Exercise functionality and check results

    @test surface_area(unit_cell) ≈ 2 * a * b + 2 * b * c + 2 * c * a
end

@testset "standardize(::OrthorhombicUnitCell): primitive, body-centered, face-centered" begin
    # --- Tests

    # ------ lattice constants already in standard form

    a = 1.0
    b = 5.0
    c = 10.0
    reference_unit_cell = OrthorhombicUnitCell(a, b, c)

    # centering = primitive_centering
    standardized_unit_cell = standardize(reference_unit_cell)

    expected_unit_cell = OrthorhombicUnitCell(a, b, c; centering=primitive_centering)
    @test standardized_unit_cell ≈ expected_unit_cell

    # centering = body_centering
    standardized_unit_cell = standardize(
        OrthorhombicUnitCell(
            lattice_constants(reference_unit_cell); centering=body_centering
        ),
    )

    expected_unit_cell = OrthorhombicUnitCell(a, b, c; centering=body_centering)
    @test standardized_unit_cell ≈ expected_unit_cell

    # centering = face_centering
    standardized_unit_cell = standardize(
        OrthorhombicUnitCell(
            lattice_constants(reference_unit_cell); centering=face_centering
        ),
    )

    expected_unit_cell = OrthorhombicUnitCell(a, b, c; centering=face_centering)
    @test standardized_unit_cell ≈ expected_unit_cell

    # ------ lattice constants not sorted

    a = 1.0
    b = 10.0
    c = 5.0
    reference_unit_cell = OrthorhombicUnitCell(a, b, c)

    # centering = primitive_centering
    standardized_unit_cell = standardize(reference_unit_cell)

    expected_unit_cell = OrthorhombicUnitCell(a, c, b; centering=primitive_centering)
    @test standardized_unit_cell ≈ expected_unit_cell

    # centering = body_centering
    standardized_unit_cell = standardize(
        OrthorhombicUnitCell(
            lattice_constants(reference_unit_cell); centering=body_centering
        ),
    )

    expected_unit_cell = OrthorhombicUnitCell(a, c, b; centering=body_centering)
    @test standardized_unit_cell ≈ expected_unit_cell

    # centering = face_centering
    standardized_unit_cell = standardize(
        OrthorhombicUnitCell(
            lattice_constants(reference_unit_cell); centering=face_centering
        ),
    )

    expected_unit_cell = OrthorhombicUnitCell(a, c, b; centering=face_centering)
    @test standardized_unit_cell ≈ expected_unit_cell
end

@testset "standardize(::OrthorhombicUnitCell): base-centered" begin
    # --- Tests

    # ------ lattice constants already in standard form

    a = 1.0
    b = 5.0
    c = 10.0
    unit_cell = OrthorhombicUnitCell(a, b, c; centering=base_centering)

    standardized_unit_cell = standardize(unit_cell)

    expected_unit_cell = OrthorhombicUnitCell(a, b, c; centering=base_centering)
    @test standardized_unit_cell ≈ expected_unit_cell

    # ------ lattice constants not sorted

    a = 5.0
    b = 1.0
    c = 10.0
    unit_cell = OrthorhombicUnitCell(a, b, c; centering=base_centering)

    standardized_unit_cell = standardize(unit_cell)

    expected_unit_cell = OrthorhombicUnitCell(b, a, c; centering=base_centering)
    @test standardized_unit_cell ≈ expected_unit_cell
end

@testset "reduced_cell(::OrthorhombicUnitCell)" begin
    # --- Preparations

    a = 2
    b = 3
    c = 5
    reference_unit_cell = OrthorhombicUnitCell(a, b, c)
    basis_a, basis_b, basis_c = basis(reference_unit_cell)

    # --- Exercise functionality and check results

    # primitive unit cell
    unit_cell = OrthorhombicUnitCell(
        lattice_constants(reference_unit_cell); centering=primitive_centering
    )

    expected_reduced_cell = reduced_cell(
        UnitCell(
            basis_a,
            basis_b,
            basis_c;
            identify_lattice_system=false,
            centering=primitive_centering,
        ),
    )

    reduced_cell_ = reduced_cell(unit_cell)
    @test reduced_cell_ isa OrthorhombicUnitCell
    @test volume(reduced_cell_) ≈ volume(reference_unit_cell)
    @test reduced_cell_ ≈ expected_reduced_cell

    # body-centered unit cell
    unit_cell = OrthorhombicUnitCell(
        lattice_constants(reference_unit_cell); centering=body_centering
    )

    expected_reduced_cell = reduced_cell(
        UnitCell(
            basis_a,
            basis_b,
            0.5 * (basis_a + basis_b + basis_c);
            identify_lattice_system=false,
            centering=primitive_centering,
        ),
    )

    reduced_cell_ = reduced_cell(unit_cell)
    @test reduced_cell_ isa TriclinicUnitCell
    @test volume(reduced_cell_) ≈ 0.5 * volume(reference_unit_cell)
    @test reduced_cell_ ≈ expected_reduced_cell

    # face-centered unit cell
    unit_cell = OrthorhombicUnitCell(
        lattice_constants(reference_unit_cell); centering=face_centering
    )

    expected_reduced_cell = reduced_cell(
        UnitCell(
            0.5 * (basis_a + basis_b),
            0.5 * (basis_a - basis_b),
            0.5 * (basis_b + basis_c);
            identify_lattice_system=false,
            centering=primitive_centering,
        ),
    )

    reduced_cell_ = reduced_cell(unit_cell)
    @test reduced_cell_ isa TriclinicUnitCell
    @test volume(reduced_cell_) ≈ 0.25 * volume(reference_unit_cell)
    @test reduced_cell_ ≈ expected_reduced_cell

    # base-centered unit cell
    unit_cell = OrthorhombicUnitCell(
        lattice_constants(reference_unit_cell); centering=base_centering
    )

    expected_reduced_cell = reduced_cell(
        UnitCell(
            basis_a,
            0.5 * (basis_a + basis_b),
            basis_c;
            identify_lattice_system=false,
            centering=primitive_centering,
        ),
    )

    reduced_cell_ = reduced_cell(unit_cell)
    @test reduced_cell_ isa MonoclinicUnitCell
    @test volume(reduced_cell_) ≈ 0.5 * volume(reference_unit_cell)
    @test reduced_cell_ ≈ expected_reduced_cell
end

@testset "is_equivalent_unit_cell(::UnitCell): orthorhombic" begin
    # --- Preparations

    a = 2
    b = 3
    c = 5
    reference_unit_cell = OrthorhombicUnitCell(a, b, c)
    basis_a, basis_b, basis_c = basis(reference_unit_cell)

    # --- Exercise functionality and check results

    # equivalent orthorhombic and triclinic unit cells
    orthorhombic_unit_cell = reference_unit_cell
    triclinic_unit_cell = UnitCell(
        basis_a,
        basis_b,
        basis_c;
        identify_lattice_system=false,
        centering=primitive_centering,
    )
    @test is_equivalent_unit_cell(orthorhombic_unit_cell, triclinic_unit_cell)

    # body-centered unit cell
    body_centering_unit_cell = OrthorhombicUnitCell(
        lattice_constants(reference_unit_cell); centering=body_centering
    )
    primitive_unit_cell = UnitCell(
        basis_a,
        basis_b,
        0.5 * (basis_a + basis_b + basis_c);
        identify_lattice_system=false,
        centering=primitive_centering,
    )
    @test is_equivalent_unit_cell(body_centering_unit_cell, primitive_unit_cell)

    # face-centered unit cell
    face_centering_unit_cell = OrthorhombicUnitCell(
        lattice_constants(reference_unit_cell); centering=face_centering
    )
    primitive_unit_cell = UnitCell(
        0.5 * (basis_a + basis_b),
        0.5 * (basis_a - basis_b),
        0.5 * (basis_b + basis_c);
        identify_lattice_system=false,
        centering=primitive_centering,
    )
    @test is_equivalent_unit_cell(face_centering_unit_cell, primitive_unit_cell)

    # base-centered unit cell
    base_centering_unit_cell = OrthorhombicUnitCell(
        lattice_constants(reference_unit_cell); centering=base_centering
    )
    primitive_unit_cell = UnitCell(
        basis_a,
        0.5 * (basis_a + basis_b),
        basis_c;
        identify_lattice_system=false,
        centering=primitive_centering,
    )
    @test is_equivalent_unit_cell(base_centering_unit_cell, primitive_unit_cell)
end

@testset "is_equivalent_unit_cell(::OrthorhombicUnitCell)" begin
    # --- Preparations

    a_ref = 2
    b_ref = 7
    c_ref = 5
    lattice_constants_ref = OrthorhombicUnitCell(a_ref, b_ref, c_ref)

    # --- Exercise functionality and check results

    # unit cells are equivalent
    lattice_constants_test = OrthorhombicUnitCell(a_ref + 1e-9, b_ref + 3e-8, c_ref - 1e-9)
    @test is_equivalent_unit_cell(lattice_constants_test, lattice_constants_ref)

    # test unit cell and reference unit cell are unrelated
    lattice_constants_test = OrthorhombicUnitCell(2 * a_ref, b_ref / 2, c_ref)
    @test !is_equivalent_unit_cell(lattice_constants_test, lattice_constants_ref)

    # test unit cell and reference unit cell are for different lattice systems
    lattice_constants_test = CubicUnitCell(1)
    @test !is_equivalent_unit_cell(lattice_constants_test, lattice_constants_ref)
end

@testset ":(==)(::OrthorhombicUnitCell)" begin
    # --- Exercise functionality and check results

    # unit cells are equal
    @test OrthorhombicUnitCell(1, 2, 3) == OrthorhombicUnitCell(1, 2, 3)

    # lattice constants are not equal
    @test OrthorhombicUnitCell(1, 2, 3) != OrthorhombicUnitCell(1, 2, 4)

    # centerings are not equal
    @test (
        OrthorhombicUnitCell(1, 2, 3; centering=primitive_centering) !=
        OrthorhombicUnitCell(1, 2, 3; centering=body_centering)
    )

    # symmetry elements are not equal
    @test (
        OrthorhombicUnitCell(1, 2, 3; symmetry_elements=Set()) !=
        OrthorhombicUnitCell(1, 2, 3; symmetry_elements=Set([a_2_1, d_perp_110]))
    )
end

@testset "isapprox(::OrthorhombicUnitCell)" begin
    # --- Preparations

    x = OrthorhombicUnitCell(1.0, 2.0, 3.0)
    y = OrthorhombicUnitCell(1.5, 2.5, 3.5)

    # --- Exercise functionality and check results

    # x ≈ (x + delta)
    @test x ≈ OrthorhombicUnitCell(1.0 + 1e-9, 2.0, 3.0)
    @test x ≈ OrthorhombicUnitCell(1.0, 2.0 + 1e-9, 3.0)
    @test x ≈ OrthorhombicUnitCell(1.0, 2.0, 3.0 - 1e-9)

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

@testset "-(::OrthorhombicUnitCell)" begin
    # --- Tests

    lattice_constants_x = (a=1, b=3, c=5)
    x = OrthorhombicUnitCell(lattice_constants_x)

    lattice_constants_y = (a=2, b=10, c=20)
    y = OrthorhombicUnitCell(lattice_constants_y)

    @test x - y == OrthorhombicUnitCellDelta(
        lattice_constants_x.a - lattice_constants_y.a,
        lattice_constants_x.b - lattice_constants_y.b,
        lattice_constants_x.c - lattice_constants_y.c,
    )
end
