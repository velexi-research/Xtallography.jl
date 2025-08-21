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
lattices_tests.jl defines tests for lattices.jl
"""
# --- Imports

# Standard library
using LinearAlgebra: qr, I
using Test

# External packages
using Combinatorics: permutations

# Xtallography package
using Xtallography

# --- Tests

# ------ Types

@testset "LatticeConstants Subtypes" begin
    expected_types = [
        TriclinicLatticeConstants,
        MonoclinicLatticeConstants,
        OrthorhombicLatticeConstants,
        TetragonalLatticeConstants,
        RhombohedralLatticeConstants,
        HexagonalLatticeConstants,
        CubicLatticeConstants,
    ]
    for type in expected_types
        @test type <: LatticeConstants
    end
end

# Helper function @testset "LatticeConstants(::Vector, ::Vector, ::Vector) constructor"
function test_basis_rotations_and_permutations(
    rotations::Vector,
    expected_lattice_constants::LatticeConstants,
    basis_a::Vector{<:Real},
    basis_b::Vector{<:Real},
    basis_c::Vector{<:Real};
    centering=nothing,
)
    # --- Preparations

    # Get type of expected lattice constants
    expected_type = typeof(expected_lattice_constants)

    # --- Test LatticeConstants(::Vector, ::Vector, ::Vector) for rotations of all basis
    #     permutations

    for rotation in rotations
        # Generate bases to check
        if centering == base_centering
            # Do not permute basis vectors for base-centering
            bases_to_test = ([basis_a, basis_b, basis_c],)
        else
            bases_to_test = permutations([basis_a, basis_b, basis_c])
        end

        # Exercise functionality and check results
        for basis in bases_to_test
            local lattice_constants

            # Construct LatticeConstants object
            if isnothing(centering)
                lattice_constants = LatticeConstants(basis[1], basis[2], basis[3])
            else
                lattice_constants = LatticeConstants(
                    basis[1], basis[2], basis[3]; centering=centering
                )
            end

            # Check results
            @test lattice_constants isa expected_type
            @test lattice_constants ≈ expected_lattice_constants
        end
    end
end

@testset "LatticeConstants(::Vector, ::Vector, ::Vector) constructor: cubic" begin
    # --- Preparations

    # Generate random rotation matrix
    rotations = [I, qr(rand(3)).Q, qr(rand(3)).Q]

    # Generate basis
    a = 5
    basis_a = [a, 0, 0]
    basis_b = [0, a, 0]
    basis_c = [0, 0, a]

    # --- Tests

    expected_lattice_constants = standardize(CubicLatticeConstants(a))

    # default centering keyword argument
    test_basis_rotations_and_permutations(
        rotations, expected_lattice_constants, basis_a, basis_b, basis_c
    )

    # centering keyword argument provided
    for centering in (primitive_centering, body_centering, face_centering)
        test_basis_rotations_and_permutations(
            rotations,
            expected_lattice_constants,
            basis_a,
            basis_b,
            basis_c;
            centering=centering,
        )
    end
end

@testset "LatticeConstants(::Vector, ::Vector, ::Vector) constructor: tetragonal" begin
    # --- Preparations

    # Generate random rotation matrix
    rotations = [I, qr(rand(3)).Q, qr(rand(3)).Q]

    # --- Tests

    # ------ a < c

    a = 5
    c = 10

    basis_a = [a, 0, 0]
    basis_b = [0, a, 0]
    basis_c = [0, 0, c]

    expected_lattice_constants = standardize(TetragonalLatticeConstants(a, c))

    # default centering keyword argument
    test_basis_rotations_and_permutations(
        rotations, expected_lattice_constants, basis_a, basis_b, basis_c
    )

    # centering keyword argument provided
    for centering in (primitive_centering, body_centering)
        test_basis_rotations_and_permutations(
            rotations,
            expected_lattice_constants,
            basis_a,
            basis_b,
            basis_c;
            centering=centering,
        )
    end

    # ------ a > c

    a = 10
    c = 5

    basis_a = [a, 0, 0]
    basis_b = [0, a, 0]
    basis_c = [0, 0, c]

    expected_lattice_constants = standardize(TetragonalLatticeConstants(a, c))

    # default centering keyword argument
    test_basis_rotations_and_permutations(
        rotations, expected_lattice_constants, basis_a, basis_b, basis_c
    )

    # centering keyword argument provided
    for centering in (primitive_centering, body_centering)
        test_basis_rotations_and_permutations(
            rotations,
            expected_lattice_constants,
            basis_a,
            basis_b,
            basis_c;
            centering=centering,
        )
    end
end

@testset "LatticeConstants(::Vector, ::Vector, ::Vector) constructor: orthorhombic" begin
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
    expected_lattice_constants = OrthorhombicLatticeConstants(a, b, c)
    test_basis_rotations_and_permutations(
        rotations, expected_lattice_constants, basis_a, basis_b, basis_c
    )

    # centering = primitive-centering, body-centering, or face-centering
    expected_lattice_constants = OrthorhombicLatticeConstants(a, b, c)
    for centering in (primitive_centering, body_centering, face_centering)
        test_basis_rotations_and_permutations(
            rotations,
            expected_lattice_constants,
            basis_a,
            basis_b,
            basis_c;
            centering=centering,
        )
    end

    # centering = base-centering; a, b < c
    a = 5
    b = 8
    c = 10

    basis_a = [a, 0, 0]
    basis_b = [0, b, 0]
    basis_c = [0, 0, c]

    expected_lattice_constants = OrthorhombicLatticeConstants(a, b, c)

    test_basis_rotations_and_permutations(
        rotations,
        expected_lattice_constants,
        basis_a,
        basis_b,
        basis_c;
        centering=base_centering,
    )

    # centering = base-centering; a < c < b
    a = 5
    b = 10
    c = 8

    basis_a = [a, 0, 0]
    basis_b = [0, b, 0]
    basis_c = [0, 0, c]

    expected_lattice_constants = OrthorhombicLatticeConstants(a, b, c)
    test_basis_rotations_and_permutations(
        rotations,
        expected_lattice_constants,
        basis_a,
        basis_b,
        basis_c;
        centering=base_centering,
    )

    # centering = base-centering; b < c < a
    a = 10
    b = 5
    c = 8

    basis_a = [a, 0, 0]
    basis_b = [0, b, 0]
    basis_c = [0, 0, c]

    expected_lattice_constants = OrthorhombicLatticeConstants(b, a, c)
    test_basis_rotations_and_permutations(
        rotations,
        expected_lattice_constants,
        basis_a,
        basis_b,
        basis_c;
        centering=base_centering,
    )
end

@testset "LatticeConstants(::Vector, ::Vector, ::Vector) constructor: rhombohedral" begin
    # --- Preparations

    # Generate random rotation matrix
    rotations = [I, qr(rand(3)).Q, qr(rand(3)).Q]

    # --- Tests

    a = 5
    α = 2π / 5

    basis_a = [a, 0, 0]
    basis_b = [a * cos(α), a * sin(α), 0]
    basis_c = [
        a * cos(α),
        a / sin(α) * (cos(α) - cos(α)^2),
        a / sin(α) * sqrt(1 - 3 * cos(α)^2 + 2 * cos(α)^3),
    ]

    expected_lattice_constants = standardize(RhombohedralLatticeConstants(a, α))

    # default centering keyword argument
    test_basis_rotations_and_permutations(
        rotations, expected_lattice_constants, basis_a, basis_b, basis_c
    )

    # centering keyword argument provided
    for centering in (primitive_centering,)
        test_basis_rotations_and_permutations(
            rotations,
            expected_lattice_constants,
            basis_a,
            basis_b,
            basis_c;
            centering=centering,
        )
    end
end

@testset "LatticeConstants(::Vector, ::Vector, ::Vector) constructor: hexagonal" begin
    # --- Preparations

    # Generate random rotation matrix
    rotations = [I, qr(rand(3)).Q, qr(rand(3)).Q]

    # --- Tests

    # ------ a < c, β = 2π / 3

    a = 3
    c = 5

    basis_a = [a, 0, 0]
    basis_b = [-0.5 * a, 0.5 * a * sqrt(3), 0]
    basis_c = [0, 0, c]

    expected_lattice_constants = standardize(HexagonalLatticeConstants(a, c))

    # default centering keyword argument
    test_basis_rotations_and_permutations(
        rotations, expected_lattice_constants, basis_a, basis_b, basis_c
    )

    # centering keyword argument provided
    for centering in (primitive_centering,)
        test_basis_rotations_and_permutations(
            rotations,
            expected_lattice_constants,
            basis_a,
            basis_b,
            basis_c;
            centering=centering,
        )
    end

    # ------ a < c, β = π / 3

    a = 3
    c = 5

    basis_a = [a, 0, 0]
    basis_b = [0.5 * a, 0.5 * a * sqrt(3), 0]
    basis_c = [0, 0, c]

    expected_lattice_constants = standardize(HexagonalLatticeConstants(a, c))

    # default centering keyword argument
    test_basis_rotations_and_permutations(
        rotations, expected_lattice_constants, basis_a, basis_b, basis_c
    )

    # centering keyword argument provided
    for centering in (primitive_centering,)
        test_basis_rotations_and_permutations(
            rotations,
            expected_lattice_constants,
            basis_a,
            basis_b,
            basis_c;
            centering=centering,
        )
    end

    # ------ a > c, β = 2π / 3

    a = 5
    c = 3

    basis_a = [a, 0, 0]
    basis_b = [-0.5 * a, 0.5 * a * sqrt(3), 0]
    basis_c = [0, 0, c]

    expected_lattice_constants = standardize(HexagonalLatticeConstants(a, c))

    # default centering keyword argument
    test_basis_rotations_and_permutations(
        rotations, expected_lattice_constants, basis_a, basis_b, basis_c
    )

    # centering keyword argument provided
    for centering in (primitive_centering,)
        test_basis_rotations_and_permutations(
            rotations,
            expected_lattice_constants,
            basis_a,
            basis_b,
            basis_c;
            centering=centering,
        )
    end

    # ------ a > c, β = π / 3

    a = 5
    c = 3

    basis_a = [a, 0, 0]
    basis_b = [0.5 * a, 0.5 * a * sqrt(3), 0]
    basis_c = [0, 0, c]

    expected_lattice_constants = standardize(HexagonalLatticeConstants(a, c))

    # default centering keyword argument
    test_basis_rotations_and_permutations(
        rotations, expected_lattice_constants, basis_a, basis_b, basis_c
    )

    # centering keyword argument provided
    for centering in (primitive_centering,)
        test_basis_rotations_and_permutations(
            rotations,
            expected_lattice_constants,
            basis_a,
            basis_b,
            basis_c;
            centering=centering,
        )
    end
end

@testset "LatticeConstants(::Vector, ::Vector, ::Vector) constructor: monoclinic" begin
    # --- Preparations

    # Generate random rotation matrix
    rotations = [I, qr(rand(3)).Q, qr(rand(3)).Q]

    # --- Tests

    a = 5
    b = 8
    c = 10
    β = 3π / 5

    basis_a = [a, 0, 0]
    basis_b = [0, b, 0]
    basis_c = [c * cos(β), 0, c * sin(β)]

    # default centering keyword argument
    expected_lattice_constants = standardize(MonoclinicLatticeConstants(a, b, c, β))
    test_basis_rotations_and_permutations(
        rotations, expected_lattice_constants, basis_a, basis_b, basis_c
    )

    # centering = primitive
    expected_lattice_constants = standardize(MonoclinicLatticeConstants(a, b, c, β))
    test_basis_rotations_and_permutations(
        rotations,
        expected_lattice_constants,
        basis_a,
        basis_b,
        basis_c;
        centering=primitive_centering,
    )

    # centering = body-centering
    expected_lattice_constants, _ = standardize(
        MonoclinicLatticeConstants(a, b, c, β), body_centering
    )
    test_basis_rotations_and_permutations(
        rotations,
        expected_lattice_constants,
        basis_a,
        basis_b,
        basis_c;
        centering=body_centering,
    )

    # centering = base-centering
    expected_lattice_constants, _ = standardize(
        MonoclinicLatticeConstants(a, b, c, β), base_centering
    )
    test_basis_rotations_and_permutations(
        rotations,
        expected_lattice_constants,
        basis_a,
        basis_b,
        basis_c;
        centering=base_centering,
    )
end

@testset "LatticeConstants(::Vector, ::Vector, ::Vector) constructor: triclinic" begin
    # --- Preparations

    # Generate random rotation matrix
    rotations = [I, qr(rand(3)).Q, qr(rand(3)).Q]

    # --- Tests

    # ------ General triclinic unit cell

    a = 5
    b = 8
    c = 10
    α = 2π / 5
    β = 3π / 5
    γ = 4π / 5

    basis_a = [a, 0, 0]
    basis_b = [b * cos(γ), b * sin(γ), 0]
    basis_c = [
        c * cos(β),
        c / sin(γ) * (cos(α) - cos(β) * cos(γ)),
        2 * c / sin(γ) * sqrt(
            sin(0.5 * (α + β + γ)) *
            sin(0.5 * (α - β + γ)) *
            sin(0.5 * (α + β - γ)) *
            sin(0.5 * (-α + β + γ)),
        ),
    ]

    expected_lattice_constants = standardize(TriclinicLatticeConstants(a, b, c, α, β, γ))

    # default centering keyword argument
    test_basis_rotations_and_permutations(
        rotations, expected_lattice_constants, basis_a, basis_b, basis_c
    )

    # centering keyword argument provided
    for centering in (primitive_centering,)
        test_basis_rotations_and_permutations(
            rotations,
            expected_lattice_constants,
            basis_a,
            basis_b,
            basis_c;
            centering=centering,
        )
    end

    # ------ Almost rhomobohedral: all three angles equal, edge lengths unequal

    a = 5
    b = 6
    c = 7
    α = 2π / 5

    basis_a = [a, 0, 0]
    basis_b = [b * cos(α), b * sin(α), 0]
    basis_c = [
        c * cos(α),
        c / sin(α) * (cos(α) - cos(α)^2),
        c / sin(α) * sqrt(1 - 3 * cos(α)^2 + 2 * cos(α)^3),
    ]

    expected_lattice_constants = standardize(TriclinicLatticeConstants(a, b, c, α, α, α))

    # default centering keyword argument
    test_basis_rotations_and_permutations(
        rotations, expected_lattice_constants, basis_a, basis_b, basis_c
    )

    # centering keyword argument provided
    for centering in (primitive_centering,)
        test_basis_rotations_and_permutations(
            rotations,
            expected_lattice_constants,
            basis_a,
            basis_b,
            basis_c;
            centering=centering,
        )
    end
end

@testset "LatticeConstants(::Vector, ::Vector, ::Vector) constructor: identify_lattice_system = false" begin
    # --- Preparations

    a = 5
    basis_a = [a, 0, 0]
    basis_b = [0, a, 0]
    basis_c = [0, 0, a]

    # --- Exercise functionality and check results

    expected_lattice_constants = TriclinicLatticeConstants(a, a, a, π / 2, π / 2, π / 2)

    lattice_constants = LatticeConstants(
        basis_a, basis_b, basis_c; identify_lattice_system=false
    )
    @test lattice_constants isa TriclinicLatticeConstants
    @test lattice_constants ≈ expected_lattice_constants
end

@testset "UnitCell constructor" begin
    # --- Exercise functionality and check results

    # ------ Cubic unit cells

    lattice_constants = CubicLatticeConstants(1)

    for centering in (primitive_centering, body_centering, face_centering)
        unit_cell = UnitCell(lattice_constants, centering)

        @test unit_cell.lattice_constants == lattice_constants
        @test unit_cell.centering == centering
    end

    # ------ tetragonal unit cells

    lattice_constants = TetragonalLatticeConstants(1, 2)

    for centering in (primitive_centering, body_centering)
        unit_cell = UnitCell(lattice_constants, centering)

        @test unit_cell.lattice_constants == lattice_constants
        @test unit_cell.centering == centering
    end

    # ------ orthorhombic unit cells

    lattice_constants = OrthorhombicLatticeConstants(1, 2, 3)

    for centering in (primitive_centering, body_centering, face_centering, base_centering)
        unit_cell = UnitCell(lattice_constants, centering)

        @test unit_cell.lattice_constants == lattice_constants
        @test unit_cell.centering == centering
    end

    # ------ hexagonal unit cells

    lattice_constants = HexagonalLatticeConstants(1, 2)

    for centering in (primitive_centering,)
        unit_cell = UnitCell(lattice_constants, centering)

        @test unit_cell.lattice_constants == lattice_constants
        @test unit_cell.centering == centering
    end

    # ------ rhombohedral unit cells

    lattice_constants = RhombohedralLatticeConstants(1, π / 3)

    for centering in (primitive_centering,)
        unit_cell = UnitCell(lattice_constants, centering)

        @test unit_cell.lattice_constants == lattice_constants
        @test unit_cell.centering == centering
    end

    # ------ monoclinic unit cells

    lattice_constants = MonoclinicLatticeConstants(1, 2, 3, 3π / 5)

    for centering in (primitive_centering, body_centering, base_centering)
        unit_cell = UnitCell(lattice_constants, centering)

        @test unit_cell.lattice_constants == lattice_constants
        @test unit_cell.centering == centering
    end

    # ------ triclinic unit cells

    lattice_constants = TriclinicLatticeConstants(1, 2, 3, 2π / 5, 3π / 5, 4π / 5)

    for centering in (primitive_centering,)
        unit_cell = UnitCell(lattice_constants, centering)

        @test unit_cell.lattice_constants == lattice_constants
        @test unit_cell.centering == centering
    end
end

# ------ lattice methods

@testset "lattice_system(::UnitCell)" begin
    # --- Tests

    @test lattice_system(UnitCell(CubicLatticeConstants(1), primitive_centering)) == cubic
    @test lattice_system(UnitCell(TetragonalLatticeConstants(1, 2), base_centering)) ==
        tetragonal
end

@testset "is_bravais_lattice(::UnitCell)" begin
    # --- Tests

    # Cubic
    @test is_bravais_lattice(UnitCell(CubicLatticeConstants(1), primitive_centering))
    @test is_bravais_lattice(UnitCell(CubicLatticeConstants(1), body_centering))
    @test is_bravais_lattice(UnitCell(CubicLatticeConstants(1), face_centering))
    @test !is_bravais_lattice(UnitCell(CubicLatticeConstants(1), base_centering))

    # Tetragonal
    @test is_bravais_lattice(
        UnitCell(TetragonalLatticeConstants(1, 2), primitive_centering)
    )
    @test is_bravais_lattice(UnitCell(TetragonalLatticeConstants(1, 2), body_centering))
    @test !is_bravais_lattice(UnitCell(TetragonalLatticeConstants(1, 2), face_centering))
    @test !is_bravais_lattice(UnitCell(TetragonalLatticeConstants(1, 2), base_centering))

    # Orthorhombic
    @test is_bravais_lattice(
        UnitCell(OrthorhombicLatticeConstants(1, 2, 3), primitive_centering)
    )
    @test is_bravais_lattice(
        UnitCell(OrthorhombicLatticeConstants(1, 2, 3), body_centering)
    )
    @test is_bravais_lattice(
        UnitCell(OrthorhombicLatticeConstants(1, 2, 3), face_centering)
    )
    @test is_bravais_lattice(
        UnitCell(OrthorhombicLatticeConstants(1, 2, 3), base_centering)
    )

    # Hexagonal
    @test is_bravais_lattice(UnitCell(HexagonalLatticeConstants(1, 2), primitive_centering))
    @test !is_bravais_lattice(UnitCell(HexagonalLatticeConstants(1, 2), body_centering))
    @test !is_bravais_lattice(UnitCell(HexagonalLatticeConstants(1, 2), face_centering))
    @test !is_bravais_lattice(UnitCell(HexagonalLatticeConstants(1, 2), base_centering))

    # Rhombohedral
    @test is_bravais_lattice(
        UnitCell(RhombohedralLatticeConstants(1, π / 3), primitive_centering)
    )
    @test !is_bravais_lattice(
        UnitCell(RhombohedralLatticeConstants(1, π / 3), body_centering)
    )
    @test !is_bravais_lattice(
        UnitCell(RhombohedralLatticeConstants(1, π / 3), face_centering)
    )
    @test !is_bravais_lattice(
        UnitCell(RhombohedralLatticeConstants(1, π / 3), base_centering)
    )

    # Monoclinic
    @test is_bravais_lattice(
        UnitCell(MonoclinicLatticeConstants(1, 2, 3, 3π / 5), primitive_centering)
    )
    @test is_bravais_lattice(
        UnitCell(MonoclinicLatticeConstants(1, 2, 3, 3π / 5), body_centering)
    )
    @test !is_bravais_lattice(
        UnitCell(MonoclinicLatticeConstants(1, 2, 3, 3π / 5), face_centering)
    )
    @test is_bravais_lattice(
        UnitCell(MonoclinicLatticeConstants(1, 2, 3, 3π / 5), base_centering)
    )

    # Triclinic
    @test is_bravais_lattice(
        UnitCell(
            TriclinicLatticeConstants(1, 2, 3, 2π / 5, 3π / 5, 4π / 5), primitive_centering
        ),
    )
    @test !is_bravais_lattice(
        UnitCell(TriclinicLatticeConstants(1, 2, 3, 2π / 5, 3π / 5, 4π / 5), body_centering)
    )
    @test !is_bravais_lattice(
        UnitCell(TriclinicLatticeConstants(1, 2, 3, 2π / 5, 3π / 5, 4π / 5), face_centering)
    )
    @test !is_bravais_lattice(
        UnitCell(TriclinicLatticeConstants(1, 2, 3, 2π / 5, 3π / 5, 4π / 5), base_centering)
    )
end

# ------ LatticeConstants methods

@testset "isapprox(::LatticeConstants): comparison between different types" begin
    # --- Tests

    # x ≈ y
    x = TetragonalLatticeConstants(1 + 1e-9, 2 - 1e-9)
    y = TetragonalLatticeConstants(1, 2)
    @test x ≈ y

    # x ≉ y
    x = CubicLatticeConstants(1)
    y = TetragonalLatticeConstants(1, 2)
    @test x ≉ y
end

@testset "standardize(::UnitCell)" begin
    # --- Tests

    # Check sequence of method calls works
    unit_cell = UnitCell(
        TriclinicLatticeConstants(1, 2, 3, 2π / 5, 3π / 5, 4π / 5), primitive_centering
    )
    standardized_unit_cell = standardize(unit_cell)
    expected_standardized_unit_cell = UnitCell(
        TriclinicLatticeConstants(1, 2, 3, 2π / 5, 2π / 5, π / 5), primitive_centering
    )
    @test standardized_unit_cell ≈ expected_standardized_unit_cell
end

@testset "standardize(::LatticeConstants)" begin
    # --- Tests

    # Check that no centering is returned
    lattice_constants = TriclinicLatticeConstants(1, 2, 3, 2π / 5, 3π / 5, 4π / 5)
    expected_standardized_lattice_constants, _ = standardize(
        TriclinicLatticeConstants(1, 2, 3, 2π / 5, 3π / 5, 4π / 5), primitive_centering
    )

    @test standardize(TriclinicLatticeConstants(1, 2, 3, 2π / 5, 3π / 5, 4π / 5)) ≈
        expected_standardized_lattice_constants
end

# ------ LatticeConstantDeltas methods

@testset "isapprox(::LatticeConstantDeltas): comparison between different types" begin
    # --- Tests

    # x ≈ y
    x = CubicLatticeConstantDeltas(1 + 1e-9)
    y = CubicLatticeConstantDeltas(1)
    @test x ≈ y

    # x ≉ y
    x = CubicLatticeConstantDeltas(1)
    y = TetragonalLatticeConstantDeltas(1, 2)
    @test x ≉ y
end

# ------ UnitCell methods

@testset "isapprox(::UnitCell)" begin
    # --- Preparations

    x = UnitCell(CubicLatticeConstants(1), primitive_centering)

    # --- Exercise functionality and check results

    # x.centering != y.centering
    y = UnitCell(CubicLatticeConstants(1), body_centering)
    @test x ≉ y

    # x.lattice_constants ≈ (x.lattice_constants + delta)
    y = UnitCell(CubicLatticeConstants(1 + 1e-8), primitive_centering)
    @test x ≈ y

    # x.lattice_constants ≉ y.lattice_constants
    y = UnitCell(CubicLatticeConstants(2), primitive_centering)
    @test x ≉ y
end

@testset "basis(::UnitCell)" begin
    # --- Tests

    a = 2
    c = 5
    unit_cell = UnitCell(HexagonalLatticeConstants(a, c), primitive_centering)
    basis_a, basis_b, basis_c = basis(unit_cell)

    # Check results
    @test basis_a ≈ [a, 0, 0]
    @test basis_b ≈ [-0.5 * a, a * sqrt(3) / 2, 0]
    @test basis_c ≈ [0, 0, c]
end

@testset "volume(::UnitCell)" begin
    # --- Tests

    unit_cell = UnitCell(OrthorhombicLatticeConstants(1, 2, 3), primitive_centering)
    @test volume(unit_cell) ≈ 6
end

@testset "surface_area(::UnitCell)" begin
    # --- Tests

    unit_cell = UnitCell(OrthorhombicLatticeConstants(1, 2, 3), primitive_centering)
    @test surface_area(unit_cell) ≈ 22
end

@testset "reduced_cell(): minimum sum of length squared not unique" begin
    # --- Preparations

    lattice_constants = TriclinicLatticeConstants(
        sqrt(6), sqrt(8), sqrt(8), π / 3, acos(sqrt(3) / 6), acos(sqrt(3) / 4)
    )
    centering = primitive_centering
    unit_cell = UnitCell(lattice_constants, centering)

    # --- Exercise functionality

    reduced_cell_ = reduced_cell(unit_cell)

    # --- Check results

    # Check lattice constants
    @test isapprox(reduced_cell_.lattice_constants.a, 2.449; atol=0.0005)
    @test isapprox(reduced_cell_.lattice_constants.b, 2.828; atol=0.0005)
    @test isapprox(reduced_cell_.lattice_constants.c, 2.828; atol=0.0005)
    @test isapprox(reduced_cell_.lattice_constants.α, 104.47 * π / 180; atol=0.0005)
    @test isapprox(reduced_cell_.lattice_constants.β, 106.78 * π / 180; atol=0.0005)
    @test isapprox(reduced_cell_.lattice_constants.γ, 115.66 * π / 180; atol=0.0005)
    @test reduced_cell_.centering == primitive_centering
end

@testset "is_equivalent_unit_cell(::UnitCell): valid arguments" begin
    # --- Tests

    # equivalent unit cells, default atol and rtol
    b = 2
    c = 3
    β = 3π / 5
    a = -2 * c * cos(β)
    unit_cell_ref = UnitCell(MonoclinicLatticeConstants(a, b, c, β), primitive_centering)
    unit_cell_test = UnitCell(
        OrthorhombicLatticeConstants(a, 2 * c * sin(β), b), base_centering
    )

    @test is_equivalent_unit_cell(unit_cell_test, unit_cell_ref)

    # nonequivalent unit cells that are equivalent when atol is sufficiently large
    a = 1
    b = 2
    c = 3
    unit_cell_ref = UnitCell(OrthorhombicLatticeConstants(a, b, c), primitive_centering)
    unit_cell_test = UnitCell(
        OrthorhombicLatticeConstants(a + 2, b - 1, c + 5), primitive_centering
    )
    @test !is_equivalent_unit_cell(unit_cell_test, unit_cell_ref)
    @test !is_equivalent_unit_cell(unit_cell_test, unit_cell_ref; atol=4)
    @test is_equivalent_unit_cell(unit_cell_test, unit_cell_ref; atol=6)

    # nonequivalent unit cells that are equivalent when rtol is sufficiently large
    a = 1
    b = 2
    c = 3
    unit_cell_ref = UnitCell(OrthorhombicLatticeConstants(a, b, c), primitive_centering)
    unit_cell_test = UnitCell(
        OrthorhombicLatticeConstants(a + 2, b - 1, c + 5), primitive_centering
    )
    @test !is_equivalent_unit_cell(unit_cell_test, unit_cell_ref)
    @test !is_equivalent_unit_cell(unit_cell_test, unit_cell_ref; atol=0, rtol=0.5)
    @test is_equivalent_unit_cell(unit_cell_test, unit_cell_ref; atol=0, rtol=1)

    # unit cells for different lattice systems
    b = 2
    c = 3
    β = 3π / 5
    a = -2 * c * cos(β)
    unit_cell_ref = UnitCell(MonoclinicLatticeConstants(a, b, c, β), primitive_centering)
    unit_cell_test = UnitCell(CubicLatticeConstants(a), face_centering)

    @test !is_equivalent_unit_cell(unit_cell_test, unit_cell_ref)
end

@testset "is_equivalent_unit_cell(::UnitCell): invalid arguments" begin
    # --- Preparations

    unit_cell_ref = UnitCell(CubicLatticeConstants(1.0), primitive_centering)
    unit_cell_test = UnitCell(CubicLatticeConstants(1.0), primitive_centering)

    # --- Exercise functionality and check results

    # atol < 0
    expected_message = "`atol` must be nonnegative"
    @test_throws DomainError(-1, expected_message) is_equivalent_unit_cell(
        unit_cell_test, unit_cell_ref; atol=-1
    )

    # rtol < 0
    expected_message = "`rtol` must be nonnegative"
    @test_throws DomainError(-1, expected_message) is_equivalent_unit_cell(
        unit_cell_test, unit_cell_ref; rtol=-1
    )
end

@testset "is_equivalent_unit_cell(::LatticeConstants): valid arguments" begin
    # --- Tests

    # ------ equivalent unit cells

    a = 1
    b = 2
    c = 3
    β = 3π / 5
    lattice_constants_ref = MonoclinicLatticeConstants(a, b, c, β)

    c_alt = sqrt(a^2 + c^2 + 2 * a * c * cos(β))
    β_alt = π - asin(sin(β) / c_alt * c)
    lattice_constants_test = MonoclinicLatticeConstants(a, b, c_alt, β_alt)

    @test is_equivalent_unit_cell(lattice_constants_test, lattice_constants_ref)

    # ------ nonequivalent unit cells

    # Case #1
    b = 2
    c = 3
    β = 3π / 5
    a = -2 * c * cos(β)
    lattice_constants_ref = MonoclinicLatticeConstants(a, b, c, β)
    lattice_constants_test = OrthorhombicLatticeConstants(a, 2 * c * sin(β), b)

    @test !is_equivalent_unit_cell(lattice_constants_test, lattice_constants_ref)

    # Case #2
    b = 2
    c = 3
    β = 3π / 5
    a = -2 * c * cos(β)
    lattice_constants_ref = UnitCell(
        MonoclinicLatticeConstants(a, b, c, β), primitive_centering
    )
    lattice_constants_test = UnitCell(CubicLatticeConstants(a), face_centering)

    @test !is_equivalent_unit_cell(lattice_constants_test, lattice_constants_ref)
end

@testset "is_equivalent_unit_cell(::LatticeConstants): invalid arguments" begin
    # --- Preparations

    lattice_constants_ref = CubicLatticeConstants(1.0)
    lattice_constants_test = CubicLatticeConstants(1.0)

    # --- Exercise functionality and check results

    # atol < 0
    expected_message = "`atol` must be nonnegative"
    @test_throws DomainError(-1, expected_message) is_equivalent_unit_cell(
        lattice_constants_test, lattice_constants_ref; atol=-1
    )

    # rtol < 0
    expected_message = "`rtol` must be nonnegative"
    @test_throws DomainError(-1, expected_message) is_equivalent_unit_cell(
        lattice_constants_test, lattice_constants_ref; rtol=-1
    )
end

@testset "is_supercell(): different LatticeConstant types" begin
    # --- Tests

    cubic_lattice_constants = CubicLatticeConstants(1.0)
    orthorhombic_lattice_constants = OrthorhombicLatticeConstants(1.0, 2.0, 3.0)

    @test !is_supercell(cubic_lattice_constants, orthorhombic_lattice_constants)
end
