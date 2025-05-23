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
rhombohedral_tests.jl defines tests for lattice/rhombohedral.jl
"""
# --- Imports

# Standard library
using Test
using LinearAlgebra: det, norm, cross

# Xtallography package
using Xtallography

# --- Tests

# ------ Types

@testset "RhombohedralLatticeConstants constructor: valid arguments" begin
    # --- Preparations

    a = 1
    α = π / 4

    # --- Tests

    # ------ basic tests

    lattice_constants = RhombohedralLatticeConstants(a, α)

    @test lattice_constants.a == a
    @test lattice_constants.α == α
end

@testset "RhombohedralLatticeConstants constructor: invalid arguments" begin
    # --- Preparations

    # Valid arguments
    a = 1
    α = π / 4

    # --- Tests

    # ------ a

    # a = 0
    expected_message = "`a` must be positive"
    @test_throws DomainError(0, expected_message) RhombohedralLatticeConstants(0, α)

    # a < 0
    expected_message = "`a` must be positive"
    @test_throws DomainError(-2.0, expected_message) RhombohedralLatticeConstants(-2.0, α)

    # ------ α

    # α = 0
    expected_message = "`α` must satisfy 0 < α < 2π / 3"
    @test_throws DomainError(0, expected_message) RhombohedralLatticeConstants(a, 0)

    # α < 0
    expected_message = "`α` must satisfy 0 < α < 2π / 3"
    @test_throws DomainError(-1.0, expected_message) RhombohedralLatticeConstants(a, -1)

    # α = 2π / 3
    expected_message = "`α` must satisfy 0 < α < 2π / 3"
    @test_throws DomainError(2π / 3, expected_message) RhombohedralLatticeConstants(
        a, 2π / 3
    )

    # α > 2π / 3
    expected_message = "`α` must satisfy 0 < α < 2π / 3"
    @test_throws DomainError(3π / 2, expected_message) RhombohedralLatticeConstants(
        a, 3π / 2
    )
end

@testset "RhombohedralLatticeConstantDeltas constructor" begin
    # --- Tests

    Δa = 1
    Δα = π / 3
    Δlattice_constants = RhombohedralLatticeConstantDeltas(Δa, Δα)

    @test Δlattice_constants.Δa == Δa
    @test Δlattice_constants.Δα == Δα
end

# ------ LatticeConstants functions

@testset "isapprox(::RhombohedralLatticeConstants)" begin
    # --- Preparations

    x = RhombohedralLatticeConstants(1.0, π / 4)
    y = RhombohedralLatticeConstants(1.5, π / 5)

    # --- Exercise functionality and check results

    # x ≈ (x + delta)
    @test x ≈ RhombohedralLatticeConstants(1.0 + 1e-9, π / 4)
    @test x ≈ RhombohedralLatticeConstants(1.0, π / 4 + 1e-9)

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

@testset "convert(::RhombohedralLatticeConstants)" begin
    # --- Tests

    x = RhombohedralLatticeConstants((rand(2) .+ 1)...)

    x_vector = convert(Vector, x)
    @test x_vector == [x.a, x.α]
end

@testset "-(::RhombohedralLatticeConstants)" begin
    # --- Tests

    x = RhombohedralLatticeConstants(1, π / 6)
    y = RhombohedralLatticeConstants(2, π / 5)
    @test x - y == RhombohedralLatticeConstantDeltas(x.a - y.a, x.α - y.α)
end

@testset "lattice_system(::RhombohedralLatticeConstants)" begin
    # --- Tests

    lattice_constants = RhombohedralLatticeConstants(1, π / 5)
    @test lattice_system(lattice_constants) === rhombohedral
end

@testset "standardize(): rhombohedral" begin
    # --- Tests

    # ------ Rhombohedral lattices have no lattice constants conventions for primitive
    #        centering

    lattice_constants = RhombohedralLatticeConstants(1.0, 2π / 5)

    standardized_lattice_constants, standardized_centering = standardize(
        lattice_constants, primitive
    )

    expected_lattice_constants = lattice_constants
    @test standardized_lattice_constants ≈ expected_lattice_constants

    # ------ Invalid centering

    for centering in (body_centered, face_centered, base_centered)
        expected_message =
            "Invalid Bravais lattice: " *
            "(lattice_system=Rhombohedral, centering=$(nameof(typeof(centering))))"

        @test_throws ArgumentError(expected_message) standardize(
            lattice_constants, centering
        )
    end
end

# ------ LatticeConstantDeltas functions

@testset "isapprox(::RhombohedralLatticeConstantDeltas)" begin
    # --- Preparations

    x = RhombohedralLatticeConstantDeltas(1.0, π / 4)
    y = RhombohedralLatticeConstantDeltas(1.5, π / 5)

    # --- Exercise functionality and check results

    # x ≈ (x + delta)
    @test x ≈ RhombohedralLatticeConstantDeltas(1.0 + 1e-9, π / 4)
    @test x ≈ RhombohedralLatticeConstantDeltas(1.0, π / 4 + 1e-9)

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

@testset "convert(::RhombohedralLatticeConstantDeltas)" begin
    # --- Tests

    x = RhombohedralLatticeConstantDeltas(rand(2)...)

    x_vector = convert(Vector, x)
    @test x_vector == [x.Δa, x.Δα]
end

@testset "lattice_system(::RhombohedralLatticeConstantDeltas)" begin
    # --- Tests

    Δlattice_constants = RhombohedralLatticeConstantDeltas(1, π / 5)
    @test lattice_system(Δlattice_constants) === rhombohedral
end

# ------ Unit cell computations

@testset "basis(::RhombohedralLatticeConstants)" begin
    # --- Preparations

    a = 2
    α = 3π / 5
    lattice_constants = RhombohedralLatticeConstants(a, α)

    # --- Exercise functionality

    basis_a, basis_b, basis_c = basis(lattice_constants)

    # --- Check results

    @test basis_a ≈ [a, 0, 0]
    @test basis_b ≈ [a * cos(α), a * sin(α), 0]
    @test basis_c ≈ [
        a * cos(α),
        a / sin(α) * (cos(α) - cos(α)^2),
        a / sin(α) * sqrt(1 - 3 * cos(α)^2 + 2 * cos(α)^3),
    ]
end

@testset "volume(::RhombohedralLatticeConstants)" begin
    # --- Preparations

    # Construct basis vectors for unit cell
    a = 6
    α = 3π / 5
    lattice_constants = RhombohedralLatticeConstants(a, α)

    basis_a, basis_b, basis_c = basis(lattice_constants)

    # --- Exercise functionality and check results

    @test volume(lattice_constants) ≈ abs(det(hcat(basis_a, basis_b, basis_c)))
end

@testset "surface_area(::RhombohedralLatticeConstants)" begin
    # --- Preparations

    # Construct basis vectors for unit cell
    a = 6
    α = 3π / 5
    lattice_constants = RhombohedralLatticeConstants(a, α)

    basis_a, basis_b, basis_c = basis(lattice_constants)

    # --- Exercise functionality and check results

    @test surface_area(lattice_constants) ≈
        2 * norm(cross(basis_a, basis_b)) +
          2 * norm(cross(basis_b, basis_c)) +
          2 * norm(cross(basis_c, basis_a))
end

@testset "reduced_cell(): rhombohedral" begin
    # --- Preparations

    a = 2
    α = 2π / 5
    lattice_constants = RhombohedralLatticeConstants(a, α)
    basis_a, basis_b, basis_c = basis(lattice_constants)

    # --- Exercise functionality and check results

    # primitive unit cell
    unit_cell = UnitCell(lattice_constants, primitive)

    expected_reduced_cell = reduced_cell(
        UnitCell(LatticeConstants(basis_a, basis_b, basis_c), primitive)
    )

    reduced_cell_ = reduced_cell(unit_cell)
    @test reduced_cell_.lattice_constants isa RhombohedralLatticeConstants
    @test volume(reduced_cell_) ≈ volume(unit_cell)
    @test reduced_cell_ ≈ expected_reduced_cell
end

@testset "is_equivalent_unit_cell(::UnitCell): rhombohedral" begin
    # --- Preparations

    a = 2
    α = 3π / 5
    lattice_constants = RhombohedralLatticeConstants(a, α)
    basis_a, basis_b, basis_c = basis(lattice_constants)

    # --- Exercise functionality and check results

    # equivalent rhombohedral and triclinic unit cells
    rhombohedral_unit_cell = UnitCell(lattice_constants, primitive)
    triclinic_unit_cell = UnitCell(
        LatticeConstants(basis_a, basis_b, basis_c; identify_lattice_system=false),
        primitive,
    )
    @test is_equivalent_unit_cell(rhombohedral_unit_cell, triclinic_unit_cell)
end

@testset "is_equivalent_unit_cell(::RhombohedralLatticeConstants)" begin
    # --- Preparations

    a_ref = 2
    α_ref = π / 4
    lattice_constants_ref = RhombohedralLatticeConstants(a_ref, α_ref)

    # --- Exercise functionality and check results

    # unit cells are equivalent
    lattice_constants_test = RhombohedralLatticeConstants(a_ref + 1e-9, α_ref - 1e-9)
    @test is_equivalent_unit_cell(lattice_constants_test, lattice_constants_ref)

    # test unit cell and reference unit cell are unrelated
    lattice_constants_test = RhombohedralLatticeConstants(2 * a_ref, α_ref)
    @test !is_equivalent_unit_cell(lattice_constants_test, lattice_constants_ref)

    # test unit cell and reference unit cell are for different lattice systems
    lattice_constants_test = CubicLatticeConstants(1)
    @test !is_equivalent_unit_cell(lattice_constants_test, lattice_constants_ref)
end
