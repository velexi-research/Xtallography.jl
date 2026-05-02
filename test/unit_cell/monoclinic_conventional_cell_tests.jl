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
Tests for the unit cell standardization methods for monoclinic lattices
"""
# --- Imports

# Standard library
using LinearAlgebra: dot, norm
using Logging
using Test

# Xtallography package
using Xtallography

# Notes
# =====
# These tests adopt the following variable conventions.
#
# - Unless otherwise noted, lattice constants and basis vectors refer to the monoclinic
#   (not orthorhombic or rhombohedral) unit cell.
#
# - Lattice constants and basis vectors for orthorhombic unit cells are indicated by
#   the "o_" prefix.
#
# - Lattice constants and basis vectors for rhombohedral unit cells are indicated by
#   the "r_" prefix.

# --- Tests

@testset "conventional_cell(::MonoclinicUnitCell): invalid arguments" begin
    # --- Preparations

    # lattice constants for unit cell
    a = 1.0
    b = 2.0
    c = 3.0
    β = 3π / 5

    # --- Tests

    # ------ Invalid centering

    for centering_ in (face_centering,)
        expected_message =
            "Invalid Bravais lattice: " *
            "(lattice_system=Monoclinic, centering=$(nameof(typeof(centering_))))"

        unit_cell = MonoclinicUnitCell(a, b, c, β; centering=centering_)
        @test_throws ArgumentError(expected_message) conventional_cell(unit_cell)
    end
end

@testset "conventional_cell(::MonoclinicUnitCell): non-limiting cases" begin
    # monoclinic unit cells that are not equivalent to any higher symmetry unit cell

    # --- centering = primitive_centering

    # Construct monoclinic unit cell
    a = 1.0
    b = 2.0
    c = 3.0
    β = 4π / 7
    unit_cell = MonoclinicUnitCell(a, b, c, β)

    # Check test conditions
    @test unit_cell isa MonoclinicUnitCell
    @test centering(unit_cell) === primitive_centering

    # Exercise functionality
    iucr_unit_cell = conventional_cell(unit_cell)

    # Check results
    expected_unit_cell = standardize(unit_cell)
    @test iucr_unit_cell == expected_unit_cell

    # --- centering = body_centering

    # Construct monoclinic unit cell
    a_ref = 1.0
    b_ref = 2.0
    c_ref = 3.0
    β_ref = 4π / 7

    a = a_ref
    b = b_ref
    c = sqrt((2 * a_ref)^2 + c_ref^2 + 2 * (2 * a_ref) * c_ref * cos(β_ref))
    β = asin(sin(β_ref) / c * c_ref)
    unit_cell = MonoclinicUnitCell(a, b, c, β; centering=body_centering)

    # Check test conditions
    @test unit_cell isa MonoclinicUnitCell
    @test centering(unit_cell) === body_centering

    # Exercise functionality
    iucr_unit_cell = conventional_cell(unit_cell)

    # Check results
    expected_unit_cell = MonoclinicUnitCell(
        a_ref, b_ref, c_ref, β_ref; centering=body_centering
    )
    @test iucr_unit_cell ≈ expected_unit_cell
end

@testset "conventional_cell(::MonoclinicUnitCell): limiting cases - mP --> oP, oC" begin
    # --- Preparations

    # Construct basis for orthorhombic unit cell
    o_a = 1.0
    o_b = 2.0
    o_c = 3.0

    primitive_orthorhombic_unit_cell = OrthorhombicUnitCell(o_a, o_b, o_c)
    o_basis_a, o_basis_b, o_basis_c = basis(primitive_orthorhombic_unit_cell)

    # --- Tests

    # ------ β = π / 2

    # Construct monoclinic unit cell
    a = o_a
    b = o_b
    c = o_c
    β = π / 2
    unit_cell = MonoclinicUnitCell(a, b, c, β)

    # Exercise functionality
    iucr_unit_cell = conventional_cell(unit_cell)

    # Check results
    expected_unit_cell = standardize(
        OrthorhombicUnitCell(o_a, o_b, o_c; centering=primitive_centering)
    )
    @test iucr_unit_cell ≈ expected_unit_cell

    # ------ a = -2 * c * cos(β)

    # Construct monoclinic unit cell
    basis_a = o_basis_a
    basis_b = 0.5 * (o_basis_a + o_basis_b)
    basis_c = o_basis_c
    unit_cell = standardize(
        UnitCell(basis_a, basis_b, basis_c; centering=primitive_centering)
    )

    # Check test conditions
    @test unit_cell isa MonoclinicUnitCell
    @test centering(unit_cell) === primitive_centering
    @test lattice_constants(unit_cell).a ≈
        -2 * lattice_constants(unit_cell).c * cos(lattice_constants(unit_cell).β)

    # Exercise functionality
    iucr_unit_cell = conventional_cell(unit_cell)

    # Check results
    expected_unit_cell = standardize(
        OrthorhombicUnitCell(o_a, o_b, o_c; centering=base_centering)
    )
    @test iucr_unit_cell ≈ expected_unit_cell

    # ------ a = c, basis_b = o_basis_c

    # Construct monoclinic unit cell
    basis_a = 0.5 * (o_basis_a + o_basis_b)
    basis_b = o_basis_c
    basis_c = 0.5 * (o_basis_a - o_basis_b)
    unit_cell = MonoclinicUnitCell(
        norm(basis_a),
        norm(basis_b),
        norm(basis_c),
        π - acos(dot(basis_a, basis_c) / norm(basis_a) / norm(basis_c)),
    )

    # Check test conditions
    @test unit_cell isa MonoclinicUnitCell
    @test centering(unit_cell) === primitive_centering
    @test lattice_constants(unit_cell).a ≈ lattice_constants(unit_cell).c

    # Exercise functionality
    iucr_unit_cell = conventional_cell(unit_cell)

    # Check results
    expected_unit_cell = standardize(
        OrthorhombicUnitCell(o_a, o_b, o_c; centering=base_centering)
    )
    @test iucr_unit_cell ≈ expected_unit_cell
end

@testset "conventional_cell(::MonoclinicUnitCell): limiting cases - mI --> oI, oC, oF" begin
    # --- Preparations

    # Construct basis for orthorhombic unit cell
    o_a = 1.0
    o_b = 2.0
    o_c = 3.0

    primitive_orthorhombic_unit_cell = OrthorhombicUnitCell(o_a, o_b, o_c)
    o_basis_a, o_basis_b, o_basis_c = basis(primitive_orthorhombic_unit_cell)

    # --- Tests

    # ------ β = π / 2

    # Construct monoclinic unit cell
    a = o_a
    b = o_b
    c = o_c
    β = π / 2
    unit_cell = MonoclinicUnitCell(a, b, c, β; centering=body_centering)

    # Exercise functionality
    iucr_unit_cell = conventional_cell(unit_cell)

    # Check results
    expected_unit_cell = standardize(
        OrthorhombicUnitCell(a, b, c; centering=body_centering)
    )
    @test iucr_unit_cell ≈ expected_unit_cell

    # ------ a = -c cos(β), basis_b = o_basis_a

    # Construct monoclinic unit cell
    basis_a = o_basis_c
    basis_b = o_basis_a
    basis_c = o_basis_c - o_basis_b
    unit_cell = MonoclinicUnitCell(
        norm(basis_a),
        norm(basis_b),
        norm(basis_c),
        π - acos(dot(basis_a, basis_c) / norm(basis_a) / norm(basis_c));
        centering=body_centering,
    )

    # Check test conditions
    @test unit_cell isa MonoclinicUnitCell
    @test centering(unit_cell) === body_centering
    @test lattice_constants(unit_cell).a ≈
        -lattice_constants(unit_cell).c * cos(lattice_constants(unit_cell).β)

    # Exercise functionality
    iucr_unit_cell = conventional_cell(unit_cell)

    # Check results
    expected_unit_cell = standardize(
        OrthorhombicUnitCell(a, b, c; centering=base_centering)
    )
    @test iucr_unit_cell ≈ expected_unit_cell

    # ------ a = -c cos(β), basis_b = o_basis_b

    # Construct monoclinic unit cell
    basis_a = o_basis_c
    basis_b = o_basis_b
    basis_c = o_basis_c - o_basis_a
    unit_cell = MonoclinicUnitCell(
        norm(basis_a),
        norm(basis_b),
        norm(basis_c),
        π - acos(dot(basis_a, basis_c) / norm(basis_a) / norm(basis_c));
        centering=body_centering,
    )

    # Check test conditions
    @test unit_cell isa MonoclinicUnitCell
    @test centering(unit_cell) === body_centering
    @test lattice_constants(unit_cell).a ≈
        -lattice_constants(unit_cell).c * cos(lattice_constants(unit_cell).β)

    # Exercise functionality
    iucr_unit_cell = conventional_cell(unit_cell)

    # Check results
    expected_unit_cell = standardize(
        OrthorhombicUnitCell(a, b, c; centering=base_centering)
    )
    @test iucr_unit_cell ≈ expected_unit_cell

    # ------ a = c

    # Construct monoclinic unit cell
    basis_a = 0.5 * (o_basis_a + o_basis_c)
    basis_b = o_basis_b
    basis_c = 0.5 * (o_basis_a - o_basis_c)
    unit_cell = MonoclinicUnitCell(
        norm(basis_a),
        norm(basis_b),
        norm(basis_c),
        π - acos(dot(basis_a, basis_c) / norm(basis_a) / norm(basis_c));
        centering=body_centering,
    )

    # Check test conditions
    @test unit_cell isa MonoclinicUnitCell
    @test centering(unit_cell) === body_centering
    @test lattice_constants(unit_cell).a ≈ lattice_constants(unit_cell).c

    # Exercise functionality
    iucr_unit_cell = conventional_cell(unit_cell)

    # Check results
    expected_unit_cell = standardize(
        OrthorhombicUnitCell(a, b, c; centering=face_centering)
    )
    @test iucr_unit_cell ≈ expected_unit_cell
end

@testset "conventional_cell(::MonoclinicUnitCell): limiting cases - mI --> hR" begin
    # --- Tests

    # ------ r_α < π/3: a^2 + b^2 = c^2, a^2 + a c cos(β) = b^2

    # Construct basis for rhombohedral unit cell
    r_a = 1.0
    r_α = π / 5

    expected_rhombohedral_unit_cell = RhombohedralUnitCell(r_a, r_α)
    r_basis_a, r_basis_b, r_basis_c = basis(expected_rhombohedral_unit_cell)

    # Construct monoclinic unit cell
    basis_a = r_basis_a
    basis_b = -r_basis_b + r_basis_c
    basis_c = r_basis_a - r_basis_b - r_basis_c
    unit_cell = MonoclinicUnitCell(
        norm(basis_a),
        norm(basis_b),
        norm(basis_c),
        acos(dot(basis_a, basis_c) / norm(basis_a) / norm(basis_c));
        centering=body_centering,
    )

    # Check test conditions
    @test unit_cell isa MonoclinicUnitCell
    @test centering(unit_cell) === body_centering
    @test lattice_constants(unit_cell).a^2 + lattice_constants(unit_cell).b^2 ≈
        lattice_constants(unit_cell).c^2
    @test lattice_constants(unit_cell).a^2 +
          lattice_constants(unit_cell).a *
          lattice_constants(unit_cell).c *
          cos(lattice_constants(unit_cell).β) ≈ lattice_constants(unit_cell).b^2

    # Exercise functionality
    iucr_unit_cell = conventional_cell(unit_cell)

    # Check results
    expected_unit_cell = standardize(RhombohedralUnitCell(r_a, r_α))
    @test iucr_unit_cell ≈ expected_unit_cell

    # ------ π/3 < r_α < π/2: a^2 + b^2 = c^2, b^2 + a c cos(β) = a^2

    # Construct basis for rhombohedral unit cell
    r_a = 1.0
    r_α = 2π / 5

    expected_rhombohedral_unit_cell = RhombohedralUnitCell(r_a, r_α)
    r_basis_a, r_basis_b, r_basis_c = basis(expected_rhombohedral_unit_cell)

    # Construct monoclinic unit cell
    basis_a = r_basis_a
    basis_b = r_basis_b - r_basis_c
    basis_c = -r_basis_a + r_basis_b + r_basis_c
    unit_cell = MonoclinicUnitCell(
        norm(basis_a),
        norm(basis_b),
        norm(basis_c),
        acos(dot(basis_a, basis_c) / norm(basis_a) / norm(basis_c));
        centering=body_centering,
    )

    # Check test conditions
    @test unit_cell isa MonoclinicUnitCell
    @test centering(unit_cell) === body_centering
    @test lattice_constants(unit_cell).a^2 + lattice_constants(unit_cell).b^2 ≈
        lattice_constants(unit_cell).c^2
    @test lattice_constants(unit_cell).b^2 +
          lattice_constants(unit_cell).a *
          lattice_constants(unit_cell).c *
          cos(lattice_constants(unit_cell).β) ≈ lattice_constants(unit_cell).a^2

    # Exercise functionality
    iucr_unit_cell = conventional_cell(unit_cell)

    # Check results
    expected_unit_cell = standardize(RhombohedralUnitCell(r_a, r_α))
    @test iucr_unit_cell ≈ expected_unit_cell

    # ------ π/2 < r_α < acos(-1/3):  a^2 + b^2 = c^2, b^2 + a c cos(β) = a^2

    # Construct basis for rhombohedral unit cell
    r_a = 1.0
    r_α = 3π / 5

    expected_rhombohedral_unit_cell = RhombohedralUnitCell(r_a, r_α)
    r_basis_a, r_basis_b, r_basis_c = basis(expected_rhombohedral_unit_cell)

    # Construct monoclinic unit cell
    basis_a = -r_basis_a
    basis_b = -r_basis_b + r_basis_c
    basis_c = r_basis_a + r_basis_b + r_basis_c
    unit_cell = MonoclinicUnitCell(
        norm(basis_a),
        norm(basis_b),
        norm(basis_c),
        acos(dot(basis_a, basis_c) / norm(basis_a) / norm(basis_c));
        centering=body_centering,
    )

    # Check test conditions
    @test unit_cell isa MonoclinicUnitCell
    @test centering(unit_cell) === body_centering
    @test lattice_constants(unit_cell).c^2 + 3 * lattice_constants(unit_cell).b^2 ≈
        9 * lattice_constants(unit_cell).a^2
    @test lattice_constants(unit_cell).c ≈
        -3 * lattice_constants(unit_cell).a * cos(lattice_constants(unit_cell).β)

    # Exercise functionality
    iucr_unit_cell = conventional_cell(unit_cell)

    # Check results
    expected_unit_cell = standardize(RhombohedralUnitCell(r_a, r_α))
    @test iucr_unit_cell ≈ expected_unit_cell

    # ------ r_α > acos(-1/3): a^2 + 3 b^2 = 9 c^2, a = -3 c cos(β)

    # Construct basis for rhombohedral unit cell
    r_a = 1.0
    r_α = 1.05 * acos(-1 / 3)

    expected_rhombohedral_unit_cell = RhombohedralUnitCell(r_a, r_α)
    r_basis_a, r_basis_b, r_basis_c = basis(expected_rhombohedral_unit_cell)

    # Construct monoclinic unit cell
    basis_a = r_basis_a + r_basis_b + r_basis_c
    basis_b = r_basis_b - r_basis_c
    basis_c = -r_basis_a
    unit_cell = MonoclinicUnitCell(
        norm(basis_a),
        norm(basis_b),
        norm(basis_c),
        acos(dot(basis_a, basis_c) / norm(basis_a) / norm(basis_c));
        centering=body_centering,
    )

    # Check test conditions
    @test unit_cell isa MonoclinicUnitCell
    @test centering(unit_cell) === body_centering
    @test lattice_constants(unit_cell).a^2 + 3 * lattice_constants(unit_cell).b^2 ≈
        9 * lattice_constants(unit_cell).c^2
    @test lattice_constants(unit_cell).a ≈
        -3 * lattice_constants(unit_cell).c * cos(lattice_constants(unit_cell).β)

    # Exercise functionality
    iucr_unit_cell = conventional_cell(unit_cell)

    # Check results
    expected_unit_cell = standardize(RhombohedralUnitCell(r_a, r_α))
    @test iucr_unit_cell ≈ expected_unit_cell
end

@testset "conventional_cell()::MonoclinicUnitCell: chain of limiting cases" begin
    # --- Exercise functionality and check results

    # primitive unit cell: aP --> mP
    a = 5
    b = 10
    c = 7
    β = 5π / 8
    monoclinic_unit_cell = MonoclinicUnitCell(a, b, c, β)
    basis_a, basis_b, basis_c = basis(monoclinic_unit_cell)

    triclinic_unit_cell = UnitCell(
        basis_a,
        basis_b,
        basis_c;
        identify_lattice_system=false,
        centering=primitive_centering,
    )
    expected_unit_cell = standardize(monoclinic_unit_cell)
    @test triclinic_unit_cell isa TriclinicUnitCell
    @test expected_unit_cell isa MonoclinicUnitCell
    @debug "chain of limiting cases: aP --> mP"
    @test conventional_cell(triclinic_unit_cell) ≈ expected_unit_cell

    # body-centered unit cell: aP --> mI
    #
    # Case #1: m_basis_a and m_basis_b in triclinic basis
    a = 5
    b = 10
    c = 7
    β = 5π / 8
    monoclinic_unit_cell = MonoclinicUnitCell(a, b, c, β; centering=body_centering)
    basis_a, basis_b, basis_c = basis(monoclinic_unit_cell)

    triclinic_unit_cell = UnitCell(
        basis_a,
        basis_b,
        0.5 * (basis_a + basis_b + basis_c);
        identify_lattice_system=false,
        centering=primitive_centering,
    )
    expected_unit_cell = standardize(monoclinic_unit_cell)
    @test triclinic_unit_cell isa TriclinicUnitCell
    @test expected_unit_cell isa MonoclinicUnitCell
    @debug "chain of limiting cases: aP --> mI"
    @test conventional_cell(triclinic_unit_cell) ≈ expected_unit_cell

    # body-centered unit cell: aP --> mI
    #
    # Case #1: m_basis_a and m_basis_c in triclinic basis
    a = 5
    b = 10
    c = 7
    β = 5π / 8
    monoclinic_unit_cell = MonoclinicUnitCell(a, b, c, β; centering=body_centering)
    basis_a, basis_b, basis_c = basis(monoclinic_unit_cell)

    triclinic_unit_cell = UnitCell(
        basis_a,
        basis_c,
        0.5 * (basis_a + basis_b + basis_c);
        identify_lattice_system=false,
        centering=primitive_centering,
    )
    expected_unit_cell = standardize(monoclinic_unit_cell)
    @test triclinic_unit_cell isa TriclinicUnitCell
    @test expected_unit_cell isa MonoclinicUnitCell
    @debug "chain of limiting cases: aP --> mI"
    @test conventional_cell(triclinic_unit_cell) ≈ expected_unit_cell

    # base-centered unit cell: aP --> mC
    #
    # Case #1: m_basis_a and m_basis_c in triclinic basis
    a = 5
    b = 10
    c = 7
    β = 5π / 8
    monoclinic_unit_cell = MonoclinicUnitCell(a, b, c, β; centering=base_centering)
    basis_a, basis_b, basis_c = basis(monoclinic_unit_cell)

    triclinic_unit_cell = UnitCell(
        basis_a,
        0.5 * (basis_a + basis_b),
        basis_c;
        identify_lattice_system=false,
        centering=primitive_centering,
    )
    expected_unit_cell = standardize(monoclinic_unit_cell)
    @test triclinic_unit_cell isa TriclinicUnitCell
    @test expected_unit_cell isa MonoclinicUnitCell
    @test conventional_cell(triclinic_unit_cell) ≈ expected_unit_cell

    # base-centered unit cell: aP --> mC
    #
    # Case #2: m_basis_a and two base-centered lattice vectors in triclinic basis
    a = 5
    b = 10
    c = 7
    β = 5π / 8
    monoclinic_unit_cell = MonoclinicUnitCell(a, b, c, β; centering=base_centering)
    basis_a, basis_b, basis_c = basis(monoclinic_unit_cell)

    triclinic_unit_cell = UnitCell(
        basis_a,
        0.5 * (basis_a + basis_b),
        0.5 * (basis_a + basis_b) + basis_c;
        identify_lattice_system=false,
        centering=primitive_centering,
    )
    expected_unit_cell = standardize(monoclinic_unit_cell)
    @test triclinic_unit_cell isa TriclinicUnitCell
    @test expected_unit_cell isa MonoclinicUnitCell
    @debug "chain of limiting cases: aP --> mC"
    @test conventional_cell(triclinic_unit_cell) ≈ expected_unit_cell
end
