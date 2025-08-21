#   Copyright 2025 Velexi Corporation
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
Tests for symmetry/symmetry_elements.jl
"""
# --- Imports

# Standard library
using Test

# External packages
using Xtallography

# --- Tests

# ------ Symmetry element types

@testset "SymmetryElement Subtypes" begin
    expected_types = [GlidePlane, ScrewAxis]

    for type in expected_types
        @test type <: SymmetryElement
    end
end

@testset "GlidePlane constructor" begin
    # --- Tests

    # ----- Valid arguments

    translation = "1,1,1"
    reflection_plane = "0,1,0"
    symmetry_element = GlidePlane(translation, reflection_plane)

    @test symmetry_element.translation == translation
    @test symmetry_element.reflection_plane == reflection_plane
end

@testset "GlidePlane constants" begin
    expected_constants = [
        b_perp_a,
        c_perp_a,
        n_perp_a,
        d_perp_a,
        a_perp_b,
        c_perp_b,
        n_perp_b,
        d_perp_b,
        a_perp_c,
        b_perp_c,
        n_perp_c,
        d_perp_c,
    ]
    for constant in expected_constants
        @test constant isa GlidePlane
    end
end

@testset "ScrewAxis constructor" begin
    # --- Tests

    # ----- Valid arguments

    axis = "1,1,1"
    n = 2
    m = 1
    symmetry_element = ScrewAxis(axis, n, m)

    @test symmetry_element.axis == axis
    @test symmetry_element.n == n
    @test symmetry_element.m == m

    # ----- Invalid arguments

    axis = "1,1,1"
    n = 2
    m = 4

    expected_message = "`m` be no greater than `n` (n=2,m=4)"
    @test_throws ArgumentError(expected_message) symmetry_element = ScrewAxis(axis, n, m)
end

@testset "ScrewAxis constants" begin
    expected_constants = [
        a_2_1,
        a_4_1,
        a_4_2,
        a_4_3,
        b_2_1,
        b_4_1,
        b_4_2,
        b_4_3,
        c_2_1,
        c_3_1,
        c_3_2,
        c_4_1,
        c_4_2,
        c_4_3,
        c_6_1,
        c_6_2,
        c_6_3,
        c_6_4,
        c_6_5,
    ]

    for constant in expected_constants
        @test constant isa ScrewAxis
    end
end

# ------ Unit cell symmetry

@testset "UnitCellSymmetry inner constructor" begin
    # --- Tests

    # empty symmetry elements
    centering = face_centering
    symmetry_elements = Vector{SymmetryElement}()

    unit_cell_symmetry = UnitCellSymmetry(centering, symmetry_elements)

    @test unit_cell_symmetry.centering == centering
    @test unit_cell_symmetry.symmetry_elements == symmetry_elements

    # non-empty symmetry elements
    centering = body_centering
    symmetry_elements = Vector{SymmetryElement}()
    push!(symmetry_elements, c_perp_a)
    push!(symmetry_elements, c_4_1)
    push!(symmetry_elements, d_perp_c)

    unit_cell_symmetry = UnitCellSymmetry(centering, symmetry_elements)

    @test unit_cell_symmetry.centering == centering
    @test Set(unit_cell_symmetry.symmetry_elements) == Set(symmetry_elements)
end

@testset "UnitCellSymmetry outer constructor: UnitCellSymmetry(::centering;::Union{Vector,Nothing})" begin
    # --- Tests

    # default symmetry elements
    centering = face_centering

    unit_cell_symmetry = UnitCellSymmetry(centering)

    @test unit_cell_symmetry.centering == centering
    @test unit_cell_symmetry.symmetry_elements isa Vector
    @test isempty(unit_cell_symmetry.symmetry_elements)

    # symmetry_elements = nothing
    centering = face_centering
    symmetry_elements = nothing

    unit_cell_symmetry = UnitCellSymmetry(centering; symmetry_elements=symmetry_elements)

    @test unit_cell_symmetry.centering == centering
    @test unit_cell_symmetry.symmetry_elements isa Vector
    @test isempty(unit_cell_symmetry.symmetry_elements)

    # non-default, non-nothing symmetry_elements
    centering = base_centering
    symmetry_elements = [c_perp_a, c_4_1, d_perp_c]

    unit_cell_symmetry = UnitCellSymmetry(centering; symmetry_elements=symmetry_elements)

    @test unit_cell_symmetry.centering == centering
    @test Set(unit_cell_symmetry.symmetry_elements) == Set(symmetry_elements)
end

@testset "UnitCellSymmetry outer constructor: UnitCellSymmetry()" begin
    unit_cell_symmetry = UnitCellSymmetry()
    @test unit_cell_symmetry == primitive_unit_cell_symmetry
end

@testset "UnitCellSymmetry constants" begin
    expected_constants = [primitive_unit_cell_symmetry]

    for constant in expected_constants
        @test constant isa UnitCellSymmetry
    end
end

@testset ":(==)(::UnitCellSymmetry)" begin
    # --- Tests

    # x == y
    centering = primitive_centering
    symmetry_elements = [c_perp_a, c_4_1, d_perp_c]
    x = UnitCellSymmetry(centering; symmetry_elements=symmetry_elements)
    y = UnitCellSymmetry(centering; symmetry_elements=symmetry_elements)

    @test x == y
    @test x !== y

    # x.centering == y.centering
    symmetry_elements = [c_perp_a, c_4_1, d_perp_c]
    x = UnitCellSymmetry(face_centering; symmetry_elements=symmetry_elements)
    y = UnitCellSymmetry(body_centering; symmetry_elements=symmetry_elements)

    @test x != y

    # x.symmetry_elements and y.symmetry_elements contain different elements
    centering = base_centering
    symmetry_elements_x = [c_perp_a, c_4_1, d_perp_c]
    symmetry_elements_y = [c_perp_a, c_4_3]
    x = UnitCellSymmetry(centering; symmetry_elements=symmetry_elements_x)
    y = UnitCellSymmetry(centering; symmetry_elements=symmetry_elements_y)

    @test x != y

    # x.symmetry_elements and y.symmetry_elements contain same elements in different order
    centering = base_centering
    symmetry_elements_x = [c_perp_a, c_4_1, d_perp_c]
    symmetry_elements_y = [d_perp_c, c_perp_a, c_4_1]
    x = UnitCellSymmetry(centering; symmetry_elements=symmetry_elements_x)
    y = UnitCellSymmetry(centering; symmetry_elements=symmetry_elements_y)

    @test x == y
    @test x !== y
end
