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
Tests for unit_cell/unit_cell_symmetry.jl
"""
# --- Imports

# Standard library
using Test

# External packages
using Xtallography

# --- Tests

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
