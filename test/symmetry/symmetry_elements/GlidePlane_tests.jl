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
Tests for GlidePlane type (defined in symmetry/symmetry_elements.jl)
"""
# --- Imports

# Standard library
using Test

# External packages
using Xtallography

# --- Tests

@testset "GlidePlane: inner constructor" begin
    # --- Valid arguments

    glide = (1, 0, 0)
    normal = (0, 1, 0)
    location = (1, 2, 3)
    symmetry_element = GlidePlane(glide, normal, location)

    @test symmetry_element.glide == glide
    @test symmetry_element.normal == normal
    @test symmetry_element.location == location

    # --- Invalid arguments

    glide = (1, 1, 1)
    normal = (0, 1, 0)
    location = (1, 2, 3)

    expected_message = "`glide` must be orthogonal to `normal` (glide=$glide,normal=$normal)"
    @test_throws ArgumentError(expected_message) GlidePlane(glide, normal, location)
end

@testset "GlidePlane: outer constructor" begin
    # --- default location 

    glide = (1, 0, 0)
    normal = (0, 1, 0)

    symmetry_element = GlidePlane(glide, normal)

    @test symmetry_element.glide == glide
    @test symmetry_element.normal == normal
    @test symmetry_element.location == (0, 0, 0)

    # --- non-default location 

    glide = (1, 0, 0)
    normal = (0, 2, 0)
    location = (1, 2, 3)

    symmetry_element = GlidePlane(glide, normal; location=location)

    @test symmetry_element.glide == glide
    @test symmetry_element.normal == normal
    @test symmetry_element.location == location
end

@testset ":(==)(::GlidePlane,::GlidePlane)" begin
    # --- Identical mirror plane

    symmetry_element_1 = GlidePlane((1, 0, 0), (0, 1, 0), (0, 0, 0))
    symmetry_element_2 = GlidePlane((1, 0, 0), (0, 1, 0), (0, 0, 0))

    @test symmetry_element_1 == symmetry_element_2

    # --- Equivalent mirror plane

    # glides differ, normals same, locations same
    symmetry_element_1 = GlidePlane((1, 0, 0), (0, 1, 0), (0, 0, 0))
    symmetry_element_2 = GlidePlane((1//2, 0, 0), (0, 1, 0), (0, 0, 0))

    @test symmetry_element_1 == symmetry_element_2

    # glides same, normals differ, locations same
    symmetry_element_1 = GlidePlane((1, 0, 0), (0, 1, 0), (0, 0, 0))
    symmetry_element_2 = GlidePlane((1//2, 0, 0), (0, 1, 0), (0, 0, 0))

    @test symmetry_element_1 == symmetry_element_2

    # glides same, normals same, locations differ
    symmetry_element_1 = GlidePlane((1, 0, 0), (0, 1, 0), (0, 0, 0))
    symmetry_element_2 = GlidePlane((1, 0, 0), (0, 1, 0), (3//4, 0, 0))

    @test symmetry_element_1 == symmetry_element_2

    # glides differ, normals differ, locations same
    symmetry_element_1 = GlidePlane((2, 0, 0), (0, 1, 0), (0, 0, 0))
    symmetry_element_2 = GlidePlane((1, 0, 0), (0, 1//2, 0), (0, 0, 0))

    @test symmetry_element_1 == symmetry_element_2

    # glides differ, normals same, locations differ
    symmetry_element_1 = GlidePlane((2, 0, 0), (0, 1, 0), (0, 0, 0))
    symmetry_element_2 = GlidePlane((1, 0, 0), (0, 1, 0), (3//2, 0, 0))

    @test symmetry_element_1 == symmetry_element_2

    # glides same, normals differ, locations differ
    symmetry_element_1 = GlidePlane((0, 0, 2//3), (2//7, 0, 0), (0, 1//2, 0))
    symmetry_element_2 = GlidePlane((0, 0, 2//3), (7//3, 0, 0), (0, -2//7, 3//2))

    @test symmetry_element_1 == symmetry_element_2

    # --- Inequivalent mirror plane

    # glides differ
    symmetry_element_1 = GlidePlane((1, 1, 0), (0, 0, 1), (0, 0, 0))
    symmetry_element_2 = GlidePlane((1, 0, 0), (0, 0, 1), (0, 0, 0))

    @test symmetry_element_1 != symmetry_element_2

    # normals differ
    symmetry_element_1 = GlidePlane((0, 0, 1), (1, 1, 0), (0, 0, 0))
    symmetry_element_2 = GlidePlane((0, 0, 1), (1, 0, 0), (0, 0, 0))

    @test symmetry_element_1 != symmetry_element_2

    # line between locations is not orthogonal to normal
    symmetry_element_1 = GlidePlane((0, 0, 1), (1//3, 0, 0), (1//6, 1, 0))
    symmetry_element_2 = GlidePlane((0, 0, 1), (1, 0, 0), (0, 0, 0))

    @test symmetry_element_1 != symmetry_element_2
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
