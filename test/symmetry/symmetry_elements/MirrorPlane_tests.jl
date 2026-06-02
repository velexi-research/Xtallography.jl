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
Tests for MirrorPlane type (defined in symmetry/symmetry_elements.jl)
"""
# --- Imports

# Standard library
using Test

# External packages
using Xtallography

# --- Tests

@testset "MirrorPlane: inner constructor" begin
    # --- Valid arguments

    normal = (0, 1, 0)
    location = (1, 2, 3)
    symmetry_element = MirrorPlane(normal, location)

    @test symmetry_element.normal == normal
    @test symmetry_element.location == location

    # --- Invalid arguments

    # normal = (0,0,0)
    normal = (0, 0, 0)
    location = (1, 2, 3)

    expected_message = "`normal` must be a nonzero vector (normal=$normal)"

    @test_throws ArgumentError(expected_message) MirrorPlane(normal, location)

    # normal does not contain exactly 3 elements
    normal = (1, 0, 0, 0)
    location = (0, 0, 0)

    @test_throws MethodError MirrorPlane(normal, location)

    # location does not contain exactly 3 elements
    normal = (1, 0, 0)
    location = (1, 0)

    @test_throws MethodError MirrorPlane(normal, location)
end

@testset "MirrorPlane: outer constructor" begin
    # --- default location

    normal = (0, 1, 0)

    symmetry_element = MirrorPlane(normal)

    @test symmetry_element.normal == normal
    @test symmetry_element.location == (0, 0, 0)

    # --- non-default location

    normal = (0, 2, 0)
    location = (1, 2, 3)

    symmetry_element = MirrorPlane(normal; location=location)

    @test symmetry_element.normal == normal
    @test symmetry_element.location == location
end

@testset ":(==)(::MirrorPlane,::MirrorPlane)" begin
    # --- Identical mirror plane

    symmetry_element_1 = MirrorPlane((1, 0, 0), (0, 0, 0))
    symmetry_element_2 = MirrorPlane((1, 0, 0), (0, 0, 0))

    @test symmetry_element_1 == symmetry_element_2

    # --- Equivalent mirror plane

    # normals differ, locations same
    symmetry_element_1 = MirrorPlane((1, 0, 0), (0, 0, 0))
    symmetry_element_2 = MirrorPlane((1//2, 0, 0), (0, 0, 0))

    @test symmetry_element_1 == symmetry_element_2

    # normals same, locations differ
    symmetry_element_1 = MirrorPlane((1, 0, 0), (0, 0, 0))
    symmetry_element_2 = MirrorPlane((1, 0, 0), (0, 3//4, 0))

    @test symmetry_element_1 == symmetry_element_2

    # normals differ, locations differ
    symmetry_element_1 = MirrorPlane((2//3, 0, 0), (0, 1//2, 0))
    symmetry_element_2 = MirrorPlane((1, 0, 0), (0, -2//7, 3//2))

    @test symmetry_element_1 == symmetry_element_2

    # --- Inequivalent mirror plane

    # normals differ
    symmetry_element_1 = MirrorPlane((1, 1, 0), (0, 0, 0))
    symmetry_element_2 = MirrorPlane((1, 0, 0), (0, 0, 0))

    @test symmetry_element_1 != symmetry_element_2

    # line between locations is not orthogonal to normal
    symmetry_element_1 = MirrorPlane((1//3, 0, 0), (1//6, 1, 0))
    symmetry_element_2 = MirrorPlane((1, 0, 0), (0, 0, 0))

    @test symmetry_element_1 != symmetry_element_2
end
