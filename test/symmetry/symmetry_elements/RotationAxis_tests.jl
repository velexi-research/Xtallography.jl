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
Tests for RotationAxis type (defined in symmetry/symmetry_elements.jl)
"""
# --- Imports

# Standard library
using Test

# External packages
using Xtallography

# --- Tests

@testset "RotationAxis: inner constructor" begin
    # --- Valid arguments

    # normal usage
    n = 3
    direction = (1, 0, 0)
    location = (0, 0, 0)

    symmetry_element = RotationAxis(n, direction, location)

    @test symmetry_element.n == n
    @test symmetry_element.direction == direction
    @test symmetry_element.location == location

    # --- Invalid arguments

    # n = 0
    n = 0
    direction = (1, 0, 0)
    location = (0, 0, 0)

    expected_message = "`n` must be positive (n=0)"
    @test_throws ArgumentError(expected_message) RotationAxis(n, direction, location)

    # n < 0
    n = -2
    direction = (1, 0, 0)
    location = (0, 0, 0)

    expected_message = "`n` must be positive (n=-2)"
    @test_throws ArgumentError(expected_message) RotationAxis(n, direction, location)

    # direction does not contain exactly 3 elements
    n = 4
    direction = (1, 0, 0, 0)
    location = (0, 0, 0)

    @test_throws MethodError RotationAxis(n, direction, location)

    # location does not contain exactly 3 elements
    n = 4
    direction = (1, 0, 0)
    location = (1, 0)

    @test_throws MethodError RotationAxis(n, direction, location)
end

@testset "RotationAxis: outer constructor" begin
    # --- default location 

    n = 3
    direction = (1, 0, 0)

    symmetry_element = RotationAxis(n, direction)

    @test symmetry_element.n == n
    @test symmetry_element.direction == direction
    @test symmetry_element.location == (0, 0, 0)

    # --- non-default location 

    n = 3
    direction = (1, 0, 0)
    location = (2, 0, 0)

    symmetry_element = RotationAxis(n, direction; location=location)

    @test symmetry_element.n == n
    @test symmetry_element.direction == direction
    @test symmetry_element.location == location
end

@testset ":(==)(::RotationAxis,::RotationAxis)" begin
    # --- Identical rotation axes

    symmetry_element_1 = RotationAxis(2, (1, 0, 0), (0, 0, 0))
    symmetry_element_2 = RotationAxis(2, (1, 0, 0), (0, 0, 0))

    @test symmetry_element_1 == symmetry_element_2

    # --- Equivalent rotation axes: directions differ, locations same

    symmetry_element_1 = RotationAxis(2, (1, 0, 0), (0, 0, 0))
    symmetry_element_2 = RotationAxis(2, (1//2, 0, 0), (0, 0, 0))

    @test symmetry_element_1 == symmetry_element_2

    # --- Equivalent rotation axes: directions same, locations differ

    symmetry_element_1 = RotationAxis(2, (1, 0, 0), (0, 0, 0))
    symmetry_element_2 = RotationAxis(2, (1, 0, 0), (3//4, 0, 0))

    @test symmetry_element_1 == symmetry_element_2

    # --- Equivalent rotation axes: directions differ, locations differ

    symmetry_element_1 = RotationAxis(2, (2//3, 0, 0), (0, 0, 0))
    symmetry_element_2 = RotationAxis(2, (1, 0, 0), (3//2, 0, 0))

    @test symmetry_element_1 == symmetry_element_2

    # --- Inequivalent rotation axes: directions differ

    symmetry_element_1 = RotationAxis(2, (1, 1, 0), (0, 0, 0))
    symmetry_element_2 = RotationAxis(2, (1, 0, 0), (0, 0, 0))

    @test symmetry_element_1 != symmetry_element_2

    # --- Inequivalent rotation axes: line between locations != direction

    symmetry_element_1 = RotationAxis(2, (1//3, 0, 0), (0, 1, 0))
    symmetry_element_2 = RotationAxis(2, (1, 0, 0), (0, 0, 0))

    @test symmetry_element_1 != symmetry_element_2

    # --- Inequivalent rotation axes: order is different

    symmetry_element_1 = RotationAxis(2, (2, 0, 0), (2, 0, 0))
    symmetry_element_2 = RotationAxis(3, (1, 0, 0), (0, 0, 0))

    @test symmetry_element_1 != symmetry_element_2
end
