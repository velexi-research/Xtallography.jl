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
Tests for RotoinversionAxis type (defined in symmetry/symmetry_elements.jl)
"""
# --- Imports

# Standard library
using Test

# External packages
using Xtallography

# --- Tests

@testset "RotoinversionAxis constructor" begin
    # --- Valid arguments

    direction = (1, 0, 0)
    location = (0, 0, 0)
    n = 3

    symmetry_element = RotoinversionAxis(direction, location, n)

    @test symmetry_element.direction == direction
    @test symmetry_element.location == location
    @test symmetry_element.n == n

    # --- Invalid arguments

    # direction does not contain exactly 3 elements
    direction = (1, 0, 0, 0)
    location = (0, 0, 0)
    n = 4

    @test_throws MethodError RotoinversionAxis(direction, location, n)

    # location does not contain exactly 3 elements
    direction = (1, 0, 0)
    location = (1, 0)
    n = 4

    @test_throws MethodError RotoinversionAxis(direction, location, n)

    # n is not an integer
    direction = (1, 0, 0)
    location = (0, 0, 0)
    n = 4.2

    @test_throws InexactError(:Int64, Int64, 4.2) RotoinversionAxis(direction, location, n)
end

@testset ":(==)(::RotoinversionAxis,::RotoinversionAxis)" begin
    # --- Identical rotoinversion axes

    symmetry_element_1 = RotoinversionAxis((1, 0, 0), (0, 0, 0), 2)
    symmetry_element_2 = RotoinversionAxis((1, 0, 0), (0, 0, 0), 2)

    @test symmetry_element_1 == symmetry_element_2

    # --- Equivalent rotoinversion axes: directions differ, locations same

    symmetry_element_1 = RotoinversionAxis((1, 0, 0), (0, 0, 0), 2)
    symmetry_element_2 = RotoinversionAxis((1//2, 0, 0), (0, 0, 0), 2)

    @test symmetry_element_1 == symmetry_element_2

    # --- Equivalent rotoinversion axes: directions same, locations differ

    symmetry_element_1 = RotoinversionAxis((1, 0, 0), (0, 0, 0), 2)
    symmetry_element_2 = RotoinversionAxis((1, 0, 0), (3//4, 0, 0), 2)

    @test symmetry_element_1 == symmetry_element_2

    # --- Equivalent rotoinversion axes: directions differ, locations differ

    symmetry_element_1 = RotoinversionAxis((2//3, 0, 0), (0, 0, 0), 2)
    symmetry_element_2 = RotoinversionAxis((1, 0, 0), (3//2, 0, 0), 2)

    @test symmetry_element_1 == symmetry_element_2

    # --- Inequivalent rotoinversion axes: directions differ

    symmetry_element_1 = RotoinversionAxis((1, 1, 0), (0, 0, 0), 2)
    symmetry_element_2 = RotoinversionAxis((1, 0, 0), (0, 0, 0), 2)

    @test symmetry_element_1 != symmetry_element_2

    # --- Inequivalent rotoinversion axes: line between locations != direction

    symmetry_element_1 = RotoinversionAxis((1//3, 0, 0), (0, 1, 0), 2)
    symmetry_element_2 = RotoinversionAxis((1, 0, 0), (0, 0, 0), 2)

    @test symmetry_element_1 != symmetry_element_2

    # --- Inequivalent rotoinversion axes: order is different

    symmetry_element_1 = RotoinversionAxis((2, 0, 0), (2, 0, 0), 3)
    symmetry_element_2 = RotoinversionAxis((1, 0, 0), (0, 0, 0), 2)

    @test symmetry_element_1 != symmetry_element_2
end
