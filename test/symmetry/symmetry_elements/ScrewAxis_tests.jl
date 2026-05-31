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
Tests for ScrewAxis type (defined in symmetry/symmetry_elements.jl)
"""
# --- Imports

# Standard library
using Test

# External packages
using Xtallography

# --- Tests

@testset "ScrewAxis: inner constructor" begin
    # --- Valid arguments

    n = 2
    m = 1
    direction = (1, 1, 1)
    location = (0, 0, 0)
    symmetry_element = ScrewAxis(n, m, direction, location)

    @test symmetry_element.n == n
    @test symmetry_element.m == m
    @test symmetry_element.direction == direction
    @test symmetry_element.location == location

    # --- Invalid arguments

    # n = 0
    n = 0
    m = 4
    direction = (1, 1, 1)
    location = (1, 0, 0)

    expected_message = "`n` must be positive (n=0)"
    @test_throws ArgumentError(expected_message) symmetry_element = ScrewAxis(
        n, m, direction, location
    )

    # n < 0
    n = -10
    m = 4
    direction = (1, 1, 1)
    location = (1, 0, 0)

    expected_message = "`n` must be positive (n=-10)"
    @test_throws ArgumentError(expected_message) symmetry_element = ScrewAxis(
        n, m, direction, location
    )

    # m = 0
    n = 2
    m = 0
    direction = (1, 1, 1)
    location = (1, 0, 0)

    expected_message = "`m` must be positive (m=0)"
    @test_throws ArgumentError(expected_message) symmetry_element = ScrewAxis(
        n, m, direction, location
    )

    # m < 0
    n = 2
    m = -1
    direction = (1, 1, 1)
    location = (1, 0, 0)

    expected_message = "`m` must be positive (m=-1)"
    @test_throws ArgumentError(expected_message) symmetry_element = ScrewAxis(
        n, m, direction, location
    )

    # m = n
    n = 3
    m = 3
    direction = (1, 1, 1)
    location = (1, 0, 0)

    expected_message = "`m` must be less than `n` (n=3,m=3)"
    @test_throws ArgumentError(expected_message) symmetry_element = ScrewAxis(
        n, m, direction, location
    )

    # m > n
    n = 3
    m = 4
    direction = (1, 1, 1)
    location = (1, 0, 0)

    expected_message = "`m` must be less than `n` (n=3,m=4)"
    @test_throws ArgumentError(expected_message) symmetry_element = ScrewAxis(
        n, m, direction, location
    )
end

@testset "ScrewAxis: outer constructor" begin
    # --- default location

    n = 3
    m = 2
    direction = (1, 0, 0)

    symmetry_element = ScrewAxis(n, m, direction)

    @test symmetry_element.n == n
    @test symmetry_element.m == m
    @test symmetry_element.direction == direction
    @test symmetry_element.location == (0, 0, 0)

    # --- non-default location

    n = 3
    m = 1
    direction = (1, 0, 0)
    location = (2, 0, 0)

    symmetry_element = ScrewAxis(n, m, direction; location=location)

    @test symmetry_element.n == n
    @test symmetry_element.m == m
    @test symmetry_element.direction == direction
    @test symmetry_element.location == location
end

@testset ":(==)(::ScrewAxis,::ScrewAxis)" begin
    # --- Identical screw axes

    symmetry_element_1 = ScrewAxis(3, 2, (1, 0, 0), (0, 0, 0))
    symmetry_element_2 = ScrewAxis(3, 2, (1, 0, 0), (0, 0, 0))

    @test symmetry_element_1 == symmetry_element_2

    # --- Equivalent screw axes

    # directions differ, locations same
    symmetry_element_1 = ScrewAxis(3, 2, (1, 0, 0), (0, 0, 0))
    symmetry_element_2 = ScrewAxis(3, 2, (1//2, 0, 0), (0, 0, 0))

    @test symmetry_element_1 == symmetry_element_2

    # directions same, locations differ
    symmetry_element_1 = ScrewAxis(3, 2, (1, 0, 0), (0, 0, 0))
    symmetry_element_2 = ScrewAxis(3, 2, (1, 0, 0), (3//4, 0, 0))

    @test symmetry_element_1 == symmetry_element_2

    # directions differ, locations differ
    symmetry_element_1 = ScrewAxis(3, 2, (2//3, 0, 0), (0, 0, 0))
    symmetry_element_2 = ScrewAxis(3, 2, (1, 0, 0), (3//2, 0, 0))

    @test symmetry_element_1 == symmetry_element_2

    # --- Inequivalent screw axes

    # rotation orders differ
    symmetry_element_1 = ScrewAxis(4, 2, (2, 0, 0), (0, 0, 0))
    symmetry_element_2 = ScrewAxis(3, 2, (1, 0, 0), (0, 0, 0))

    @test symmetry_element_1 != symmetry_element_2

    # translation steps differ
    symmetry_element_1 = ScrewAxis(3, 2, (2, 0, 0), (0, 0, 0))
    symmetry_element_2 = ScrewAxis(3, 1, (1, 0, 0), (0, 0, 0))

    @test symmetry_element_1 != symmetry_element_2

    # directions differ
    symmetry_element_1 = ScrewAxis(3, 2, (1, 1, 0), (0, 0, 0))
    symmetry_element_2 = ScrewAxis(3, 2, (1, 0, 0), (0, 0, 0))

    @test symmetry_element_1 != symmetry_element_2

    # line between locations != direction
    symmetry_element_1 = ScrewAxis(3, 2, (1//3, 0, 0), (0, 1, 0))
    symmetry_element_2 = ScrewAxis(3, 2, (1, 0, 0), (0, 0, 0))

    @test symmetry_element_1 != symmetry_element_2
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
