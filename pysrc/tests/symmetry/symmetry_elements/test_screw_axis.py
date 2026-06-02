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
Unit tests for `ScrewAxis` class
"""
# --- Imports

# Standard library
from dataclasses import FrozenInstanceError
from fractions import Fraction
import unittest

# External packages
import juliacall
import pytest
from xtallography.symmetry import ScrewAxis

# Local packages/modules


# --- Test Suites


class test_ScrewAxis(unittest.TestCase):
    """
    Test suite for the `ScrewAxis` class
    """

    # --- Fixtures

    def setUp(self):
        """
        Prepare for test.
        """
        # Initialize JuliaCall
        self.jl = juliacall.newmodule("PyXtallographyTest")
        self.jl.seval("using Xtallography")

    def tearDown(self):
        """
        Clean up after test.
        """

    # --- Tests

    @staticmethod
    def test_init():
        """
        Test `__init__()`.
        """
        # --- basic usage

        # default location
        n = 6
        m = 4
        direction = (1, 0, 0)
        screw_axis = ScrewAxis(n, m, direction)

        assert screw_axis.n == n
        assert screw_axis.m == m
        assert screw_axis.direction == direction
        assert all([isinstance(x, Fraction)] for x in screw_axis.direction)
        assert screw_axis.location == (0, 0, 0)
        assert all([isinstance(x, Fraction)] for x in screw_axis.location)

        # --- keyword arguments

        # non-default location
        n = 6
        m = 4
        direction = (1, 0, 0)
        location = (0, 1, 0)
        screw_axis = ScrewAxis(n, m, direction, location=location)

        assert screw_axis.n == n
        assert screw_axis.m == m
        assert screw_axis.direction == direction
        assert all([isinstance(x, Fraction)] for x in screw_axis.direction)
        assert screw_axis.location == location
        assert all([isinstance(x, Fraction)] for x in screw_axis.location)

    @staticmethod
    def test_init_invalid_arguments():
        """
        Test argument checks for `__init__()`
        """
        # --- Tests

        # n = 0
        with pytest.raises(ValueError) as exception_info:
            ScrewAxis(0, 4, (1, 0, 0))

        expected_error = "`n` must be positive (n=0)"
        assert expected_error in str(exception_info)

        # n < 0
        with pytest.raises(ValueError) as exception_info:
            ScrewAxis(-5, 4, (1, 0, 0))

        expected_error = "`n` must be positive (n=-5)"
        assert expected_error in str(exception_info)

        # m = 0
        with pytest.raises(ValueError) as exception_info:
            ScrewAxis(6, 0, (1, 0, 0))

        expected_error = "`m` must be positive (m=0)"
        assert expected_error in str(exception_info)

        # m < 0
        with pytest.raises(ValueError) as exception_info:
            ScrewAxis(6, -4, (1, 0, 0))

        expected_error = "`m` must be positive (m=-4)"
        assert expected_error in str(exception_info)

        # m = n
        with pytest.raises(ValueError) as exception_info:
            ScrewAxis(6, 6, (1, 0, 0))

        expected_error = "`m` must be less than `n` (n=6,m=6)"
        assert expected_error in str(exception_info)

        # m > n
        with pytest.raises(ValueError) as exception_info:
            ScrewAxis(6, 7, (1, 0, 0))

        expected_error = "`m` must be less than `n` (n=6,m=7)"
        assert expected_error in str(exception_info)

        # --- direction

        # direction is not a 3-tuple
        with pytest.raises(ValueError) as exception_info:
            ScrewAxis(5, 3, (1, 2, 3, 4))

        expected_error = "`direction` must be a 3-tuple (direction=(1, 2, 3, 4))"
        assert expected_error in str(exception_info)

        # some element of direction is not convertible to a Fraction object
        invalid_direction = (1, "b", 3)
        with pytest.raises(ValueError) as exception_info:
            ScrewAxis(5, 3, invalid_direction)

        expected_error = (
            "all elements of `direction` must convertible to Fraction objects "
            f"(direction={invalid_direction}). "
            "[caused by"
        )

        assert expected_error in str(exception_info)

        # --- location

        # location is not a 3-tuple
        with pytest.raises(ValueError) as exception_info:
            ScrewAxis(5, 3, (1, 2, 3), location=(1, 2))

        expected_error = "`location` must be a 3-tuple (location=(1, 2))"
        assert expected_error in str(exception_info)

        # some element of location is not convertible to a Fraction object
        invalid_location = (1, 3, "b")
        with pytest.raises(ValueError) as exception_info:
            ScrewAxis(5, 3, (1, 2, 3), location=invalid_location)

        expected_error = (
            "all elements of `location` must convertible to Fraction objects "
            f"(location={invalid_location}). "
            "[caused by"
        )

        assert expected_error in str(exception_info)

    @staticmethod
    def test_frozen():
        """
        Test that ScrewAxis objects are immutable.
        """
        # --- Preparations

        screw_axis = ScrewAxis(6, 4, (1, 0, 0))

        # --- Tests

        # attempt to change `n` field
        with pytest.raises(FrozenInstanceError) as exception_info:
            screw_axis.n = 10

        expected_error = "cannot assign to field 'n'"
        assert expected_error in str(exception_info)

        # attempt to change `m` field
        with pytest.raises(FrozenInstanceError) as exception_info:
            screw_axis.m = 5

        expected_error = "cannot assign to field 'm'"
        assert expected_error in str(exception_info)

        # attempt to change `direction` field
        with pytest.raises(FrozenInstanceError) as exception_info:
            screw_axis.direction = (0, 1, 0)

        expected_error = "cannot assign to field 'direction'"
        assert expected_error in str(exception_info)

        # attempt to change `location` field
        with pytest.raises(FrozenInstanceError) as exception_info:
            screw_axis.location = (0, 1, 0)

        expected_error = "cannot assign to field 'location'"
        assert expected_error in str(exception_info)

    def test_to_julia(self):
        """
        Test `to_julia()`.
        """
        # --- Preparations

        n = 6
        m = 4
        direction = (1, 0, 0)
        screw_axis = ScrewAxis(n, m, direction)

        # --- Tests

        screw_axis_jl = screw_axis.to_julia()
        assert self.jl.isa(screw_axis_jl, self.jl.ScrewAxis)
        assert screw_axis_jl.n == n
        assert screw_axis_jl.m == m
        assert screw_axis_jl.direction == direction

    def test_from_julia(self):
        """
        Test `from_julia()`.
        """
        # --- default location

        n = 6
        m = 4
        direction = (1, 2, 3)
        screw_axis_jl = self.jl.ScrewAxis(n, m, direction)
        screw_axis = ScrewAxis.from_julia(screw_axis_jl)

        assert screw_axis == ScrewAxis(n, m, direction, location=(0, 0, 0))

        # --- non-default location

        n = 6
        m = 3
        direction = (1, 2, 3)
        location = (4, 5, 6)
        screw_axis_jl = self.jl.ScrewAxis(n, m, direction, location=location)
        screw_axis = ScrewAxis.from_julia(screw_axis_jl)

        assert screw_axis == ScrewAxis(n, m, direction, location=location)

    def test_from_julia_invalid_args(self):
        """
        Test argument checks for `from_julia()`
        """
        screw_axis_jl = 10
        with pytest.raises(ValueError) as exception_info:
            ScrewAxis.from_julia(screw_axis_jl)

        expected_error = (
            "`screw_axis_jl` must be a Julia `ScrewAxis` object. "
            f"(screw_axis_jl={screw_axis_jl})."
        )
        assert expected_error in str(exception_info)

    def test_repr(self):
        """
        Test `__repr__()`.
        """
        # default location
        n = 6
        m = 4
        direction = (1, 0, 0)
        screw_axis = ScrewAxis(n, m, direction)

        assert str(screw_axis) == (
            "ScrewAxis(n=6,m=4,"
            "direction=(Fraction(1, 1), Fraction(0, 1), Fraction(0, 1)),"
            "location=(Fraction(0, 1), Fraction(0, 1), Fraction(0, 1)))"
        )

        # non-default location
        n = 6
        m = 3
        direction = (1, 0, 0)
        location = (1, 0, 1)
        screw_axis = ScrewAxis(n, m, direction, location=location)

        assert str(screw_axis) == (
            "ScrewAxis(n=6,m=3,"
            "direction=(Fraction(1, 1), Fraction(0, 1), Fraction(0, 1)),"
            "location=(Fraction(1, 1), Fraction(0, 1), Fraction(1, 1)))"
        )

    def test_eq(self):
        """
        Test `__eq__()`.
        """
        # --- Identical rotation axes

        symmetry_element_1 = ScrewAxis(6, 4, (1, 0, 0), (0, 0, 0))
        symmetry_element_2 = ScrewAxis(6, 4, (1, 0, 0), (0, 0, 0))

        assert symmetry_element_1 == symmetry_element_2

        # --- Equivalent rotation axes

        # directions differ, locations same
        symmetry_element_1 = ScrewAxis(6, 4, (1, 0, 0), (0, 0, 0))
        symmetry_element_2 = ScrewAxis(6, 4, (0.5, 0, 0), (0, 0, 0))

        assert symmetry_element_1 == symmetry_element_2

        # directions same, locations differ
        symmetry_element_1 = ScrewAxis(6, 4, (1, 0, 0), (0, 0, 0))
        symmetry_element_2 = ScrewAxis(6, 4, (1, 0, 0), (0.75, 0, 0))

        assert symmetry_element_1 == symmetry_element_2

        # directions differ, locations differ
        symmetry_element_1 = ScrewAxis(6, 4, (Fraction(2, 3), 0, 0), (0, 0, 0))
        symmetry_element_2 = ScrewAxis(6, 4, (1, 0, 0), (Fraction(3, 2), 0, 0))

        assert symmetry_element_1 == symmetry_element_2

        # --- Inequivalent rotation axes

        # orders differ
        symmetry_element_1 = ScrewAxis(6, 4, (2, 0, 0), (2, 0, 0))
        symmetry_element_2 = ScrewAxis(5, 4, (1, 0, 0), (0, 0, 0))

        assert symmetry_element_1 != symmetry_element_2

        # number of translation steps differ
        symmetry_element_1 = ScrewAxis(6, 4, (2, 0, 0), (2, 0, 0))
        symmetry_element_2 = ScrewAxis(6, 3, (1, 0, 0), (0, 0, 0))

        assert symmetry_element_1 != symmetry_element_2

        # directions differ
        symmetry_element_1 = ScrewAxis(6, 4, (1, 1, 0), (0, 0, 0))
        symmetry_element_2 = ScrewAxis(6, 4, (1, 0, 0), (0, 0, 0))

        assert symmetry_element_1 != symmetry_element_2

        # line between locations != direction
        symmetry_element_1 = ScrewAxis(6, 4, (Fraction(1, 3), 0, 0), (0, 1, 0))
        symmetry_element_2 = ScrewAxis(6, 4, (1, 0, 0), (0, 0, 0))

        assert symmetry_element_1 != symmetry_element_2
