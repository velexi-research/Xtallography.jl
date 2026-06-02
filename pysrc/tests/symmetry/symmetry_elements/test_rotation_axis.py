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
Unit tests for `RotationAxis` class
"""
# --- Imports

# Standard library
from dataclasses import FrozenInstanceError
from fractions import Fraction
import unittest

# External packages
import juliacall
import pytest
from xtallography.symmetry import RotationAxis

# Local packages/modules


# --- Test Suites


class test_xtallography_symmetry_symmetry_elements_RotationAxis(unittest.TestCase):
    """
    Test suite for the `RotationAxis` class
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
        direction = (1, 0, 0)
        rotation_axis = RotationAxis(n, direction)

        assert rotation_axis.n == n
        assert rotation_axis.direction == direction
        assert all([isinstance(x, Fraction)] for x in rotation_axis.direction)
        assert rotation_axis.location == (0, 0, 0)
        assert all([isinstance(x, Fraction)] for x in rotation_axis.location)

        # --- keyword arguments

        # non-default location
        n = 6
        direction = (1, 0, 0)
        location = (0, 1, 0)
        rotation_axis = RotationAxis(n, direction, location=location)

        assert rotation_axis.n == n
        assert rotation_axis.direction == direction
        assert all([isinstance(x, Fraction)] for x in rotation_axis.direction)
        assert rotation_axis.location == location
        assert all([isinstance(x, Fraction)] for x in rotation_axis.location)

    @staticmethod
    def test_init_invalid_arguments():
        """
        Test argument checks for `__init__()`
        """
        # --- n

        # n = 0
        with pytest.raises(ValueError) as exception_info:
            RotationAxis(0, (1, 0, 0))

        expected_error = "`n` must be positive (n=0)"
        assert expected_error in str(exception_info)

        # n < 0
        with pytest.raises(ValueError) as exception_info:
            RotationAxis(-5, (1, 0, 0))

        expected_error = "`n` must be positive (n=-5)"
        assert expected_error in str(exception_info)

        # --- direction

        # direction is not a 3-tuple
        with pytest.raises(ValueError) as exception_info:
            RotationAxis(5, (1, 2, 3, 4))

        expected_error = "`direction` must be a 3-tuple (direction=(1, 2, 3, 4))"
        assert expected_error in str(exception_info)

        # some element of direction is not convertible to a Fraction object
        invalid_direction = (1, "b", 3)
        with pytest.raises(ValueError) as exception_info:
            RotationAxis(5, invalid_direction)

        expected_error = (
            "all elements of `direction` must convertible to Fraction objects "
            f"(direction={invalid_direction}). "
            "[caused by"
        )

        assert expected_error in str(exception_info)

        # --- location

        # location is not a 3-tuple
        with pytest.raises(ValueError) as exception_info:
            RotationAxis(5, (1, 2, 3), location=(1, 2))

        expected_error = "`location` must be a 3-tuple (location=(1, 2))"
        assert expected_error in str(exception_info)

        # some element of location is not convertible to a Fraction object
        invalid_location = (1, "b", 3)
        with pytest.raises(ValueError) as exception_info:
            RotationAxis(5, (1, 2, 3), location=invalid_location)

        expected_error = (
            "all elements of `location` must convertible to Fraction objects "
            f"(location={invalid_location}). "
            "[caused by"
        )

        assert expected_error in str(exception_info)

    @staticmethod
    def test_frozen():
        """
        Test that RotationAxis objects are immutable.
        """
        # --- Preparations

        rotation_axis = RotationAxis(3, (1, 0, 0))

        # --- Tests

        # attempt to change `n` field
        with pytest.raises(FrozenInstanceError) as exception_info:
            rotation_axis.n = 10

        expected_error = "cannot assign to field 'n'"
        assert expected_error in str(exception_info)

        # attempt to change `direction` field
        with pytest.raises(FrozenInstanceError) as exception_info:
            rotation_axis.direction = (0, 0, 1)

        expected_error = "cannot assign to field 'direction'"
        assert expected_error in str(exception_info)

        # attempt to change `location` field
        with pytest.raises(FrozenInstanceError) as exception_info:
            rotation_axis.location = (0, 0, 1)

        expected_error = "cannot assign to field 'location'"
        assert expected_error in str(exception_info)

    def test_to_julia(self):
        """
        Test `to_julia()`.
        """
        # --- Preparations

        n = 6
        direction = (1, 0, 0)

        # --- Tests

        # default location
        rotation_axis = RotationAxis(n, direction)
        rotation_axis_jl = rotation_axis.to_julia()

        assert self.jl.isa(rotation_axis_jl, self.jl.RotationAxis)
        assert rotation_axis_jl.n == n
        assert rotation_axis_jl.direction == direction
        assert rotation_axis_jl.location == (0, 0, 0)

        # non-default location
        location = (0, 0, 1)
        rotation_axis = RotationAxis(n, direction, location=location)
        rotation_axis_jl = rotation_axis.to_julia()

        assert self.jl.isa(rotation_axis_jl, self.jl.RotationAxis)
        assert rotation_axis_jl.n == n
        assert rotation_axis_jl.direction == direction
        assert rotation_axis_jl.location == location

    def test_from_julia_invalid_args(self):
        """
        Test argument checks for `from_julia()`
        """
        rotation_axis_jl = 10
        with pytest.raises(ValueError) as exception_info:
            RotationAxis.from_julia(rotation_axis_jl)

        expected_error = (
            "`rotation_axis_jl` must be a Julia `RotationAxis` object. "
            f"(rotation_axis_jl={rotation_axis_jl})."
        )
        assert expected_error in str(exception_info)

    def test_repr(self):
        """
        Test `__repr__()`.
        """
        # default location
        n = 6
        direction = (1, 0, 0)
        rotation_axis = RotationAxis(n, direction)

        assert str(rotation_axis) == (
            "RotationAxis(n=6,"
            "direction=(Fraction(1, 1), Fraction(0, 1), Fraction(0, 1)),"
            "location=(Fraction(0, 1), Fraction(0, 1), Fraction(0, 1)))"
        )

        # non-default location
        n = 6
        direction = (1, 0, 0)
        location = (1, 0, 1)
        rotation_axis = RotationAxis(n, direction, location=location)

        assert str(rotation_axis) == (
            "RotationAxis(n=6,"
            "direction=(Fraction(1, 1), Fraction(0, 1), Fraction(0, 1)),"
            "location=(Fraction(1, 1), Fraction(0, 1), Fraction(1, 1)))"
        )

    def test_eq(self):
        """
        Test `__eq__()`.
        """
        # --- Identical rotation axes

        symmetry_element_1 = RotationAxis(2, (1, 0, 0), (0, 0, 0))
        symmetry_element_2 = RotationAxis(2, (1, 0, 0), (0, 0, 0))

        assert symmetry_element_1 == symmetry_element_2

        # --- Equivalent rotation axes

        # directions differ, locations same
        symmetry_element_1 = RotationAxis(2, (1, 0, 0), (0, 0, 0))
        symmetry_element_2 = RotationAxis(2, (0.5, 0, 0), (0, 0, 0))

        assert symmetry_element_1 == symmetry_element_2

        # directions same, locations differ
        symmetry_element_1 = RotationAxis(2, (1, 0, 0), (0, 0, 0))
        symmetry_element_2 = RotationAxis(2, (1, 0, 0), (0.75, 0, 0))

        assert symmetry_element_1 == symmetry_element_2

        # directions differ, locations differ
        symmetry_element_1 = RotationAxis(2, (Fraction(2, 3), 0, 0), (0, 0, 0))
        symmetry_element_2 = RotationAxis(2, (1, 0, 0), (Fraction(3, 2), 0, 0))

        assert symmetry_element_1 == symmetry_element_2

        # --- Inequivalent rotation axes

        # orders differ
        symmetry_element_1 = RotationAxis(2, (2, 0, 0), (2, 0, 0))
        symmetry_element_2 = RotationAxis(3, (1, 0, 0), (0, 0, 0))

        assert symmetry_element_1 != symmetry_element_2

        # directions differ
        symmetry_element_1 = RotationAxis(2, (1, 1, 0), (0, 0, 0))
        symmetry_element_2 = RotationAxis(2, (1, 0, 0), (0, 0, 0))

        assert symmetry_element_1 != symmetry_element_2

        # line between locations != direction
        symmetry_element_1 = RotationAxis(2, (Fraction(1, 3), 0, 0), (0, 1, 0))
        symmetry_element_2 = RotationAxis(2, (1, 0, 0), (0, 0, 0))

        assert symmetry_element_1 != symmetry_element_2
