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
Unit tests for `RotoinversionAxis` class
"""
# --- Imports

# Standard library
from dataclasses import FrozenInstanceError
from fractions import Fraction
import unittest

# External packages
import juliacall
import pytest
from xtallography.symmetry import RotoinversionAxis

# Local packages/modules


# --- Test Suites


class test_RotoinversionAxis(unittest.TestCase):
    """
    Test suite for the `RotoinversionAxis` class
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

        # default center
        n = 6
        direction = (1, 0, 0)
        rotoinversion_axis = RotoinversionAxis(n, direction)

        assert rotoinversion_axis.n == n
        assert rotoinversion_axis.direction == direction
        assert all([isinstance(x, Fraction)] for x in rotoinversion_axis.direction)
        assert rotoinversion_axis.center == (0, 0, 0)
        assert all([isinstance(x, Fraction)] for x in rotoinversion_axis.center)

        # --- keyword arguments

        # non-default center
        n = 6
        direction = (1, 0, 0)
        center = (0, 1, 0)
        rotoinversion_axis = RotoinversionAxis(n, direction, center=center)

        assert rotoinversion_axis.n == n
        assert rotoinversion_axis.direction == direction
        assert all([isinstance(x, Fraction)] for x in rotoinversion_axis.direction)
        assert rotoinversion_axis.center == center
        assert all([isinstance(x, Fraction)] for x in rotoinversion_axis.center)

    @staticmethod
    def test_init_invalid_arguments():
        """
        Test argument checks for `__init__()`
        """
        # --- n

        # n = 0
        with pytest.raises(ValueError) as exception_info:
            RotoinversionAxis(0, (1, 0, 0))

        expected_error = "`n` must be positive (n=0)"
        assert expected_error in str(exception_info)

        # n < 0
        with pytest.raises(ValueError) as exception_info:
            RotoinversionAxis(-5, (1, 0, 0))

        expected_error = "`n` must be positive (n=-5)"
        assert expected_error in str(exception_info)

        # --- direction

        # direction is not a 3-tuple
        with pytest.raises(ValueError) as exception_info:
            RotoinversionAxis(5, (1, 2, 3, 4))

        expected_error = "`direction` must be a 3-tuple (direction=(1, 2, 3, 4))"
        assert expected_error in str(exception_info)

        # some element of direction is not convertible to a Fraction object
        invalid_direction = (1, "b", 3)
        with pytest.raises(ValueError) as exception_info:
            RotoinversionAxis(5, invalid_direction)

        expected_error = (
            "all elements of `direction` must convertible to Fraction objects "
            f"(direction={invalid_direction}). "
            "[caused by"
        )

        assert expected_error in str(exception_info)

        # --- center

        # center is not a 3-tuple
        with pytest.raises(ValueError) as exception_info:
            RotoinversionAxis(5, (1, 2, 3), center=(1, 2))

        expected_error = "`center` must be a 3-tuple (center=(1, 2))"
        assert expected_error in str(exception_info)

        # some element of center is not convertible to a Fraction object
        invalid_center = (1, "b", 3)
        with pytest.raises(ValueError) as exception_info:
            RotoinversionAxis(5, (1, 2, 3), invalid_center)

        expected_error = (
            "all elements of `center` must convertible to Fraction objects "
            f"(center={invalid_center}). "
            "[caused by"
        )

        assert expected_error in str(exception_info)

    @staticmethod
    def test_frozen():
        """
        Test that RotoinversionAxis objects are immutable.
        """
        # --- Preparations

        rotoinversion_axis = RotoinversionAxis(3, (1, 0, 0))

        # --- Tests

        # attempt to change `n` field
        with pytest.raises(FrozenInstanceError) as exception_info:
            rotoinversion_axis.n = 10

        expected_error = "cannot assign to field 'n'"
        assert expected_error in str(exception_info)

        # attempt to change `direction` field
        with pytest.raises(FrozenInstanceError) as exception_info:
            rotoinversion_axis.direction = (0, 0, 1)

        expected_error = "cannot assign to field 'direction'"
        assert expected_error in str(exception_info)

        # attempt to change `center` field
        with pytest.raises(FrozenInstanceError) as exception_info:
            rotoinversion_axis.center = (0, 0, 1)

        expected_error = "cannot assign to field 'center'"
        assert expected_error in str(exception_info)

    def test_to_julia(self):
        """
        Test `to_julia()`.
        """
        # --- Preparations

        n = 6
        direction = (1, 0, 0)

        # --- Tests

        # default center
        rotoinversion_axis = RotoinversionAxis(n, direction)
        rotoinversion_axis_jl = rotoinversion_axis.to_julia()

        assert self.jl.isa(rotoinversion_axis_jl, self.jl.RotoinversionAxis)
        assert rotoinversion_axis_jl.n == n
        assert rotoinversion_axis_jl.direction == direction
        assert rotoinversion_axis_jl.center == (0, 0, 0)

        # non-default center
        center = (0, 0, 1)
        rotoinversion_axis = RotoinversionAxis(n, direction, center=center)
        rotoinversion_axis_jl = rotoinversion_axis.to_julia()

        assert self.jl.isa(rotoinversion_axis_jl, self.jl.RotoinversionAxis)
        assert rotoinversion_axis_jl.n == n
        assert rotoinversion_axis_jl.direction == direction
        assert rotoinversion_axis_jl.center == center

    def test_from_julia(self):
        """
        Test `from_julia()`.
        """
        # --- default center

        n = 4
        direction = (1, 2, 3)
        rotoinversion_axis_jl = self.jl.RotoinversionAxis(n, direction)
        rotoinversion_axis = RotoinversionAxis.from_julia(rotoinversion_axis_jl)

        assert rotoinversion_axis == RotoinversionAxis(n, direction, center=(0, 0, 0))

        # --- non-default center

        n = 6
        direction = (1, 2, 3)
        center = (4, 5, 6)
        rotoinversion_axis_jl = self.jl.RotoinversionAxis(n, direction, center=center)
        rotoinversion_axis = RotoinversionAxis.from_julia(rotoinversion_axis_jl)

        assert rotoinversion_axis == RotoinversionAxis(n, direction, center=center)

    def test_from_julia_invalid_args(self):
        """
        Test argument checks for `from_julia()`
        """
        rotoinversion_axis_jl = 10
        with pytest.raises(ValueError) as exception_info:
            RotoinversionAxis.from_julia(rotoinversion_axis_jl)

        expected_error = (
            "`rotoinversion_axis_jl` must be a Julia `RotoinversionAxis` object. "
            f"(rotoinversion_axis_jl={rotoinversion_axis_jl})."
        )
        assert expected_error in str(exception_info)

    def test_repr(self):
        """
        Test `__repr__()`.
        """
        # default center
        n = 6
        direction = (1, 0, 0)
        rotoinversion_axis = RotoinversionAxis(n, direction)

        assert str(rotoinversion_axis) == (
            "RotoinversionAxis(n=6,"
            "direction=(Fraction(1, 1), Fraction(0, 1), Fraction(0, 1)),"
            "center=(Fraction(0, 1), Fraction(0, 1), Fraction(0, 1)))"
        )

        # non-default center
        n = 6
        direction = (1, 0, 0)
        center = (1, 0, 1)
        rotoinversion_axis = RotoinversionAxis(n, direction, center=center)

        assert str(rotoinversion_axis) == (
            "RotoinversionAxis(n=6,"
            "direction=(Fraction(1, 1), Fraction(0, 1), Fraction(0, 1)),"
            "center=(Fraction(1, 1), Fraction(0, 1), Fraction(1, 1)))"
        )

    def test_eq(self):
        """
        Test `__eq__()`.
        """
        # --- Identical rotation axes

        symmetry_element_1 = RotoinversionAxis(2, (1, 0, 0), (0, 0, 0))
        symmetry_element_2 = RotoinversionAxis(2, (1, 0, 0), (0, 0, 0))

        assert symmetry_element_1 == symmetry_element_2

        # --- Equivalent rotation axes

        # directions differ, center same
        symmetry_element_1 = RotoinversionAxis(2, (1, 0, 0), (0, 0, 0))
        symmetry_element_2 = RotoinversionAxis(2, (0.5, 0, 0), (0, 0, 0))

        assert symmetry_element_1 == symmetry_element_2

        # --- Inequivalent rotation axes

        # order is different
        symmetry_element_1 = RotoinversionAxis(2, (2, 0, 0), (2, 0, 0))
        symmetry_element_2 = RotoinversionAxis(3, (1, 0, 0), (0, 0, 0))

        assert symmetry_element_1 != symmetry_element_2

        # directions differ
        symmetry_element_1 = RotoinversionAxis(2, (1, 1, 0), (0, 0, 0))
        symmetry_element_2 = RotoinversionAxis(2, (1, 0, 0), (0, 0, 0))

        assert symmetry_element_1 != symmetry_element_2

        # centers differ
        symmetry_element_1 = RotoinversionAxis(2, (Fraction(1, 3), 0, 0), (0, 1, 0))
        symmetry_element_2 = RotoinversionAxis(2, (1, 0, 0), (0, 0, 0))

        assert symmetry_element_1 != symmetry_element_2
