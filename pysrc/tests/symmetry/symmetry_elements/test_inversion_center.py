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
Unit tests for `InversionCenter` class
"""
# --- Imports

# Standard library
from dataclasses import FrozenInstanceError
from fractions import Fraction
import unittest

# External packages
import juliacall
import pytest
from xtallography.symmetry import InversionCenter, MirrorPlane

# Local packages/modules


# --- Test Suites


class test_InversionCenter(unittest.TestCase):
    """
    Test suite for the `InversionCenter` class
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
        inversion_center = InversionCenter()

        assert inversion_center.center == (0, 0, 0)
        assert all([isinstance(x, Fraction)] for x in inversion_center.center)

        # --- keyword arguments

        # non-default center
        center = (1, 0, 0)
        inversion_center = InversionCenter(center=center)

        assert inversion_center.center == center
        assert all([isinstance(x, Fraction)] for x in inversion_center.center)

    @staticmethod
    def test_init_invalid_arguments():
        """
        Test argument checks for `__init__()`
        """
        # --- center

        # center is not a 3-tuple
        with pytest.raises(ValueError) as exception_info:
            InversionCenter((1, 2, 3, 4))

        expected_error = "`center` must be a 3-tuple (center=(1, 2, 3, 4))"
        assert expected_error in str(exception_info)

        # some element of center is not convertible to a Fraction object
        invalid_center = (1, "b", 3)
        with pytest.raises(ValueError) as exception_info:
            InversionCenter(invalid_center)

        expected_error = (
            "all elements of `center` must convertible to Fraction objects "
            f"(center={invalid_center}). "
            "[caused by"
        )

        assert expected_error in str(exception_info)

    @staticmethod
    def test_frozen():
        """
        Test that InversionCenter objects are immutable.
        """
        # --- Preparations

        inversion_center = InversionCenter()

        # --- Tests

        # attempt to change `center` field
        with pytest.raises(FrozenInstanceError) as exception_info:
            inversion_center.center = (0, 0, 1)

        expected_error = "cannot assign to field 'center'"
        assert expected_error in str(exception_info)

    def test_to_julia(self):
        """
        Test `to_julia()`.
        """
        # --- default center
        inversion_center = InversionCenter()
        inversion_center_jl = inversion_center.to_julia()

        assert self.jl.isa(inversion_center_jl, self.jl.InversionCenter)
        assert inversion_center_jl.center == (0, 0, 0)

        # non-default center
        center = (0, 0, 1)
        inversion_center = InversionCenter(center=center)
        inversion_center_jl = inversion_center.to_julia()

        assert self.jl.isa(inversion_center_jl, self.jl.InversionCenter)
        assert inversion_center_jl.center == center

    def test_from_julia(self):
        """
        Test `from_julia()`.
        """
        center = (1, 2, 3)
        inversion_center_jl = self.jl.InversionCenter(center)
        inversion_center = InversionCenter.from_julia(inversion_center_jl)

        assert inversion_center == InversionCenter(center)

    def test_from_julia_invalid_args(self):
        """
        Test argument checks for `from_julia()`
        """
        inversion_center_jl = 10
        with pytest.raises(ValueError) as exception_info:
            InversionCenter.from_julia(inversion_center_jl)

        expected_error = (
            "`inversion_center_jl` must be a Julia `InversionCenter` object. "
            f"(inversion_center_jl={inversion_center_jl})."
        )
        assert expected_error in str(exception_info)

    def test_repr(self):
        """
        Test `__repr__()`.
        """
        # default center
        inversion_center = InversionCenter()

        assert str(inversion_center) == (
            "InversionCenter(center=(Fraction(0, 1), Fraction(0, 1), Fraction(0, 1)))"
        )

        # non-default center
        center = (1, 0, 1)
        inversion_center = InversionCenter(center=center)

        assert str(inversion_center) == (
            "InversionCenter(center=(Fraction(1, 1), Fraction(0, 1), Fraction(1, 1)))"
        )

    def test_eq(self):
        """
        Test `__eq__()`.
        """
        # --- Comparison with non-InversionCenter

        symmetry_element_1 = InversionCenter((1, 0, 0))
        symmetry_element_2 = MirrorPlane((0, 0, 1), (0, 0, 0))

        assert symmetry_element_1 != symmetry_element_2

        # --- Identical inversion centers

        symmetry_element_1 = InversionCenter((1, 0, 0))
        symmetry_element_2 = InversionCenter((1, 0, 0))

        assert symmetry_element_1 == symmetry_element_2

        # --- Inequivalent inversion centers

        # centers differ
        symmetry_element_1 = InversionCenter((1, 1, 0))
        symmetry_element_2 = InversionCenter((1, 0, 0))

        assert symmetry_element_1 != symmetry_element_2
