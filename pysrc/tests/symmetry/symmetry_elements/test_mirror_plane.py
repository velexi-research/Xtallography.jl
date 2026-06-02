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
Unit tests for `MirrorPlane` class
"""
# --- Imports

# Standard library
from dataclasses import FrozenInstanceError
from fractions import Fraction
import unittest

# External packages
import juliacall
import pytest
from xtallography.symmetry import MirrorPlane

# Local packages/modules


# --- Test Suites


class test_MirrorPlane(unittest.TestCase):
    """
    Test suite for the `MirrorPlane` class
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
        normal = (1, 0, 0)
        mirror_plane = MirrorPlane(normal)

        assert mirror_plane.normal == normal
        assert all([isinstance(x, Fraction)] for x in mirror_plane.normal)
        assert mirror_plane.location == (0, 0, 0)
        assert all([isinstance(x, Fraction)] for x in mirror_plane.location)

        # --- keyword arguments

        # non-default location
        normal = (1, 0, 0)
        location = (0, 1, 0)
        mirror_plane = MirrorPlane(normal, location=location)

        assert mirror_plane.normal == normal
        assert all([isinstance(x, Fraction)] for x in mirror_plane.normal)
        assert mirror_plane.location == location
        assert all([isinstance(x, Fraction)] for x in mirror_plane.location)

    @staticmethod
    def test_init_invalid_arguments():
        """
        Test argument checks for `__init__()`
        """
        # --- normal

        # normal = (0,0,0)
        with pytest.raises(ValueError) as exception_info:
            MirrorPlane((0, 0, 0))

        expected_error = "`normal` must be a nonzero vector (normal=(0, 0, 0))"
        assert expected_error in str(exception_info)

        # normal is not a 3-tuple
        with pytest.raises(ValueError) as exception_info:
            MirrorPlane((1, 2, 3, 4))

        expected_error = "`normal` must be a 3-tuple (normal=(1, 2, 3, 4))"
        assert expected_error in str(exception_info)

        # some element of normal is not convertible to a Fraction object
        invalid_normal = (1, "b", 3)
        with pytest.raises(ValueError) as exception_info:
            MirrorPlane(invalid_normal)

        expected_error = (
            "all elements of `normal` must convertible to Fraction objects "
            f"(normal={invalid_normal}). "
            "[caused by"
        )

        assert expected_error in str(exception_info)

        # --- location

        # location is not a 3-tuple
        with pytest.raises(ValueError) as exception_info:
            MirrorPlane((1, 2, 3), location=(1, 2))

        expected_error = "`location` must be a 3-tuple (location=(1, 2))"
        assert expected_error in str(exception_info)

        # some element of location is not convertible to a Fraction object
        invalid_location = (1, "b", 3)
        with pytest.raises(ValueError) as exception_info:
            MirrorPlane((1, 2, 3), invalid_location)

        expected_error = (
            "all elements of `location` must convertible to Fraction objects "
            f"(location={invalid_location}). "
            "[caused by"
        )

        assert expected_error in str(exception_info)

    @staticmethod
    def test_frozen():
        """
        Test that MirrorPlane objects are immutable.
        """
        # --- Preparations

        mirror_plane = MirrorPlane((1, 0, 0))

        # --- Tests

        # attempt to change `normal` field
        with pytest.raises(FrozenInstanceError) as exception_info:
            mirror_plane.normal = (0, 0, 1)

        expected_error = "cannot assign to field 'normal'"
        assert expected_error in str(exception_info)

        # attempt to change `location` field
        with pytest.raises(FrozenInstanceError) as exception_info:
            mirror_plane.location = (0, 0, 1)

        expected_error = "cannot assign to field 'location'"
        assert expected_error in str(exception_info)

    def test_to_julia(self):
        """
        Test `to_julia()`.
        """
        # --- Preparations

        normal = (1, 0, 0)

        # --- Tests

        # default location
        mirror_plane = MirrorPlane(normal)
        mirror_plane_jl = mirror_plane.to_julia()

        assert self.jl.isa(mirror_plane_jl, self.jl.MirrorPlane)
        assert mirror_plane_jl.normal == normal
        assert mirror_plane_jl.location == (0, 0, 0)

        # non-default location
        location = (0, 0, 1)
        mirror_plane = MirrorPlane(normal, location=location)
        mirror_plane_jl = mirror_plane.to_julia()

        assert self.jl.isa(mirror_plane_jl, self.jl.MirrorPlane)
        assert mirror_plane_jl.normal == normal
        assert mirror_plane_jl.location == location

    def test_from_julia(self):
        """
        Test `from_julia()`.
        """
        # --- default location

        normal = (1, 2, 3)
        mirror_plane_jl = self.jl.MirrorPlane(normal)
        mirror_plane = MirrorPlane.from_julia(mirror_plane_jl)

        assert mirror_plane == MirrorPlane(normal)

        # --- non-default location

        normal = (1, 2, 3)
        location = (4, 5, 6)
        mirror_plane_jl = self.jl.MirrorPlane(normal, location=location)
        mirror_plane = MirrorPlane.from_julia(mirror_plane_jl)

        assert mirror_plane == MirrorPlane(normal, location=location)

    def test_from_julia_invalid_args(self):
        """
        Test argument checks for `from_julia()`
        """
        mirror_plane_jl = 10
        with pytest.raises(ValueError) as exception_info:
            MirrorPlane.from_julia(mirror_plane_jl)

        expected_error = (
            "`mirror_plane_jl` must be a Julia `MirrorPlane` object. "
            f"(mirror_plane_jl={mirror_plane_jl})."
        )
        assert expected_error in str(exception_info)

    def test_repr(self):
        """
        Test `__repr__()`.
        """
        # default location
        normal = (1, 0, 0)
        mirror_plane = MirrorPlane(normal)

        assert str(mirror_plane) == (
            "MirrorPlane(normal=(Fraction(1, 1), Fraction(0, 1), Fraction(0, 1)),"
            "location=(Fraction(0, 1), Fraction(0, 1), Fraction(0, 1)))"
        )

        # non-default location
        normal = (1, 0, 0)
        location = (1, 0, 1)
        mirror_plane = MirrorPlane(normal, location=location)

        assert str(mirror_plane) == (
            "MirrorPlane(normal=(Fraction(1, 1), Fraction(0, 1), Fraction(0, 1)),"
            "location=(Fraction(1, 1), Fraction(0, 1), Fraction(1, 1)))"
        )

    def test_eq(self):
        """
        Test `__eq__()`.
        """
        # --- Identical mirror planes

        symmetry_element_1 = MirrorPlane((1, 0, 0), (0, 0, 0))
        symmetry_element_2 = MirrorPlane((1, 0, 0), (0, 0, 0))

        assert symmetry_element_1 == symmetry_element_2

        # --- Equivalent mirror planes

        # normals differ, locations same
        symmetry_element_1 = MirrorPlane((1, 0, 0), (0, 0, 0))
        symmetry_element_2 = MirrorPlane((0.5, 0, 0), (0, 0, 0))

        assert symmetry_element_1 == symmetry_element_2

        # normals same, locations differ
        symmetry_element_1 = MirrorPlane((1, 0, 0), (0, 0, 0))
        symmetry_element_2 = MirrorPlane((1, 0, 0), (0, 0.75, 0))

        assert symmetry_element_1 == symmetry_element_2

        # normals differ, locations differ
        symmetry_element_1 = MirrorPlane((Fraction(2, 3), 0, 0), (0, 0, 0))
        symmetry_element_2 = MirrorPlane((1, 0, 0), (0, Fraction(3, 2), 0))

        assert symmetry_element_1 == symmetry_element_2

        # --- Inequivalent mirror planes

        # normals differ
        symmetry_element_1 = MirrorPlane((1, 1, 0), (0, 0, 0))
        symmetry_element_2 = MirrorPlane((1, 0, 0), (0, 0, 0))

        assert symmetry_element_1 != symmetry_element_2

        # line between locations is not orthogonal to normal
        symmetry_element_1 = MirrorPlane((Fraction(1, 3), 0, 0), (1, 1, 0))
        symmetry_element_2 = MirrorPlane((1, 0, 0), (0, 0, 0))

        assert symmetry_element_1 != symmetry_element_2
