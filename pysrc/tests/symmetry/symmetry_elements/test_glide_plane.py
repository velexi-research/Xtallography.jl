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
Unit tests for `GlidePlane` class
"""
# --- Imports

# Standard library
from dataclasses import FrozenInstanceError
from fractions import Fraction
import unittest

# External packages
import juliacall
import pytest
from xtallography.symmetry import GlidePlane

# Local packages/modules


# --- Test Suites


class test_GlidePlane(unittest.TestCase):
    """
    Test suite for the `GlidePlane` class
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
        glide = (0, 0, 1)
        normal = (1, 0, 0)
        glide_plane = GlidePlane(glide, normal)

        assert glide_plane.glide == glide
        assert all([isinstance(x, Fraction)] for x in glide_plane.glide)
        assert glide_plane.normal == normal
        assert all([isinstance(x, Fraction)] for x in glide_plane.normal)
        assert glide_plane.location == (0, 0, 0)
        assert all([isinstance(x, Fraction)] for x in glide_plane.location)

        # --- keyword arguments

        # non-default location
        glide = (0, 0, 1)
        normal = (1, 0, 0)
        location = (0, 1, 0)
        glide_plane = GlidePlane(glide, normal, location=location)

        assert glide_plane.glide == glide
        assert all([isinstance(x, Fraction)] for x in glide_plane.glide)
        assert glide_plane.normal == normal
        assert all([isinstance(x, Fraction)] for x in glide_plane.normal)
        assert glide_plane.location == location
        assert all([isinstance(x, Fraction)] for x in glide_plane.location)

    @staticmethod
    def test_init_invalid_arguments():
        """
        Test argument checks for `__init__()`
        """
        # --- glide

        # glide = (0,0,0)
        with pytest.raises(ValueError) as exception_info:
            GlidePlane((0, 0, 0), (1, 0, 0))

        expected_error = "`glide` must be a nonzero vector (glide=(0, 0, 0))"
        assert expected_error in str(exception_info)

        # glide is not a 3-tuple
        with pytest.raises(ValueError) as exception_info:
            GlidePlane((1, 2, 3, 4), (1, 0, 0))

        expected_error = "`glide` must be a 3-tuple (glide=(1, 2, 3, 4))"
        assert expected_error in str(exception_info)

        # some element of glide is not convertible to a Fraction object
        invalid_glide = (1, "b", 3)
        with pytest.raises(ValueError) as exception_info:
            GlidePlane(invalid_glide, (1, 0, 0))

        expected_error = (
            "all elements of `glide` must convertible to Fraction objects "
            f"(glide={invalid_glide}). "
            "[caused by"
        )

        assert expected_error in str(exception_info)

        # --- normal

        # normal = (0,0,0)
        with pytest.raises(ValueError) as exception_info:
            GlidePlane((1, 0, 0), (0, 0, 0))

        expected_error = "`normal` must be a nonzero vector (normal=(0, 0, 0))"
        assert expected_error in str(exception_info)

        # normal is not a 3-tuple
        with pytest.raises(ValueError) as exception_info:
            GlidePlane((1, 0, 0), (1, 2, 3, 4))

        expected_error = "`normal` must be a 3-tuple (normal=(1, 2, 3, 4))"
        assert expected_error in str(exception_info)

        # some element of normal is not convertible to a Fraction object
        invalid_normal = (1, "b", 3)
        with pytest.raises(ValueError) as exception_info:
            GlidePlane((1, 0, 0), invalid_normal)

        expected_error = (
            "all elements of `normal` must convertible to Fraction objects "
            f"(normal={invalid_normal}). "
            "[caused by"
        )

        assert expected_error in str(exception_info)

        # --- location

        # location is not a 3-tuple
        with pytest.raises(ValueError) as exception_info:
            GlidePlane((1, 0, 0), (1, 2, 3), location=(1, 2))

        expected_error = "`location` must be a 3-tuple (location=(1, 2))"
        assert expected_error in str(exception_info)

        # some element of location is not convertible to a Fraction object
        invalid_location = (1, "b", 3)
        with pytest.raises(ValueError) as exception_info:
            GlidePlane((1, 0, 0), (1, 2, 3), location=invalid_location)

        expected_error = (
            "all elements of `location` must convertible to Fraction objects "
            f"(location={invalid_location}). "
            "[caused by"
        )

        assert expected_error in str(exception_info)

        # --- glide and normal are not orthogonal

        glide = (1, 0, 0)
        normal = (1, 1, 0)
        with pytest.raises(ValueError) as exception_info:
            GlidePlane(glide, normal)

        expected_error = (
            "`glide` must be orthogonal to `normal` "
            "(glide=(Fraction(1, 1), Fraction(0, 1), Fraction(0, 1)),"
            "normal=(Fraction(1, 1), Fraction(1, 1), Fraction(0, 1)))"
        )

        assert expected_error in str(exception_info)

    @staticmethod
    def test_frozen():
        """
        Test that GlidePlane objects are immutable.
        """
        # --- Preparations

        glide_plane = GlidePlane((1, 0, 0), (0, 1, 1))

        # --- Tests

        # attempt to change `glide` field
        with pytest.raises(FrozenInstanceError) as exception_info:
            glide_plane.glide = (0, 1, 1)

        expected_error = "cannot assign to field 'glide'"
        assert expected_error in str(exception_info)

        # attempt to change `normal` field
        with pytest.raises(FrozenInstanceError) as exception_info:
            glide_plane.normal = (0, 0, 1)

        expected_error = "cannot assign to field 'normal'"
        assert expected_error in str(exception_info)

        # attempt to change `location` field
        with pytest.raises(FrozenInstanceError) as exception_info:
            glide_plane.location = (1, 1, 1)

        expected_error = "cannot assign to field 'location'"
        assert expected_error in str(exception_info)

    def test_to_julia(self):
        """
        Test `to_julia()`.
        """
        # --- Preparations

        glide = (1, 0, 0)
        normal = (0, 0, 1)

        # --- Tests

        # default location
        glide_plane = GlidePlane(glide, normal)
        glide_plane_jl = glide_plane.to_julia()

        assert self.jl.isa(glide_plane_jl, self.jl.GlidePlane)
        assert glide_plane_jl.glide == glide
        assert glide_plane_jl.normal == normal
        assert glide_plane_jl.location == (0, 0, 0)

        # non-default location
        location = (0, 0, 1)
        glide_plane = GlidePlane(glide, normal, location=location)
        glide_plane_jl = glide_plane.to_julia()

        assert self.jl.isa(glide_plane_jl, self.jl.GlidePlane)
        assert glide_plane_jl.glide == glide
        assert glide_plane_jl.normal == normal
        assert glide_plane_jl.location == location

    def test_from_julia(self):
        """
        Test `from_julia()`.
        """
        # --- default location

        glide = (1, 0, 0)
        normal = (0, 0, 1)
        glide_plane_jl = self.jl.GlidePlane(glide, normal)
        glide_plane = GlidePlane.from_julia(glide_plane_jl)

        assert glide_plane == GlidePlane(glide, normal)

        # --- non-default location

        glide = (1, 0, 0)
        normal = (0, 0, 1)
        location = (0, 0, 1)
        glide_plane_jl = self.jl.GlidePlane(glide, normal, location=location)
        glide_plane = GlidePlane.from_julia(glide_plane_jl)

        assert glide_plane == GlidePlane(glide, normal, location=location)

    def test_from_julia_invalid_args(self):
        """
        Test argument checks for `from_julia()`
        """
        glide_plane_jl = 10
        with pytest.raises(ValueError) as exception_info:
            GlidePlane.from_julia(glide_plane_jl)

        expected_error = (
            "`glide_plane_jl` must be a Julia `GlidePlane` object. "
            f"(glide_plane_jl={glide_plane_jl})."
        )
        assert expected_error in str(exception_info)

    def test_repr(self):
        """
        Test `__repr__()`.
        """
        # default location
        glide = (1, 0, 0)
        normal = (0, 0, 1)
        glide_plane = GlidePlane(glide, normal)

        assert str(glide_plane) == (
            "GlidePlane(glide=(Fraction(1, 1), Fraction(0, 1), Fraction(0, 1)),"
            "normal=(Fraction(0, 1), Fraction(0, 1), Fraction(1, 1)),"
            "location=(Fraction(0, 1), Fraction(0, 1), Fraction(0, 1)))"
        )

        # non-default location
        glide = (1, 0, 0)
        normal = (0, 0, 1)
        location = (1, 0, 1)
        glide_plane = GlidePlane(glide, normal, location=location)

        assert str(glide_plane) == (
            "GlidePlane(glide=(Fraction(1, 1), Fraction(0, 1), Fraction(0, 1)),"
            "normal=(Fraction(0, 1), Fraction(0, 1), Fraction(1, 1)),"
            "location=(Fraction(1, 1), Fraction(0, 1), Fraction(1, 1)))"
        )

    def test_eq(self):
        """
        Test `__eq__()`.
        """
        # --- Identical mirror planes

        symmetry_element_1 = GlidePlane((1, 0, 0), (0, 0, 1), (0, 0, 0))
        symmetry_element_2 = GlidePlane((1, 0, 0), (0, 0, 1), (0, 0, 0))

        assert symmetry_element_1 == symmetry_element_2

        # --- Equivalent mirror planes

        # glides same, normals differ, locations same
        symmetry_element_1 = GlidePlane((1, 0, 0), (0, 0, 1), (0, 0, 0))
        symmetry_element_2 = GlidePlane((1, 0, 0), (0, 0, 0.5), (0, 0, 0))

        assert symmetry_element_1 == symmetry_element_2

        # glides same, normals same, locations differ
        symmetry_element_1 = GlidePlane((1, 0, 0), (0, 0, 1), (0, 0, 0))
        symmetry_element_2 = GlidePlane((1, 0, 0), (0, 0, 1), (0, 0.75, 0))

        assert symmetry_element_1 == symmetry_element_2

        # glides same, normals differ, locations differ
        symmetry_element_1 = GlidePlane((1, 0, 0), (0, 0, Fraction(2, 3)), (0, 0, 0))
        symmetry_element_2 = GlidePlane((1, 0, 0), (0, 0, 1), (0, Fraction(3, 2), 0))

        assert symmetry_element_1 == symmetry_element_2

        # --- Inequivalent mirror planes

        # glides differ
        symmetry_element_1 = GlidePlane((1, 1, 0), (0, 0, 1), (0, 0, 0))
        symmetry_element_2 = GlidePlane((1, 0, 0), (0, 0, 1), (0, 0, 0))

        assert symmetry_element_1 != symmetry_element_2

        # glides differ in scale
        symmetry_element_1 = GlidePlane((1, 0, 0), (0, 0, 1), (0, 0, 0))
        symmetry_element_2 = GlidePlane((0.5, 0, 0), (0, 0, 1), (0, 0, 0))

        assert symmetry_element_1 != symmetry_element_2

        # glides differ, normals same, locations differ
        symmetry_element_1 = GlidePlane((0.75, 0, 0), (0, 0, 1), (1, 1, 0))
        symmetry_element_2 = GlidePlane((1, 0, 0), (0, 0, 1), (0, 0, 0))

        assert symmetry_element_1 != symmetry_element_2

        # normals differ
        symmetry_element_1 = GlidePlane((1, 1, 0), (0, 0, 1))
        symmetry_element_2 = GlidePlane((1, 0, 0), (0, 0, 1))

        assert symmetry_element_1 != symmetry_element_2

        # line between locations is not orthogonal to normal
        symmetry_element_1 = GlidePlane((0, 0, 1), (Fraction(1, 3), 0, 0), (1, 1, 0))
        symmetry_element_2 = GlidePlane((0, 0, 1), (1, 0, 0), (0, 0, 0))

        assert symmetry_element_1 != symmetry_element_2
