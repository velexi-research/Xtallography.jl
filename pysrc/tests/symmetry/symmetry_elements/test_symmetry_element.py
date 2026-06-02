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
Unit tests for `SymmetryElement` class
"""
# --- Imports

# Standard library
import unittest

# External packages
import juliacall
import pytest
from xtallography.symmetry import SymmetryElement
from xtallography.symmetry import (
    GlidePlane,
    InversionCenter,
    MirrorPlane,
    RotationAxis,
    RotoinversionAxis,
    ScrewAxis,
)

# Local packages/modules


# --- Test Suites


class test_SymmetryElement(unittest.TestCase):
    """
    Test suite for the `SymmetryElement` class
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

    def test_from_julia(self):
        """
        Test `SymmetryElement.from_julia()`.

        Note: only tests that the correctness of the class of the returned object.
        Correctness of the returned object is tested by the from_julia() tests for
        concrete subclasses of SymmetryElement.
        """
        # --- GlidePlane

        glide = (1, 2, 0)
        normal = (0, 0, 3)
        symmetry_element_jl = self.jl.GlidePlane(glide, normal)
        symmetry_element = SymmetryElement.from_julia(symmetry_element_jl)

        assert isinstance(symmetry_element, GlidePlane)

        # --- InversionCenter

        center = (1, 2, 0)
        symmetry_element_jl = self.jl.InversionCenter(center)
        symmetry_element = SymmetryElement.from_julia(symmetry_element_jl)

        assert isinstance(symmetry_element, InversionCenter)

        # --- MirrorPlane

        normal = (1, 1, 1)
        symmetry_element_jl = self.jl.MirrorPlane(normal)
        symmetry_element = SymmetryElement.from_julia(symmetry_element_jl)

        assert isinstance(symmetry_element, MirrorPlane)

        # --- RotationAxis

        n = 4
        direction = (1, 2, 0)
        symmetry_element_jl = self.jl.RotationAxis(n, direction)
        symmetry_element = SymmetryElement.from_julia(symmetry_element_jl)

        assert isinstance(symmetry_element, RotationAxis)

        # --- RotoinversionAxis

        n = 4
        direction = (1, 2, 0)
        center = (1, 2, 0)
        symmetry_element_jl = self.jl.RotoinversionAxis(n, direction, center=center)
        symmetry_element = SymmetryElement.from_julia(symmetry_element_jl)

        assert isinstance(symmetry_element, RotoinversionAxis)

        # --- ScrewAxis

        n = 6
        m = 4
        direction = (1, 2, 3)
        symmetry_element_jl = self.jl.ScrewAxis(n, m, direction)
        symmetry_element = SymmetryElement.from_julia(symmetry_element_jl)

        assert isinstance(symmetry_element, ScrewAxis)

    def test_from_julia_invalid_args(self):
        """
        Test argument checks for `from_julia()`
        """
        symmetry_element_jl = 10
        with pytest.raises(ValueError) as exception_info:
            SymmetryElement.from_julia(symmetry_element_jl)

        expected_error = (
            "`symmetry_element_jl` must be a Julia `SymmetryElement` object. "
            f"(symmetry_element_jl={symmetry_element_jl})."
        )
        assert expected_error in str(exception_info)
