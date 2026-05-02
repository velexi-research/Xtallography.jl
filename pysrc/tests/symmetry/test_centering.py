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
Unit tests for `xtallography.symmetry.centering` module
"""
# --- Imports

# Standard library
import unittest

# External packages
import juliacall
import pytest
from xtallography.symmetry import Centering

# Local packages/modules


# --- Test Suites


class test_xtallography_symmetry_centering_Centering(unittest.TestCase):
    """
    Test suite for the `Centering` class
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

    def test_to_julia(self):
        """
        Test `to_julia()` method.
        """
        # --- Tests

        # primitive_centering
        centering = Centering.PRIMITIVE
        centering_jl = centering.to_julia()
        assert self.jl.isa(centering_jl, self.jl.PrimitiveCentering)

        # base_centering
        centering = Centering.BASE
        centering_jl = centering.to_julia()
        assert self.jl.isa(centering_jl, self.jl.BaseCentering)

        # body_centering
        centering = Centering.BODY
        centering_jl = centering.to_julia()
        assert self.jl.isa(centering_jl, self.jl.BodyCentering)

        # face_centering
        centering = Centering.FACE
        centering_jl = centering.to_julia()
        assert self.jl.isa(centering_jl, self.jl.FaceCentering)

    def test_from_julia(self):
        """
        Test `from_julia()` method.
        """
        # --- Tests

        # primitive_centering
        centering_jl = self.jl.primitive_centering
        centering = Centering.from_julia(centering_jl)
        assert centering == Centering.PRIMITIVE

        # base_centering
        centering_jl = self.jl.base_centering
        centering = Centering.from_julia(centering_jl)
        assert centering == Centering.BASE

        # body_centering
        centering_jl = self.jl.body_centering
        centering = Centering.from_julia(centering_jl)
        assert centering == Centering.BODY

        # face_centering
        centering_jl = self.jl.face_centering
        centering = Centering.from_julia(centering_jl)
        assert centering == Centering.FACE

    def test_from_julia_invalid_args(self):
        """
        Test argument checks for `from_julia()`.
        """
        # --- Tests

        # ------ `centering_jl` is not a Julia Centering object

        centering_jl_invalid = "not a Julia Centering object"
        with pytest.raises(ValueError) as exception_info:
            Centering.from_julia(centering_jl_invalid)

        expected_error = (
            "`centering_jl` must be a Julia `Centering` object. "
            f"(centering_jl={centering_jl_invalid})."
        )
        assert expected_error in str(exception_info)

        # ------ `centering_jl` is an unsupported Julia Centering object

        self.jl.seval("struct UnsupportedCentering <: Centering end")
        centering_jl_invalid = self.jl.UnsupportedCentering()
        with pytest.raises(ValueError) as exception_info:
            Centering.from_julia(centering_jl_invalid)

        expected_error = (
            "Unsupported Centering type. (centering_jl=UnsupportedCentering)"
        )
        assert expected_error in str(exception_info)

    @staticmethod
    def test_values():
        """
        Test values() method.
        """
        # --- Tests

        # Get full list of values
        values = set(Centering.values())
        assert values == {
            "primitive",
            "base",
            "body",
            "face",
        }
