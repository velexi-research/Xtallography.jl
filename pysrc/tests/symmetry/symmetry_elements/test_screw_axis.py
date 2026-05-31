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
import unittest

# External packages
import juliacall
import pytest
from xtallography.symmetry import ScrewAxis

# Local packages/modules


# --- Test Suites


class test_xtallography_symmetry_symmetry_elements_ScrewAxis(unittest.TestCase):
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
        # --- Tests

        axis = "1,0,0"
        n = 6
        m = 4
        screw_axis = ScrewAxis(axis, n, m)

        assert screw_axis.axis == axis
        assert screw_axis.n == n
        assert screw_axis.m == m

    @staticmethod
    def test_init_invalid_arguments():
        """
        Test argument checks for `__init__()`
        """
        # --- Tests

        # n = 0
        with pytest.raises(ValueError) as exception_info:
            ScrewAxis("1,0,0", 0, 4)

        expected_error = "`n` must be positive (n=0)"
        assert expected_error in str(exception_info)

        # n < 0
        with pytest.raises(ValueError) as exception_info:
            ScrewAxis("1,0,0", -5, 4)

        expected_error = "`n` must be positive (n=-5)"
        assert expected_error in str(exception_info)

        # m = 0
        with pytest.raises(ValueError) as exception_info:
            ScrewAxis("1,0,0", 6, 0)

        expected_error = "`m` must be positive (m=0)"
        assert expected_error in str(exception_info)

        # m < 0
        with pytest.raises(ValueError) as exception_info:
            ScrewAxis("1,0,0", 6, -4)

        expected_error = "`m` must be positive (m=-4)"
        assert expected_error in str(exception_info)

        # m = n
        with pytest.raises(ValueError) as exception_info:
            ScrewAxis("1,0,0", 6, 6)

        expected_error = "`m` must be less than `n` (n=6,m=6)"
        assert expected_error in str(exception_info)

        # m > n
        with pytest.raises(ValueError) as exception_info:
            ScrewAxis("1,0,0", 6, 7)

        expected_error = "`m` must be less than `n` (n=6,m=7)"
        assert expected_error in str(exception_info)

    @staticmethod
    def test_frozen():
        """
        Test that ScrewAxis objects are immutable.
        """
        # --- Preparations

        screw_axis = ScrewAxis("1,0,0", 6, 4)

        # --- Tests

        # attempt to change `axis` field
        with pytest.raises(FrozenInstanceError) as exception_info:
            screw_axis.axis = "0,1,0"

        expected_error = "cannot assign to field 'axis'"
        assert expected_error in str(exception_info)

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

    def test_to_julia(self):
        """
        Test `to_julia()`.
        """
        # --- Preparations

        axis = "1,0,0"
        n = 6
        m = 4
        screw_axis = ScrewAxis(axis, n, m)

        # --- Tests

        screw_axis_jl = screw_axis.to_julia()
        assert self.jl.isa(screw_axis_jl, self.jl.ScrewAxis)
        assert screw_axis_jl.axis == axis
        assert screw_axis_jl.n == n
        assert screw_axis_jl.m == m

    def test_from_julia(self):
        """
        Test `from_julia()`.
        """
        # --- Tests

        axis = "1,0,0"
        n = 6
        m = 4
        screw_axis_jl = self.jl.ScrewAxis(axis, n, m)

        screw_axis = ScrewAxis.from_julia(screw_axis_jl)

        assert screw_axis.axis == axis
        assert screw_axis.n == n
        assert screw_axis.m == m

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
        # --- Tests

        axis = "1,0,0"
        n = 6
        m = 4
        screw_axis = ScrewAxis(axis, n, m)

        assert str(screw_axis) == "ScrewAxis(axis='1,0,0',n=6,m=4)"
