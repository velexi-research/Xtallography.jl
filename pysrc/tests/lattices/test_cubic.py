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
Unit tests for `xtallography.lattices.cubic` module
"""
# --- Imports

# Standard library
import unittest

# External packages
import pytest

# Local packages/modules
from xtallography import _JL
from xtallography.lattices import LatticeSystem, Centering, CubicUnitCell


# --- Test Suites


class test_xtallography_lattice_cubic(unittest.TestCase):
    """
    Test suite for the `xtallography.lattice.cubic` module
    """

    # --- Fixtures

    def setUp(self):
        """
        Prepare for test.
        """

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
        # --- Preparations

        a = 1

        # --- Tests

        # default centering
        unit_cell = CubicUnitCell(a)

        assert unit_cell.a == a
        assert isinstance(unit_cell.lattice_system, LatticeSystem)
        assert unit_cell.lattice_system == LatticeSystem.CUBIC
        assert isinstance(unit_cell.centering, Centering)
        assert unit_cell.centering == Centering.PRIMITIVE

        # centering = primitive
        unit_cell = CubicUnitCell(a, centering=Centering.PRIMITIVE)

        assert unit_cell.a == a
        assert unit_cell.lattice_system == LatticeSystem.CUBIC
        assert unit_cell.centering == Centering.PRIMITIVE

        # centering = body-centered
        unit_cell = CubicUnitCell(a, centering=Centering.BODY_CENTERED)

        assert unit_cell.a == a
        assert unit_cell.lattice_system == LatticeSystem.CUBIC
        assert unit_cell.centering == Centering.BODY_CENTERED

        # centering = face-centered
        unit_cell = CubicUnitCell(a, centering=Centering.FACE_CENTERED)

        assert unit_cell.a == a
        assert unit_cell.lattice_system == LatticeSystem.CUBIC
        assert unit_cell.centering == Centering.FACE_CENTERED

        # centering argument is a str
        unit_cell = CubicUnitCell(a, centering="primitive")

        assert unit_cell.a == a
        assert isinstance(unit_cell.lattice_system, LatticeSystem)
        assert unit_cell.lattice_system == LatticeSystem.CUBIC
        assert isinstance(unit_cell.centering, Centering)
        assert unit_cell.centering == Centering.PRIMITIVE

    @staticmethod
    def test_init_invalid_args():
        """
        Test `__init__()`: invalid arguments.
        """
        # --- Tests

        # ------ Invalid `a`

        # a < 0
        invalid_value = -1
        with pytest.raises(ValueError) as exception_info:
            CubicUnitCell(invalid_value)

        expected_error = f"`a` must be positive. (a={invalid_value})"
        assert expected_error in str(exception_info)

        # a = 0
        invalid_value = 0
        with pytest.raises(ValueError) as exception_info:
            CubicUnitCell(invalid_value)

        expected_error = f"`a` must be positive. (a={invalid_value})"
        assert expected_error in str(exception_info)

        # ------ Invalid `centering`

        # base-centered
        invalid_value = "base-centered"
        with pytest.raises(ValueError) as exception_info:
            CubicUnitCell(1.0, centering=invalid_value)

        expected_error = "('cubic', 'base-centered') is not a valid Bravais lattice."
        assert expected_error in str(exception_info)

    @staticmethod
    def test_to_julia():
        """
        Test `to_julia()`.
        """
        # --- Preparations

        a = 1

        unit_cell = CubicUnitCell(a)

        # --- Tests

        unit_cell_jl = unit_cell.to_julia()
        assert _JL.isa(unit_cell_jl, _JL.CubicLatticeConstants)
