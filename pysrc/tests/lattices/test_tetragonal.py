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
Unit tests for `xtallography.lattices.tetragonal` module
"""
# --- Imports

# Standard library
import unittest

# External packages
import pytest

# Local packages/modules
from xtallography import _JL
from xtallography.lattices import LatticeSystem, Centering, TetragonalUnitCell


# --- Test Suites


class test_xtallography_lattice_tetragonal(unittest.TestCase):
    """
    Test suite for the `xtallography.lattice.tetragonal` module
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
        c = 3

        # --- Tests

        # default centering
        unit_cell = TetragonalUnitCell(a, c)

        assert unit_cell.a == a
        assert unit_cell.c == c
        assert unit_cell.lattice_system == LatticeSystem.TETRAGONAL
        assert unit_cell.centering == Centering.PRIMITIVE

        # centering = primitive
        unit_cell = TetragonalUnitCell(a, c, centering=Centering.PRIMITIVE)
        assert unit_cell.a == a
        assert unit_cell.c == c
        assert unit_cell.lattice_system == LatticeSystem.TETRAGONAL
        assert unit_cell.centering == Centering.PRIMITIVE

        # centering = body-centered
        unit_cell = TetragonalUnitCell(a, c, centering=Centering.BODY_CENTERED)
        assert unit_cell.a == a
        assert unit_cell.c == c
        assert unit_cell.lattice_system == LatticeSystem.TETRAGONAL
        assert unit_cell.centering == Centering.BODY_CENTERED

        # centering argument is a str
        unit_cell = TetragonalUnitCell(a, c, centering="primitive")

        assert unit_cell.a == a
        assert unit_cell.c == c
        assert isinstance(unit_cell.lattice_system, LatticeSystem)
        assert unit_cell.lattice_system == LatticeSystem.TETRAGONAL
        assert isinstance(unit_cell.centering, Centering)
        assert unit_cell.centering == Centering.PRIMITIVE

    @staticmethod
    def test_init_invalid_args():
        """
        Test `__init__()`: invalid arguments.
        """
        # --- Preparations

        a = 1
        c = 3

        # --- Tests

        # ------ Invalid `a`

        # a < 0
        invalid_value = -1
        with pytest.raises(ValueError) as exception_info:
            TetragonalUnitCell(invalid_value, c)

        expected_error = f"`a` must be positive. (a={invalid_value})"
        assert expected_error in str(exception_info)

        # a = 0
        invalid_value = 0
        with pytest.raises(ValueError) as exception_info:
            TetragonalUnitCell(invalid_value, c)

        expected_error = f"`a` must be positive. (a={invalid_value})"
        assert expected_error in str(exception_info)

        # ------ Invalid `c`

        # c < 0
        invalid_value = -1
        with pytest.raises(ValueError) as exception_info:
            TetragonalUnitCell(a, invalid_value)

        expected_error = f"`c` must be positive. (c={invalid_value})"
        assert expected_error in str(exception_info)

        # c = 0
        invalid_value = 0
        with pytest.raises(ValueError) as exception_info:
            TetragonalUnitCell(a, invalid_value)

        expected_error = f"`c` must be positive. (c={invalid_value})"
        assert expected_error in str(exception_info)

        # ------ Invalid `centering`

        # base-centered
        invalid_value = "base-centered"
        with pytest.raises(ValueError) as exception_info:
            TetragonalUnitCell(a, c, centering=invalid_value)

        expected_error = (
            "('tetragonal', 'base-centered') is not a valid Bravais lattice."
        )
        assert expected_error in str(exception_info)

        # face-centered
        invalid_value = "face-centered"
        with pytest.raises(ValueError) as exception_info:
            TetragonalUnitCell(a, c, centering=invalid_value)

        expected_error = (
            "('tetragonal', 'face-centered') is not a valid Bravais lattice."
        )
        assert expected_error in str(exception_info)

    @staticmethod
    def test_to_julia():
        """
        Test `to_julia()`.
        """
        # --- Preparations

        a = 1
        c = 3

        unit_cell = TetragonalUnitCell(a, c)

        # --- Tests

        unit_cell_jl = unit_cell.to_julia()
        assert _JL.isa(unit_cell_jl, _JL.TetragonalLatticeConstants)
