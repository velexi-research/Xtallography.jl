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
Unit tests for `xtallography.lattices.hexagonal` module
"""
# --- Imports

# Standard library
import unittest

# External packages
import pytest

# Local packages/modules
from xtallography import _JL
from xtallography.lattices import LatticeSystem, Centering
from xtallography.lattices import HexagonalUnitCell, TetragonalUnitCell


# --- Test Suites


class test_xtallography_lattice_hexagonal(unittest.TestCase):
    """
    Test suite for the `xtallography.lattice.hexagonal` module
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

        unit_cell = HexagonalUnitCell(a, c)
        assert unit_cell.a == a
        assert unit_cell.c == c
        assert unit_cell.lattice_system == LatticeSystem.HEXAGONAL
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
            HexagonalUnitCell(invalid_value, c)

        expected_error = f"`a` must be positive. (a={invalid_value})"
        assert expected_error in str(exception_info)

        # a = 0
        invalid_value = 0
        with pytest.raises(ValueError) as exception_info:
            HexagonalUnitCell(invalid_value, c)

        expected_error = f"`a` must be positive. (a={invalid_value})"
        assert expected_error in str(exception_info)

        # c < 0
        invalid_value = -1
        with pytest.raises(ValueError) as exception_info:
            HexagonalUnitCell(a, invalid_value)

        expected_error = f"`c` must be positive. (c={invalid_value})"
        assert expected_error in str(exception_info)

        # c = 0
        invalid_value = 0
        with pytest.raises(ValueError) as exception_info:
            HexagonalUnitCell(a, invalid_value)

        expected_error = f"`c` must be positive. (c={invalid_value})"
        assert expected_error in str(exception_info)

    @staticmethod
    def test_to_julia():
        """
        Test `to_julia()`.
        """
        # --- Preparations

        a = 1
        c = 4

        unit_cell = HexagonalUnitCell(a, c)

        # --- Tests

        unit_cell_jl = unit_cell.to_julia()
        assert _JL.isa(unit_cell_jl, _JL.UnitCell)
        assert _JL.isa(unit_cell_jl.lattice_constants, _JL.HexagonalLatticeConstants)

    @staticmethod
    def test_from_julia():
        """
        Test `from_julia()`.
        """
        # --- Preparations

        # lattice constants
        a = 1
        c = 3

        # --- Tests

        # basic usage
        unit_cell_jl = _JL.UnitCell(_JL.HexagonalLatticeConstants(a, c), _JL.primitive)
        unit_cell = HexagonalUnitCell.from_julia(unit_cell_jl)
        assert unit_cell == HexagonalUnitCell(a, c)

    @staticmethod
    def test_repr():
        """
        Test `__repr__()`.
        """
        # --- Preparations

        # lattice constants
        a = 1
        c = 3

        # --- Tests

        # centering = primitive
        unit_cell = HexagonalUnitCell(a, c)
        assert str(unit_cell) == f"HexagonalUnitCell(a={a},c={c})"

    @staticmethod
    def test_eq():
        """
        Test `__eq__()`.
        """
        # --- Preparations

        # lattice constants
        a = 1
        c = 3

        # --- Tests

        # lattice constants are the same
        unit_cell_1 = HexagonalUnitCell(a, c)
        unit_cell_2 = HexagonalUnitCell(a, c)
        assert unit_cell_1 == unit_cell_2

        # lattice constants are the different
        unit_cell_1 = HexagonalUnitCell(a + 1, c)
        unit_cell_2 = HexagonalUnitCell(a, c)
        assert unit_cell_1 != unit_cell_2

        # types are different
        unit_cell_1 = HexagonalUnitCell(a, c)
        unit_cell_2 = TetragonalUnitCell(a, c)
        assert unit_cell_1 != unit_cell_2
