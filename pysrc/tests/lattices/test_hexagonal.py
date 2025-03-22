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
    def test_from_julia_invalid_arguments():
        """
        Test `from_julia()`.
        """
        # --- Tests

        # unit_cell_jl not a Julia UnitCell object
        unit_cell_jl_invalid = "not Julia UnitCell object"
        with pytest.raises(ValueError) as exception_info:
            HexagonalUnitCell.from_julia(unit_cell_jl_invalid)

        expected_error = (
            "`unit_cell_jl` must be a Julia `UnitCell` object. "
            f"(unit_cell_jl={unit_cell_jl_invalid})."
        )
        assert expected_error in str(exception_info)

        # unit_cell_jl is not for a hexagonal unit cell
        unit_cell_jl_invalid = _JL.UnitCell(_JL.CubicLatticeConstants(1), _JL.primitive)
        with pytest.raises(ValueError) as exception_info:
            HexagonalUnitCell.from_julia(unit_cell_jl_invalid)

        expected_error = (
            "`unit_cell_jl` must be a Julia `UnitCell` object for hexagonal "
            f"unit cell. (unit_cell_jl={unit_cell_jl_invalid})."
        )
        assert expected_error in str(exception_info)

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

        # ------ types differ

        unit_cell_1 = HexagonalUnitCell(a, c)
        unit_cell_2 = TetragonalUnitCell(a, c)
        assert unit_cell_1 != unit_cell_2

        # ------ lattice constants are the same

        unit_cell_1 = HexagonalUnitCell(a, c)
        unit_cell_2 = HexagonalUnitCell(a, c)
        assert unit_cell_1 == unit_cell_2

        # ------ lattice constants differ

        # `a` values differ
        unit_cell_1 = HexagonalUnitCell(a, c)
        unit_cell_2 = HexagonalUnitCell(a + 1, c)
        assert unit_cell_1 != unit_cell_2

        # `c` values differ
        unit_cell_1 = HexagonalUnitCell(a, c)
        unit_cell_2 = HexagonalUnitCell(a, c + 2)
        assert unit_cell_1 != unit_cell_2

    @staticmethod
    def test_isclose():
        """
        Test `isclose()`.
        """
        # --- Preparations

        # lattice constants
        a = 1
        c = 3

        # --- Tests

        # ------ types differ

        unit_cell_1 = HexagonalUnitCell(a, c)
        unit_cell_2 = TetragonalUnitCell(a, a + 1)
        assert not unit_cell_1.isclose(unit_cell_2)

        # ------ `a`

        # `a` values the equal to within tolerance
        unit_cell_1 = HexagonalUnitCell(a, c)
        unit_cell_2 = HexagonalUnitCell(a + 0.1, c)
        assert unit_cell_1.isclose(unit_cell_2, atol=0.2)

        # `a` values differ by more than tolerance
        unit_cell_1 = HexagonalUnitCell(a, c)
        unit_cell_2 = HexagonalUnitCell(a + 1, c)
        assert not unit_cell_1.isclose(unit_cell_2, atol=0.2)

        # ------ `c`

        # `c` values the equal to within tolerance
        unit_cell_1 = HexagonalUnitCell(a, c)
        unit_cell_2 = HexagonalUnitCell(a, c - 0.1)
        assert unit_cell_1.isclose(unit_cell_2, atol=0.2)

        # `c` values differ by more than tolerance
        unit_cell_1 = HexagonalUnitCell(a, c)
        unit_cell_2 = HexagonalUnitCell(a, c + 0.3)
        assert not unit_cell_1.isclose(unit_cell_2, atol=0.2)
