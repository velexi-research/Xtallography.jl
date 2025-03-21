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
Unit tests for `xtallography.lattices.rhombohedral` module
"""
# --- Imports

# Standard library
import math
import unittest

# External packages
import pytest

# Local packages/modules
from xtallography import _JL
from xtallography.lattices import LatticeSystem, Centering
from xtallography.lattices import RhombohedralUnitCell, TetragonalUnitCell


# --- Test Suites


class test_xtallography_lattice_rhombohedral(unittest.TestCase):
    """
    Test suite for the `xtallography.lattice.rhombohedral` module
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
        alpha = 0.1

        # --- Tests

        unit_cell = RhombohedralUnitCell(a, alpha)

        assert unit_cell.a == a
        assert unit_cell.alpha == alpha
        assert unit_cell.lattice_system == LatticeSystem.RHOMBOHEDRAL
        assert unit_cell.centering == Centering.PRIMITIVE

    @staticmethod
    def test_init_invalid_args():
        """
        Test `__init__()`: invalid arguments.
        """
        # --- Preparations

        a = 1
        alpha = 0.1

        # --- Tests

        # ------ Invalid `a`

        # a < 0
        invalid_value = -1
        with pytest.raises(ValueError) as exception_info:
            RhombohedralUnitCell(invalid_value, alpha)

        expected_error = f"`a` must be positive. (a={invalid_value})"
        assert expected_error in str(exception_info)

        # a = 0
        invalid_value = 0
        with pytest.raises(ValueError) as exception_info:
            RhombohedralUnitCell(invalid_value, alpha)

        expected_error = f"`a` must be positive. (a={invalid_value})"
        assert expected_error in str(exception_info)

        # ------ Invalid `alpha`

        # alpha < 0
        invalid_value = -1
        with pytest.raises(ValueError) as exception_info:
            RhombohedralUnitCell(a, invalid_value)

        expected_error = (
            f"`alpha` must lie in the interval (0, 2 pi). (alpha={invalid_value})"
        )
        assert expected_error in str(exception_info)

        # alpha = 0
        invalid_value = 0
        with pytest.raises(ValueError) as exception_info:
            RhombohedralUnitCell(a, invalid_value)

        expected_error = (
            f"`alpha` must lie in the interval (0, 2 pi). (alpha={invalid_value})"
        )
        assert expected_error in str(exception_info)

        # alpha = 2 pi
        invalid_value = 2 * math.pi
        with pytest.raises(ValueError) as exception_info:
            RhombohedralUnitCell(a, invalid_value)

        expected_error = (
            f"`alpha` must lie in the interval (0, 2 pi). (alpha={invalid_value})"
        )
        assert expected_error in str(exception_info)

        # alpha > 2 pi
        invalid_value = 7
        with pytest.raises(ValueError) as exception_info:
            RhombohedralUnitCell(a, invalid_value)

        expected_error = (
            f"`alpha` must lie in the interval (0, 2 pi). (alpha={invalid_value})"
        )
        assert expected_error in str(exception_info)

    @staticmethod
    def test_to_julia():
        """
        Test `to_julia()`.
        """
        # --- Preparations

        a = 1
        alpha = 0.1

        unit_cell = RhombohedralUnitCell(a, alpha)

        # --- Tests

        unit_cell_jl = unit_cell.to_julia()
        assert _JL.isa(unit_cell_jl, _JL.UnitCell)
        assert _JL.isa(unit_cell_jl.lattice_constants, _JL.RhombohedralLatticeConstants)

    @staticmethod
    def test_from_julia():
        """
        Test `from_julia()`.
        """
        # --- Preparations

        # lattice constants
        a = 1
        alpha = 0.1

        # --- Tests

        # basic usage
        unit_cell_jl = _JL.UnitCell(
            _JL.RhombohedralLatticeConstants(a, alpha), _JL.primitive
        )
        unit_cell = RhombohedralUnitCell.from_julia(unit_cell_jl)
        assert unit_cell == RhombohedralUnitCell(a, alpha)

    @staticmethod
    def test_repr():
        """
        Test `__repr__()`.
        """
        # --- Preparations

        # lattice constants
        a = 1
        alpha = 0.1

        # --- Tests

        # centering = primitive
        unit_cell = RhombohedralUnitCell(a, alpha)
        assert str(unit_cell) == f"RhombohedralUnitCell(a={a},alpha={alpha})"

    @staticmethod
    def test_eq():
        """
        Test `__eq__()`.
        """
        # --- Preparations

        # lattice constants
        a = 1
        alpha = 0.1

        # --- Tests

        # lattice constants are the same
        unit_cell_1 = RhombohedralUnitCell(a, alpha)
        unit_cell_2 = RhombohedralUnitCell(a, alpha)
        assert unit_cell_1 == unit_cell_2

        # lattice constants are the different
        unit_cell_1 = RhombohedralUnitCell(a + 1, alpha)
        unit_cell_2 = RhombohedralUnitCell(a, alpha)
        assert unit_cell_1 != unit_cell_2

        # types are different
        unit_cell_1 = RhombohedralUnitCell(a, alpha)
        unit_cell_2 = TetragonalUnitCell(a, a + 1)
        assert unit_cell_1 != unit_cell_2
