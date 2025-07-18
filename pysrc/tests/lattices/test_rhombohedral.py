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
import juliacall
import pytest

# Local packages/modules
from xtallography.lattices import LatticeSystem, Centering
from xtallography.lattices import RhombohedralUnitCell, TetragonalUnitCell


# --- Test Suites


class test_xtallography_lattice_rhombohedral(unittest.TestCase):
    """
    Test suite for the `RhombohedralUnitCell` class
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

    def test_to_julia(self):
        """
        Test `to_julia()`.
        """
        # --- Preparations

        a = 1
        alpha = 0.1

        unit_cell = RhombohedralUnitCell(a, alpha)

        # --- Tests

        unit_cell_jl = unit_cell.to_julia()
        assert self.jl.isa(unit_cell_jl, self.jl.UnitCell)
        assert self.jl.isa(
            unit_cell_jl.lattice_constants, self.jl.RhombohedralLatticeConstants
        )

    def test_from_julia(self):
        """
        Test `from_julia()`.
        """
        # --- Preparations

        # lattice constants
        a = 1
        alpha = 0.1

        # --- Tests

        # basic usage
        unit_cell_jl = self.jl.UnitCell(
            self.jl.RhombohedralLatticeConstants(a, alpha), self.jl.primitive
        )
        unit_cell = RhombohedralUnitCell.from_julia(unit_cell_jl)
        assert unit_cell == RhombohedralUnitCell(a, alpha)

    def test_from_julia_invalid_arguments(self):
        """
        Test `from_julia()`.
        """
        # --- Tests

        # unit_cell_jl not a Julia UnitCell object
        unit_cell_jl_invalid = "not Julia UnitCell object"
        with pytest.raises(ValueError) as exception_info:
            RhombohedralUnitCell.from_julia(unit_cell_jl_invalid)

        expected_error = (
            "`unit_cell_jl` must be a Julia `UnitCell` object. "
            f"(unit_cell_jl={unit_cell_jl_invalid})."
        )
        assert expected_error in str(exception_info)

        # unit_cell_jl is not for a rhombohedral unit cell
        unit_cell_jl_invalid = self.jl.UnitCell(
            self.jl.CubicLatticeConstants(1), self.jl.primitive
        )
        with pytest.raises(ValueError) as exception_info:
            RhombohedralUnitCell.from_julia(unit_cell_jl_invalid)

        expected_error = (
            "`unit_cell_jl` must be a Julia `UnitCell` object for rhombohedral "
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

        # ------ types differ

        unit_cell_1 = RhombohedralUnitCell(a, alpha)
        unit_cell_2 = TetragonalUnitCell(a, a + 2)
        assert unit_cell_1 != unit_cell_2

        # ------ lattice constants are the same

        unit_cell_1 = RhombohedralUnitCell(a, alpha)
        unit_cell_2 = RhombohedralUnitCell(a, alpha)
        assert unit_cell_1 == unit_cell_2

        # ------ lattice constants differ

        # `a` values differ
        unit_cell_1 = RhombohedralUnitCell(a, alpha)
        unit_cell_2 = RhombohedralUnitCell(a + 1, alpha)
        assert unit_cell_1 != unit_cell_2

        # `alpha` values differ
        unit_cell_1 = RhombohedralUnitCell(a, alpha)
        unit_cell_2 = RhombohedralUnitCell(a, alpha + 0.1)
        assert unit_cell_1 != unit_cell_2

    @staticmethod
    def test_isclose():
        """
        Test `isclose()`.
        """
        # --- Preparations

        # lattice constants
        a = 1
        alpha = 0.1

        # --- Tests

        # ------ types differ

        unit_cell_1 = RhombohedralUnitCell(a, alpha)
        unit_cell_2 = TetragonalUnitCell(a, a + 1)
        assert not unit_cell_1.isclose(unit_cell_2)

        # ------ `a`

        # `a` values the equal to within tolerance
        unit_cell_1 = RhombohedralUnitCell(a, alpha)
        unit_cell_2 = RhombohedralUnitCell(a + 0.1, alpha)
        assert unit_cell_1.isclose(unit_cell_2, atol=0.2)

        # `a` values differ by more than tolerance
        unit_cell_1 = RhombohedralUnitCell(a, alpha)
        unit_cell_2 = RhombohedralUnitCell(a + 1, alpha)
        assert not unit_cell_1.isclose(unit_cell_2, atol=0.2)

        # ------ `alpha`

        # `alpha` values the equal to within tolerance
        unit_cell_1 = RhombohedralUnitCell(a, alpha)
        unit_cell_2 = RhombohedralUnitCell(a, alpha + 0.1)
        assert unit_cell_1.isclose(unit_cell_2, atol=0.2)

        # `alpha` values differ by more than tolerance
        unit_cell_1 = RhombohedralUnitCell(a, alpha)
        unit_cell_2 = RhombohedralUnitCell(a, alpha + 0.3)
        assert not unit_cell_1.isclose(unit_cell_2, atol=0.2)
