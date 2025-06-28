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
Unit tests for `xtallography.lattices.triclinic` module
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
from xtallography.lattices import TriclinicUnitCell, TetragonalUnitCell


# --- Test Suites


class test_xtallography_lattice_triclinic(unittest.TestCase):
    """
    Test suite for the `TriclinicUnitCell` class
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
        b = 2
        c = 3
        alpha = 0.1
        beta = 0.2
        gamma = 0.3

        # --- Tests

        lattice_constants = TriclinicUnitCell(a, b, c, alpha, beta, gamma)

        assert lattice_constants.a == a
        assert lattice_constants.b == b
        assert lattice_constants.c == c
        assert lattice_constants.alpha == alpha
        assert lattice_constants.beta == beta
        assert lattice_constants.gamma == gamma
        assert lattice_constants.lattice_system == LatticeSystem.TRICLINIC
        assert lattice_constants.centering == Centering.PRIMITIVE

    @staticmethod
    def test_init_invalid_args():
        """
        Test `__init__()`: invalid arguments.
        """
        # --- Preparations

        a = 1
        b = 2
        c = 3
        alpha = 0.1
        beta = 0.2
        gamma = 0.3

        # --- Tests

        # ------ Invalid `a`

        # a < 0
        invalid_value = -1
        with pytest.raises(ValueError) as exception_info:
            TriclinicUnitCell(invalid_value, b, c, alpha, beta, gamma)

        expected_error = f"`a` must be positive. (a={invalid_value})"
        assert expected_error in str(exception_info)

        # a = 0
        invalid_value = 0
        with pytest.raises(ValueError) as exception_info:
            TriclinicUnitCell(invalid_value, b, c, alpha, beta, gamma)

        expected_error = f"`a` must be positive. (a={invalid_value})"
        assert expected_error in str(exception_info)

        # ------ Invalid `b`

        # b < 0
        invalid_value = -1
        with pytest.raises(ValueError) as exception_info:
            TriclinicUnitCell(a, invalid_value, c, alpha, beta, gamma)

        expected_error = f"`b` must be positive. (b={invalid_value})"
        assert expected_error in str(exception_info)

        # b = 0
        invalid_value = 0
        with pytest.raises(ValueError) as exception_info:
            TriclinicUnitCell(a, invalid_value, c, alpha, beta, gamma)

        expected_error = f"`b` must be positive. (b={invalid_value})"
        assert expected_error in str(exception_info)

        # ------ Invalid `c`

        # c < 0
        invalid_value = -1
        with pytest.raises(ValueError) as exception_info:
            TriclinicUnitCell(a, b, invalid_value, alpha, beta, gamma)

        expected_error = f"`c` must be positive. (c={invalid_value})"
        assert expected_error in str(exception_info)

        # c = 0
        invalid_value = 0
        with pytest.raises(ValueError) as exception_info:
            TriclinicUnitCell(a, b, invalid_value, alpha, beta, gamma)

        expected_error = f"`c` must be positive. (c={invalid_value})"
        assert expected_error in str(exception_info)

        # ------ Invalid `alpha`

        # alpha < 0
        invalid_value = -1
        with pytest.raises(ValueError) as exception_info:
            TriclinicUnitCell(a, b, c, invalid_value, beta, gamma)

        expected_error = (
            f"`alpha` must lie in the interval (0, 2 pi). (alpha={invalid_value})"
        )
        assert expected_error in str(exception_info)

        # alpha = 0
        invalid_value = 0
        with pytest.raises(ValueError) as exception_info:
            TriclinicUnitCell(a, b, c, invalid_value, beta, gamma)

        expected_error = (
            f"`alpha` must lie in the interval (0, 2 pi). (alpha={invalid_value})"
        )
        assert expected_error in str(exception_info)

        # alpha = 2 pi
        invalid_value = 2 * math.pi
        with pytest.raises(ValueError) as exception_info:
            TriclinicUnitCell(a, b, c, invalid_value, beta, gamma)

        expected_error = (
            f"`alpha` must lie in the interval (0, 2 pi). (alpha={invalid_value})"
        )
        assert expected_error in str(exception_info)

        # alpha > 2 pi
        invalid_value = 7
        with pytest.raises(ValueError) as exception_info:
            TriclinicUnitCell(a, b, c, invalid_value, beta, gamma)

        expected_error = (
            f"`alpha` must lie in the interval (0, 2 pi). (alpha={invalid_value})"
        )
        assert expected_error in str(exception_info)

        # ------ Invalid `beta`

        # beta < 0
        invalid_value = -1
        with pytest.raises(ValueError) as exception_info:
            TriclinicUnitCell(a, b, c, alpha, invalid_value, gamma)

        expected_error = (
            f"`beta` must lie in the interval (0, 2 pi). (beta={invalid_value})"
        )
        assert expected_error in str(exception_info)

        # beta = 0
        invalid_value = 0
        with pytest.raises(ValueError) as exception_info:
            TriclinicUnitCell(a, b, c, alpha, invalid_value, gamma)

        expected_error = (
            f"`beta` must lie in the interval (0, 2 pi). (beta={invalid_value})"
        )
        assert expected_error in str(exception_info)

        # beta = 2 pi
        invalid_value = 2 * math.pi
        with pytest.raises(ValueError) as exception_info:
            TriclinicUnitCell(a, b, c, alpha, invalid_value, gamma)

        expected_error = (
            f"`beta` must lie in the interval (0, 2 pi). (beta={invalid_value})"
        )
        assert expected_error in str(exception_info)

        # beta > 2 pi
        invalid_value = 7
        with pytest.raises(ValueError) as exception_info:
            TriclinicUnitCell(a, b, c, alpha, invalid_value, gamma)

        expected_error = (
            f"`beta` must lie in the interval (0, 2 pi). (beta={invalid_value})"
        )
        assert expected_error in str(exception_info)

        # ------ Invalid `gamma`

        # gamma < 0
        invalid_value = -1
        with pytest.raises(ValueError) as exception_info:
            TriclinicUnitCell(a, b, c, alpha, beta, invalid_value)

        expected_error = (
            f"`gamma` must lie in the interval (0, 2 pi). (gamma={invalid_value})"
        )
        assert expected_error in str(exception_info)

        # gamma = 0
        invalid_value = 0
        with pytest.raises(ValueError) as exception_info:
            TriclinicUnitCell(a, b, c, alpha, beta, invalid_value)

        expected_error = (
            f"`gamma` must lie in the interval (0, 2 pi). (gamma={invalid_value})"
        )
        assert expected_error in str(exception_info)

        # gamma = 2 pi
        invalid_value = 2 * math.pi
        with pytest.raises(ValueError) as exception_info:
            TriclinicUnitCell(a, b, c, alpha, beta, invalid_value)

        expected_error = (
            f"`gamma` must lie in the interval (0, 2 pi). (gamma={invalid_value})"
        )
        assert expected_error in str(exception_info)

        # gamma > 2 pi
        invalid_value = 7
        with pytest.raises(ValueError) as exception_info:
            TriclinicUnitCell(a, b, c, alpha, beta, invalid_value)

        expected_error = (
            f"`gamma` must lie in the interval (0, 2 pi). (gamma={invalid_value})"
        )
        assert expected_error in str(exception_info)

    def test_to_julia(self):
        """
        Test `to_julia()`.
        """
        # --- Preparations

        # lattice constants
        a = 1
        b = 2
        c = 3
        alpha = 0.1
        beta = 0.2
        gamma = 0.3

        unit_cell = TriclinicUnitCell(a, b, c, alpha, beta, gamma)

        # --- Tests

        unit_cell_jl = unit_cell.to_julia()
        assert self.jl.isa(unit_cell_jl, self.jl.UnitCell)
        assert self.jl.isa(
            unit_cell_jl.lattice_constants, self.jl.TriclinicLatticeConstants
        )

    def test_from_julia(self):
        """
        Test `from_julia()`.
        """
        # --- Preparations

        # lattice constants
        a = 1
        b = 2
        c = 3
        alpha = 0.1
        beta = 0.2
        gamma = 0.3

        # --- Tests

        # basic usage
        unit_cell_jl = self.jl.UnitCell(
            self.jl.TriclinicLatticeConstants(a, b, c, alpha, beta, gamma),
            self.jl.primitive,
        )
        unit_cell = TriclinicUnitCell.from_julia(unit_cell_jl)
        assert unit_cell == TriclinicUnitCell(a, b, c, alpha, beta, gamma)

    def test_from_julia_invalid_arguments(self):
        """
        Test `from_julia()`.
        """
        # --- Tests

        # unit_cell_jl not a Julia UnitCell object
        unit_cell_jl_invalid = "not Julia UnitCell object"
        with pytest.raises(ValueError) as exception_info:
            TriclinicUnitCell.from_julia(unit_cell_jl_invalid)

        expected_error = (
            "`unit_cell_jl` must be a Julia `UnitCell` object. "
            f"(unit_cell_jl={unit_cell_jl_invalid})."
        )
        assert expected_error in str(exception_info)

        # unit_cell_jl is not for a triclinic unit cell
        unit_cell_jl_invalid = self.jl.UnitCell(
            self.jl.CubicLatticeConstants(1), self.jl.primitive
        )
        with pytest.raises(ValueError) as exception_info:
            TriclinicUnitCell.from_julia(unit_cell_jl_invalid)

        expected_error = (
            "`unit_cell_jl` must be a Julia `UnitCell` object for triclinic "
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
        b = 2
        c = 3
        alpha = 0.1
        beta = 0.2
        gamma = 0.3

        # --- Tests

        # centering = primitive
        unit_cell = TriclinicUnitCell(a, b, c, alpha, beta, gamma)
        assert str(unit_cell) == (
            f"TriclinicUnitCell(a={a},b={b},c={c},"
            f"alpha={alpha},beta={beta},gamma={gamma})"
        )

    @staticmethod
    def test_eq():
        """
        Test `__eq__()`.
        """
        # --- Preparations

        # lattice constants
        a = 1
        b = 2
        c = 3
        alpha = 0.1
        beta = 0.2
        gamma = 0.3

        # --- Tests

        # ------ types differ

        unit_cell_1 = TriclinicUnitCell(a, b, c, alpha, beta, gamma)
        unit_cell_2 = TetragonalUnitCell(a, c)
        assert unit_cell_1 != unit_cell_2

        # ------ lattice constants are the same

        unit_cell_1 = TriclinicUnitCell(a, b, c, alpha, beta, gamma)
        unit_cell_2 = TriclinicUnitCell(a, b, c, alpha, beta, gamma)
        assert unit_cell_1 == unit_cell_2

        # ------ lattice constants differ

        # `a` values differ
        unit_cell_1 = TriclinicUnitCell(a, b, c, alpha, beta, gamma)
        unit_cell_2 = TriclinicUnitCell(a + 1, b, c, alpha, beta, gamma)
        assert unit_cell_1 != unit_cell_2

        # `b` values differ
        unit_cell_1 = TriclinicUnitCell(a, b, c, alpha, beta, gamma)
        unit_cell_2 = TriclinicUnitCell(a, b + 2, c, alpha, beta, gamma)
        assert unit_cell_1 != unit_cell_2

        # `c` values differ
        unit_cell_1 = TriclinicUnitCell(a, b, c, alpha, beta, gamma)
        unit_cell_2 = TriclinicUnitCell(a, b, c + 3, alpha, beta, gamma)
        assert unit_cell_1 != unit_cell_2

        # `alpha` values differ
        unit_cell_1 = TriclinicUnitCell(a, b, c, alpha, beta, gamma)
        unit_cell_2 = TriclinicUnitCell(a, b, c, alpha + 0.1, beta, gamma)
        assert unit_cell_1 != unit_cell_2

        # `beta` values differ
        unit_cell_1 = TriclinicUnitCell(a, b, c, alpha, beta, gamma)
        unit_cell_2 = TriclinicUnitCell(a, b, c, alpha, beta + 0.2, gamma)
        assert unit_cell_1 != unit_cell_2

        # `gamma` values differ
        unit_cell_1 = TriclinicUnitCell(a, b, c, alpha, beta, gamma)
        unit_cell_2 = TriclinicUnitCell(a, b, c, alpha, beta, gamma + 0.3)
        assert unit_cell_1 != unit_cell_2

    @staticmethod
    def test_isclose():
        """
        Test `isclose()`.
        """
        # --- Preparations

        # lattice constants
        a = 1
        b = 2
        c = 3
        alpha = 0.1
        beta = 0.2
        gamma = 0.3

        # --- Tests

        # ------ types differ

        unit_cell_1 = TriclinicUnitCell(a, b, c, alpha, beta, gamma)
        unit_cell_2 = TetragonalUnitCell(a, a + 1)
        assert not unit_cell_1.isclose(unit_cell_2)

        # ------ `a`

        # `a` values the equal to within tolerance
        unit_cell_1 = TriclinicUnitCell(a, b, c, alpha, beta, gamma)
        unit_cell_2 = TriclinicUnitCell(a + 0.1, b, c, alpha, beta, gamma)
        assert unit_cell_1.isclose(unit_cell_2, atol=0.2)

        # `a` values differ by more than tolerance
        unit_cell_1 = TriclinicUnitCell(a, b, c, alpha, beta, gamma)
        unit_cell_2 = TriclinicUnitCell(a + 1, b, c, alpha, beta, gamma)
        assert not unit_cell_1.isclose(unit_cell_2, atol=0.2)

        # ------ `b`

        # `b` values the equal to within tolerance
        unit_cell_1 = TriclinicUnitCell(a, b, c, alpha, beta, gamma)
        unit_cell_2 = TriclinicUnitCell(a, b + 0.1, c, alpha, beta, gamma)
        assert unit_cell_1.isclose(unit_cell_2, atol=0.2)

        # `b` values differ by more than tolerance
        unit_cell_1 = TriclinicUnitCell(a, b, c, alpha, beta, gamma)
        unit_cell_2 = TriclinicUnitCell(a, b + 1, c, alpha, beta, gamma)
        assert not unit_cell_1.isclose(unit_cell_2, atol=0.2)

        # ------ `c`

        # `c` values the equal to within tolerance
        unit_cell_1 = TriclinicUnitCell(a, b, c, alpha, beta, gamma)
        unit_cell_2 = TriclinicUnitCell(a, b, c + 0.15, alpha, beta, gamma)
        assert unit_cell_1.isclose(unit_cell_2, atol=0.2)

        # `c` values differ by more than tolerance
        unit_cell_1 = TriclinicUnitCell(a, b, c, alpha, beta, gamma)
        unit_cell_2 = TriclinicUnitCell(a, b, c + 1, alpha, beta, gamma)
        assert not unit_cell_1.isclose(unit_cell_2, atol=0.2)

        # ------ `alpha`

        # `alpha` values the equal to within tolerance
        unit_cell_1 = TriclinicUnitCell(a, b, c, alpha, beta, gamma)
        unit_cell_2 = TriclinicUnitCell(a, b, c, alpha + 0.1, beta, gamma)
        assert unit_cell_1.isclose(unit_cell_2, atol=0.2)

        # `alpha` values differ by more than tolerance
        unit_cell_1 = TriclinicUnitCell(a, b, c, alpha, beta, gamma)
        unit_cell_2 = TriclinicUnitCell(a, b, c, alpha + 0.3, beta, gamma)
        assert not unit_cell_1.isclose(unit_cell_2, atol=0.2)

        # ------ `beta`

        # `beta` values the equal to within tolerance
        unit_cell_1 = TriclinicUnitCell(a, b, c, alpha, beta, gamma)
        unit_cell_2 = TriclinicUnitCell(a, b, c, alpha, beta - 0.1, gamma)
        assert unit_cell_1.isclose(unit_cell_2, atol=0.2)

        # `beta` values differ by more than tolerance
        unit_cell_1 = TriclinicUnitCell(a, b, c, alpha, beta, gamma)
        unit_cell_2 = TriclinicUnitCell(a, b, c, alpha, beta + 0.3, gamma)
        assert not unit_cell_1.isclose(unit_cell_2, atol=0.2)

        # ------ `gamma`

        # `gamma` values the equal to within tolerance
        unit_cell_1 = TriclinicUnitCell(a, b, c, alpha, beta, gamma)
        unit_cell_2 = TriclinicUnitCell(a, b, c, alpha, beta, gamma - 0.1)
        assert unit_cell_1.isclose(unit_cell_2, atol=0.2)

        # `gamma` values differ by more than tolerance
        unit_cell_1 = TriclinicUnitCell(a, b, c, alpha, beta, gamma)
        unit_cell_2 = TriclinicUnitCell(a, b, c, alpha, beta, gamma - 0.21)
        assert not unit_cell_1.isclose(unit_cell_2, atol=0.2)
