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
Unit tests for `xtallography.lattices.monoclinic` module
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
from xtallography.lattices import MonoclinicUnitCell, TetragonalUnitCell


# --- Test Suites


class test_xtallography_lattice_monoclinic(unittest.TestCase):
    """
    Test suite for the `xtallography.lattice.monoclinic` module
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
        beta = 0.2

        # --- Tests

        # default centering
        unit_cell = MonoclinicUnitCell(a, b, c, beta)
        assert unit_cell.a == a
        assert unit_cell.b == b
        assert unit_cell.c == c
        assert unit_cell.beta == beta
        assert unit_cell.lattice_system == LatticeSystem.MONOCLINIC
        assert unit_cell.centering == Centering.PRIMITIVE

        # centering = primitive
        unit_cell = MonoclinicUnitCell(a, b, c, beta, centering=Centering.PRIMITIVE)
        assert unit_cell.a == a
        assert unit_cell.b == b
        assert unit_cell.c == c
        assert unit_cell.beta == beta
        assert unit_cell.lattice_system == LatticeSystem.MONOCLINIC
        assert unit_cell.centering == Centering.PRIMITIVE

        # centering = body-centered
        unit_cell = MonoclinicUnitCell(a, b, c, beta, centering=Centering.BODY_CENTERED)
        assert unit_cell.a == a
        assert unit_cell.b == b
        assert unit_cell.c == c
        assert unit_cell.beta == beta
        assert unit_cell.lattice_system == LatticeSystem.MONOCLINIC
        assert unit_cell.centering == Centering.BODY_CENTERED

        # centering = base-centered
        unit_cell = MonoclinicUnitCell(a, b, c, beta, centering=Centering.BASE_CENTERED)
        assert unit_cell.a == a
        assert unit_cell.b == b
        assert unit_cell.c == c
        assert unit_cell.beta == beta
        assert unit_cell.lattice_system == LatticeSystem.MONOCLINIC
        assert unit_cell.centering == Centering.BASE_CENTERED

        # centering argument is a str
        unit_cell = MonoclinicUnitCell(a, b, c, beta, centering="primitive")

        assert unit_cell.a == a
        assert unit_cell.b == b
        assert unit_cell.c == c
        assert unit_cell.beta == beta
        assert isinstance(unit_cell.lattice_system, LatticeSystem)
        assert unit_cell.lattice_system == LatticeSystem.MONOCLINIC
        assert isinstance(unit_cell.centering, Centering)
        assert unit_cell.centering == Centering.PRIMITIVE

    @staticmethod
    def test_init_invalid_args():
        """
        Test `__init__()`: invalid arguments.
        """
        # --- Preparations

        a = 1
        b = 2
        c = 3
        beta = 0.2

        # --- Tests

        # ------ Invalid `a`

        # a < 0
        invalid_value = -1
        with pytest.raises(ValueError) as exception_info:
            MonoclinicUnitCell(invalid_value, b, c, beta)

        expected_error = f"`a` must be positive. (a={invalid_value})"
        assert expected_error in str(exception_info)

        # a = 0
        invalid_value = 0
        with pytest.raises(ValueError) as exception_info:
            MonoclinicUnitCell(invalid_value, b, c, beta)

        expected_error = f"`a` must be positive. (a={invalid_value})"
        assert expected_error in str(exception_info)

        # ------ Invalid `b`

        # b < 0
        invalid_value = -1
        with pytest.raises(ValueError) as exception_info:
            MonoclinicUnitCell(a, invalid_value, c, beta)

        expected_error = f"`b` must be positive. (b={invalid_value})"
        assert expected_error in str(exception_info)

        # b = 0
        invalid_value = 0
        with pytest.raises(ValueError) as exception_info:
            MonoclinicUnitCell(a, invalid_value, c, beta)

        expected_error = f"`b` must be positive. (b={invalid_value})"
        assert expected_error in str(exception_info)

        # ------ Invalid `c`

        # c < 0
        invalid_value = -1
        with pytest.raises(ValueError) as exception_info:
            MonoclinicUnitCell(a, b, invalid_value, beta)

        expected_error = f"`c` must be positive. (c={invalid_value})"
        assert expected_error in str(exception_info)

        # c = 0
        invalid_value = 0
        with pytest.raises(ValueError) as exception_info:
            MonoclinicUnitCell(a, b, invalid_value, beta)

        expected_error = f"`c` must be positive. (c={invalid_value})"
        assert expected_error in str(exception_info)

        # ------ Invalid `beta`

        # beta < 0
        invalid_value = -1
        with pytest.raises(ValueError) as exception_info:
            MonoclinicUnitCell(a, b, c, invalid_value)

        expected_error = (
            f"`beta` must lie in the interval (0, pi). (beta={invalid_value})"
        )
        assert expected_error in str(exception_info)

        # beta = 0
        invalid_value = 0
        with pytest.raises(ValueError) as exception_info:
            MonoclinicUnitCell(a, b, c, invalid_value)

        expected_error = (
            f"`beta` must lie in the interval (0, pi). (beta={invalid_value})"
        )
        assert expected_error in str(exception_info)

        # beta = pi
        invalid_value = math.pi
        with pytest.raises(ValueError) as exception_info:
            MonoclinicUnitCell(a, b, c, invalid_value)

        expected_error = (
            f"`beta` must lie in the interval (0, pi). (beta={invalid_value})"
        )
        assert expected_error in str(exception_info)

        # beta > pi
        invalid_value = math.pi + 0.1
        with pytest.raises(ValueError) as exception_info:
            MonoclinicUnitCell(a, b, c, invalid_value)

        expected_error = (
            f"`beta` must lie in the interval (0, pi). (beta={invalid_value})"
        )
        assert expected_error in str(exception_info)

        # ------ Invalid `centering`

        # face-centered
        invalid_value = "face-centered"
        with pytest.raises(ValueError) as exception_info:
            MonoclinicUnitCell(a, b, c, beta, centering=invalid_value)

        expected_error = (
            "('monoclinic', 'face-centered') is not a valid Bravais lattice."
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
        beta = 0.1

        unit_cell = MonoclinicUnitCell(a, b, c, beta)

        # --- Tests

        unit_cell_jl = unit_cell.to_julia()
        assert self.jl.isa(unit_cell_jl, self.jl.UnitCell)
        assert self.jl.isa(
            unit_cell_jl.lattice_constants, self.jl.MonoclinicLatticeConstants
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
        beta = 0.1

        # --- Tests

        # centering = primitive
        unit_cell_jl = self.jl.UnitCell(
            self.jl.MonoclinicLatticeConstants(a, b, c, beta), self.jl.primitive
        )
        unit_cell = MonoclinicUnitCell.from_julia(unit_cell_jl)
        assert unit_cell == MonoclinicUnitCell(
            a, b, c, beta, centering=Centering.PRIMITIVE
        )

        # centering = body_centered
        unit_cell_jl = self.jl.UnitCell(
            self.jl.MonoclinicLatticeConstants(a, b, c, beta), self.jl.body_centered
        )
        unit_cell = MonoclinicUnitCell.from_julia(unit_cell_jl)
        assert unit_cell == MonoclinicUnitCell(
            a, b, c, beta, centering=Centering.BODY_CENTERED
        )

        # centering = base_centered
        unit_cell_jl = self.jl.UnitCell(
            self.jl.MonoclinicLatticeConstants(a, b, c, beta), self.jl.base_centered
        )
        unit_cell = MonoclinicUnitCell.from_julia(unit_cell_jl)
        assert unit_cell == MonoclinicUnitCell(
            a, b, c, beta, centering=Centering.BASE_CENTERED
        )

    def test_from_julia_invalid_arguments(self):
        """
        Test `from_julia()`.
        """
        # --- Tests

        # unit_cell_jl not a Julia UnitCell object
        unit_cell_jl_invalid = "not Julia UnitCell object"
        with pytest.raises(ValueError) as exception_info:
            MonoclinicUnitCell.from_julia(unit_cell_jl_invalid)

        expected_error = (
            "`unit_cell_jl` must be a Julia `UnitCell` object. "
            f"(unit_cell_jl={unit_cell_jl_invalid})."
        )
        assert expected_error in str(exception_info)

        # unit_cell_jl is not for a monoclinic unit cell
        unit_cell_jl_invalid = self.jl.UnitCell(
            self.jl.CubicLatticeConstants(1), self.jl.primitive
        )
        with pytest.raises(ValueError) as exception_info:
            MonoclinicUnitCell.from_julia(unit_cell_jl_invalid)

        expected_error = (
            "`unit_cell_jl` must be a Julia `UnitCell` object for monoclinic "
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
        beta = 0.1

        # --- Tests

        # centering = primitive
        unit_cell = MonoclinicUnitCell(a, b, c, beta, centering=Centering.PRIMITIVE)
        assert str(unit_cell) == (
            f"MonoclinicUnitCell(a={a},b={b},c={c},beta={beta},"
            f"centering='primitive')"
        )

        # centering = body_centered
        unit_cell = MonoclinicUnitCell(a, b, c, beta, centering=Centering.BODY_CENTERED)
        assert str(unit_cell) == (
            f"MonoclinicUnitCell(a={a},b={b},c={c},beta={beta},"
            f"centering='body_centered')"
        )

        # centering = base_centered
        unit_cell = MonoclinicUnitCell(a, b, c, beta, centering=Centering.BASE_CENTERED)
        assert str(unit_cell) == (
            f"MonoclinicUnitCell(a={a},b={b},c={c},beta={beta},"
            f"centering='base_centered')"
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
        beta = 0.1

        # --- Tests

        # ------ types differ

        unit_cell_1 = MonoclinicUnitCell(a, b, c, beta)
        unit_cell_2 = TetragonalUnitCell(a, a + 1)
        assert unit_cell_1 != unit_cell_2

        # ------ lattice constants are the same

        # centerings are the same
        unit_cell_1 = MonoclinicUnitCell(a, b, c, beta)
        unit_cell_2 = MonoclinicUnitCell(a, b, c, beta, centering="primitive")
        assert unit_cell_1 == unit_cell_2

        # centerings are different
        unit_cell_1 = MonoclinicUnitCell(a, b, c, beta)
        unit_cell_2 = MonoclinicUnitCell(a, b, c, beta, centering="body_centered")
        assert unit_cell_1 != unit_cell_2

        # ------ lattice constants differ, centerings are the same

        # `a` values differ
        unit_cell_1 = MonoclinicUnitCell(a, b, c, beta)
        unit_cell_2 = MonoclinicUnitCell(a + 1, b, c, beta, centering="primitive")
        assert unit_cell_1 != unit_cell_2

        # `b` values differ
        unit_cell_1 = MonoclinicUnitCell(a, b, c, beta)
        unit_cell_2 = MonoclinicUnitCell(a, b + 2, c, beta, centering="primitive")
        assert unit_cell_1 != unit_cell_2

        # `c` values differ
        unit_cell_1 = MonoclinicUnitCell(a, b, c, beta)
        unit_cell_2 = MonoclinicUnitCell(a, b, c + 3, beta, centering="primitive")
        assert unit_cell_1 != unit_cell_2

        # `beta` values differ
        unit_cell_1 = MonoclinicUnitCell(a, b, c, beta)
        unit_cell_2 = MonoclinicUnitCell(a, b, c, beta + 0.1, centering="primitive")
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
        beta = 0.1

        # --- Tests

        # ------ types differ

        unit_cell_1 = MonoclinicUnitCell(a, b, c, beta)
        unit_cell_2 = TetragonalUnitCell(a, a + 1)
        assert not unit_cell_1.isclose(unit_cell_2)

        # ------ `a`

        # `a` values differ by less than tolerance, centerings are the same
        unit_cell_1 = MonoclinicUnitCell(a, b, c, beta)
        unit_cell_2 = MonoclinicUnitCell(a + 0.1, b, c, beta, centering="primitive")
        assert unit_cell_1.isclose(unit_cell_2, atol=0.2)

        # `a` values differ by less than tolerance, centerings differ
        unit_cell_1 = MonoclinicUnitCell(a, b, c, beta)
        unit_cell_2 = MonoclinicUnitCell(a + 0.1, b, c, beta, centering="body_centered")
        assert not unit_cell_1.isclose(unit_cell_2, atol=0.2)

        # `a` values differ by more than tolerance, centerings are the same
        unit_cell_1 = MonoclinicUnitCell(a, b, c, beta)
        unit_cell_2 = MonoclinicUnitCell(a + 1, b, c, beta, centering="primitive")
        assert not unit_cell_1.isclose(unit_cell_2, atol=0.2)

        # ------ `b`

        # `b` values differ by less than tolerance, centerings are the same
        unit_cell_1 = MonoclinicUnitCell(a, b, c, beta)
        unit_cell_2 = MonoclinicUnitCell(a, b + 0.1, c, beta, centering="primitive")
        assert unit_cell_1.isclose(unit_cell_2, atol=0.2)

        # `b` values differ by less than tolerance, centerings differ
        unit_cell_1 = MonoclinicUnitCell(a, b, c, beta)
        unit_cell_2 = MonoclinicUnitCell(
            a, b - 0.15, c, beta, centering="body_centered"
        )
        assert not unit_cell_1.isclose(unit_cell_2, atol=0.2)

        # `b` values differ by more than tolerance, centerings are the same
        unit_cell_1 = MonoclinicUnitCell(a, b, c, beta)
        unit_cell_2 = MonoclinicUnitCell(a, b + 2, c, beta, centering="primitive")
        assert not unit_cell_1.isclose(unit_cell_2, atol=0.2)

        # ------ `c`

        # `c` values differ by less than tolerance, centerings are the same
        unit_cell_1 = MonoclinicUnitCell(a, b, c, beta)
        unit_cell_2 = MonoclinicUnitCell(a, b, c - 0.05, beta, centering="primitive")
        assert unit_cell_1.isclose(unit_cell_2, atol=0.2)

        # `c` values differ by less than tolerance, centerings differ
        unit_cell_1 = MonoclinicUnitCell(a, b, c, beta)
        unit_cell_2 = MonoclinicUnitCell(a, b, c + 0.1, beta, centering="body_centered")
        assert not unit_cell_1.isclose(unit_cell_2, atol=0.2)

        # `c` values differ by more than tolerance, centerings are the same
        unit_cell_1 = MonoclinicUnitCell(a, b, c, beta)
        unit_cell_2 = MonoclinicUnitCell(a, b, c + 3, beta, centering="primitive")
        assert not unit_cell_1.isclose(unit_cell_2, atol=0.2)

        # ------ `beta`

        # `beta` values differ by less than tolerance, centerings are the same
        unit_cell_1 = MonoclinicUnitCell(a, b, c, beta)
        unit_cell_2 = MonoclinicUnitCell(a, b, c, beta + 0.05, centering="primitive")
        assert unit_cell_1.isclose(unit_cell_2, atol=0.2)

        # `beta` values differ by less than tolerance, centerings differ
        unit_cell_1 = MonoclinicUnitCell(a, b, c, beta)
        unit_cell_2 = MonoclinicUnitCell(
            a, b, c, beta - 0.05, centering="body_centered"
        )
        assert not unit_cell_1.isclose(unit_cell_2, atol=0.2)

        # `beta` values differ by more than tolerance, centerings are the same
        unit_cell_1 = MonoclinicUnitCell(a, b, c, beta)
        unit_cell_2 = MonoclinicUnitCell(a, b, c, beta + 0.25, centering="primitive")
        assert not unit_cell_1.isclose(unit_cell_2, atol=0.2)
