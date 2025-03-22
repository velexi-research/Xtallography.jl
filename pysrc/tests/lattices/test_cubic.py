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
from xtallography.lattices import LatticeSystem, Centering
from xtallography.lattices import CubicUnitCell, TetragonalUnitCell


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
        assert _JL.isa(unit_cell_jl, _JL.UnitCell)
        assert _JL.isa(unit_cell_jl.lattice_constants, _JL.CubicLatticeConstants)

    @staticmethod
    def test_from_julia():
        """
        Test `from_julia()`.
        """
        # --- Preparations

        # lattice constants
        a = 1

        # --- Tests

        # centering = primitive
        unit_cell_jl = _JL.UnitCell(_JL.CubicLatticeConstants(a), _JL.primitive)
        unit_cell = CubicUnitCell.from_julia(unit_cell_jl)
        assert unit_cell == CubicUnitCell(a, centering=Centering.PRIMITIVE)

        # centering = body_centered
        unit_cell_jl = _JL.UnitCell(_JL.CubicLatticeConstants(a), _JL.body_centered)
        unit_cell = CubicUnitCell.from_julia(unit_cell_jl)
        assert unit_cell == CubicUnitCell(a, centering=Centering.BODY_CENTERED)

        # centering = face_centered
        unit_cell_jl = _JL.UnitCell(_JL.CubicLatticeConstants(a), _JL.face_centered)
        unit_cell = CubicUnitCell.from_julia(unit_cell_jl)
        assert unit_cell == CubicUnitCell(a, centering=Centering.FACE_CENTERED)

    @staticmethod
    def test_from_julia_invalid_arguments():
        """
        Test `from_julia()`.
        """
        # --- Tests

        # unit_cell_jl not a Julia UnitCell object
        unit_cell_jl_invalid = "not Julia UnitCell object"
        with pytest.raises(ValueError) as exception_info:
            CubicUnitCell.from_julia(unit_cell_jl_invalid)

        expected_error = (
            "`unit_cell_jl` must be a Julia `UnitCell` object. "
            f"(unit_cell_jl={unit_cell_jl_invalid})."
        )
        assert expected_error in str(exception_info)

        # unit_cell_jl is not for a cubic unit cell
        unit_cell_jl_invalid = _JL.UnitCell(
            _JL.TetragonalLatticeConstants(1, 2), _JL.primitive
        )
        with pytest.raises(ValueError) as exception_info:
            CubicUnitCell.from_julia(unit_cell_jl_invalid)

        expected_error = (
            "`unit_cell_jl` must be a Julia `UnitCell` object for cubic "
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

        # --- Tests

        # centering = primitive
        unit_cell = CubicUnitCell(a, centering=Centering.PRIMITIVE)
        assert str(unit_cell) == f"CubicUnitCell(a={a},centering='primitive')"

        # centering = body_centered
        unit_cell = CubicUnitCell(a, centering=Centering.BODY_CENTERED)
        assert str(unit_cell) == f"CubicUnitCell(a={a},centering='body_centered')"

        # centering = face_centered
        unit_cell = CubicUnitCell(a, centering=Centering.FACE_CENTERED)
        assert str(unit_cell) == f"CubicUnitCell(a={a},centering='face_centered')"

    @staticmethod
    def test_eq():
        """
        Test `__eq__()`.
        """
        # --- Preparations

        # lattice constants
        a = 1

        # --- Tests

        # ------ types differ

        unit_cell_1 = CubicUnitCell(a)
        unit_cell_2 = TetragonalUnitCell(a, a + 1)
        assert unit_cell_1 != unit_cell_2

        # ------ lattice constants are the same

        # centerings are the same
        unit_cell_1 = CubicUnitCell(a)
        unit_cell_2 = CubicUnitCell(a, centering="primitive")
        assert unit_cell_1 == unit_cell_2

        # centerings are different
        unit_cell_1 = CubicUnitCell(a)
        unit_cell_2 = CubicUnitCell(a, centering="face_centered")
        assert unit_cell_1 != unit_cell_2

        # ------ lattice constants differ, centerings are the same

        # `a` values differ
        unit_cell_1 = CubicUnitCell(a)
        unit_cell_2 = CubicUnitCell(a + 1, centering="primitive")
        assert unit_cell_1 != unit_cell_2

    @staticmethod
    def test_isclose():
        """
        Test `isclose()`.
        """
        # --- Preparations

        # lattice constants
        a = 1

        # --- Tests

        # ------ types differ

        unit_cell_1 = CubicUnitCell(a)
        unit_cell_2 = TetragonalUnitCell(a, a + 1)
        assert not unit_cell_1.isclose(unit_cell_2)

        # ------ `a`

        # `a` values differ by less than tolerance, centerings are the same
        unit_cell_1 = CubicUnitCell(a)
        unit_cell_2 = CubicUnitCell(a + 0.1, centering="primitive")
        assert unit_cell_1.isclose(unit_cell_2, atol=0.2)

        # `a` values differ by less than tolerance, centerings differ
        unit_cell_1 = CubicUnitCell(a)
        unit_cell_2 = CubicUnitCell(a + 0.1, centering="face_centered")
        assert not unit_cell_1.isclose(unit_cell_2, atol=0.2)

        # `a` values differ by more than tolerance, centerings are the same
        unit_cell_1 = CubicUnitCell(a)
        unit_cell_2 = CubicUnitCell(a + 1, centering="primitive")
        assert not unit_cell_1.isclose(unit_cell_2, atol=0.2)
