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
Unit tests for the `xtallography.lattices.core.UnitCell` class
"""
# --- Imports

# Standard library
import math
import unittest
import sys

# External packages
import juliacall
import pytest
import xtallography
from xtallography import _JL
from xtallography.lattices import LatticeSystem, Centering
from xtallography.lattices import (
    UnitCell,
    TriclinicUnitCell,
    MonoclinicUnitCell,
    OrthorhombicUnitCell,
    TetragonalUnitCell,
    RhombohedralUnitCell,
    HexagonalUnitCell,
    CubicUnitCell,
)

# Local packages/modules


# --- Test Suites


class test_xtallography_lattices_core_UnitCell(unittest.TestCase):
    """
    Test suite for the `xtallography.lattices.core.UnitCell` class
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
    def test_UnitCell_argument_checks():
        """
        Tests for `UnitCell()` constructor argument checks.
        """
        # --- Preparations

        # Define subclass that calls UnitCell constructor with string arguments
        class TestUnitCell(UnitCell):
            def __init__(self):
                super().__init__("cubic", "primitive")

            def to_julia(self):
                return None

            def from_julia(self):
                return None

            def __repr__(self):
                return "TestUnitCell()"

        # --- Tests

        try:
            unit_cell = TestUnitCell()
        except Exception:
            pytest.fail("Unexpected test failure.")

        assert isinstance(unit_cell.lattice_system, LatticeSystem)
        assert isinstance(unit_cell.centering, Centering)

    @staticmethod
    def test_isclose_default_args():
        """
        Test default argument cases for `UnitCell.isclose()`.
        """
        # --- Tests

        # ------ default `atol`

        # lattice constants differ by less than sqrt(eps)
        unit_cell_1 = CubicUnitCell(1)
        unit_cell_2 = CubicUnitCell(1 + 0.5 * math.sqrt(sys.float_info.epsilon))
        assert unit_cell_1.isclose(unit_cell_2)

        # lattice constants differ by more than sqrt(eps)
        unit_cell_1 = CubicUnitCell(1)
        unit_cell_2 = CubicUnitCell(1 - 2 * math.sqrt(sys.float_info.epsilon))
        assert not unit_cell_1.isclose(unit_cell_2)

        # ------ `rtol`

        # atol > 0, default rtol, lattice constants differ by less than atol
        unit_cell_1 = CubicUnitCell(1)
        unit_cell_2 = CubicUnitCell(1 + 0.05)
        assert unit_cell_1.isclose(unit_cell_2, atol=0.1)

        # atol > 0, default rtol, lattice constants differ by more than atol
        unit_cell_1 = CubicUnitCell(1)
        unit_cell_2 = CubicUnitCell(1 + 0.2)
        assert not unit_cell_1.isclose(unit_cell_2, atol=0.1)

    @staticmethod
    def test_isclose_invalid_args():
        """
        Test argument checks for `UnitCell.isclose()`.
        """
        # --- Tests

        # ------ `atol`

        # atol < 0
        unit_cell_1 = CubicUnitCell(1)
        unit_cell_2 = CubicUnitCell(1)
        atol_invalid = -0.1
        with pytest.raises(ValueError) as exception_info:
            unit_cell_1.isclose(unit_cell_2, atol=atol_invalid)

        expected_error = f"`atol` must be nonnegative. (atol={atol_invalid})"
        assert expected_error in str(exception_info)

        # ------ `rtol`

        # rtol < 0
        unit_cell_1 = CubicUnitCell(1)
        unit_cell_2 = CubicUnitCell(1)
        rtol_invalid = -0.1
        with pytest.raises(ValueError) as exception_info:
            unit_cell_1.isclose(unit_cell_2, rtol=rtol_invalid)

        expected_error = f"`rtol` must be nonnegative. (rtol={rtol_invalid})"
        assert expected_error in str(exception_info)

    @staticmethod
    def test_from_julia():
        """
        Test `UnitCell.from_julia()`.
        """
        # --- Tests

        # triclinic
        unit_cell_jl = TriclinicUnitCell(1, 2, 3, 0.1, 0.2, 0.3).to_julia()
        unit_cell = UnitCell.from_julia(unit_cell_jl)
        assert isinstance(unit_cell, TriclinicUnitCell)

        # monoclinic
        unit_cell_jl = MonoclinicUnitCell(1, 2, 3, 1.8).to_julia()
        unit_cell = UnitCell.from_julia(unit_cell_jl)
        assert isinstance(unit_cell, MonoclinicUnitCell)

        # orthorhombic
        unit_cell_jl = OrthorhombicUnitCell(1, 2, 3).to_julia()
        unit_cell = UnitCell.from_julia(unit_cell_jl)
        assert isinstance(unit_cell, OrthorhombicUnitCell)

        # tetragonal
        unit_cell_jl = TetragonalUnitCell(1, 2).to_julia()
        unit_cell = UnitCell.from_julia(unit_cell_jl)
        assert isinstance(unit_cell, TetragonalUnitCell)

        # rhombohedral
        unit_cell_jl = RhombohedralUnitCell(1, 1.8).to_julia()
        unit_cell = UnitCell.from_julia(unit_cell_jl)
        assert isinstance(unit_cell, RhombohedralUnitCell)

        # hexagonal
        unit_cell_jl = HexagonalUnitCell(1, 2).to_julia()
        unit_cell = UnitCell.from_julia(unit_cell_jl)
        assert isinstance(unit_cell, HexagonalUnitCell)

        # cubic
        unit_cell_jl = CubicUnitCell(1).to_julia()
        unit_cell = UnitCell.from_julia(unit_cell_jl)
        assert isinstance(unit_cell, CubicUnitCell)

    @staticmethod
    def test_from_julia_invalid_args():
        """
        Test argument checks for `UnitCell.from_julia()`.
        """
        # --- Tests

        # ------ `unit_cell_jl` is not a Julia LatticeSystem object

        unit_cell_jl_invalid = "not a Julia UnitCell object"
        with pytest.raises(ValueError) as exception_info:
            xtallography.lattices.UnitCell.from_julia(unit_cell_jl_invalid)

        expected_error = (
            "`unit_cell_jl` must be a Julia `UnitCell` object. "
            f"(unit_cell_jl={unit_cell_jl_invalid})."
        )
        assert expected_error in str(exception_info)

        # ------ `unit_cell_jl` is an unsupported Julia UnitCell object

        _JL.seval("import Xtallography: lattice_system")
        _JL.seval("struct UnsupportedLatticeSystem <: LatticeSystem end")
        _JL.seval("struct UnsupportedLatticeConstants <: LatticeConstants{Cubic} end")
        _JL.seval(
            "lattice_system(::UnsupportedLatticeConstants) = UnsupportedLatticeSystem()"
        )
        unit_cell_jl_invalid = _JL.UnitCell(
            _JL.UnsupportedLatticeConstants(), _JL.primitive
        )
        with pytest.raises(ValueError) as exception_info:
            xtallography.lattices.UnitCell.from_julia(unit_cell_jl_invalid)

        expected_error = (
            "Unsupported LatticeSystem type. "
            "(lattice_system_jl=UnsupportedLatticeSystem)"
        )
        assert expected_error in str(exception_info)
