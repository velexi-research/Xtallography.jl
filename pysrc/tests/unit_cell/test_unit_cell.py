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
Unit tests for the `xtallography.unit_cell.unit_cell` module
"""
# --- Imports

# Standard library
import math
import unittest
import sys

# External packages
import juliacall
import pytest
from xtallography.symmetry import Centering, LatticeSystem
from xtallography.unit_cell import (
    LatticeConstants,
    UnitCell,
    UnitCellSymmetry,
    TriclinicLatticeConstants,
    TriclinicUnitCell,
    MonoclinicUnitCell,
    OrthorhombicUnitCell,
    TetragonalUnitCell,
    RhombohedralUnitCell,
    HexagonalUnitCell,
    CubicLatticeConstants,
    CubicUnitCell,
    PRIMITIVE_UNIT_CELL_SYMMETRY,
)

# Local packages/modules


# --- Test Suites

# Lattice-specific tests for the following methods are located in subclass unit test files
#
# __init__()
# to_julia()
# from_julia()
# __repr__()
# __eq__()
# isclose()


class test_xtallography_unit_cell_UnitCell(unittest.TestCase):
    """
    Test suite for the `UnitCell` class
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
    def test_init_arg_checks_lattice_system():
        """
        Test argument checks for `UnitCell.__init__()`: lattice_system
        """
        # --- Preparations

        # define concrete subclass of UnitCell for testing purposes
        class TestUnitCell(UnitCell):
            def __init__(
                self, lattice_system: LatticeSystem, lattice_constants: LatticeConstants
            ):
                super().__init__(lattice_system, lattice_constants)

            def to_julia(self):
                return None

            def __repr__(self):
                return ""

        # --- Tests

        # ------ lattice_system is str

        # correct capitalization
        lattice_system = "cubic"
        lattice_constants = CubicLatticeConstants(1)

        unit_cell = TestUnitCell(lattice_system, lattice_constants)
        assert unit_cell.lattice_system == LatticeSystem.CUBIC

        # incorrect capitalization
        lattice_system = "Cubic"
        lattice_constants = CubicLatticeConstants(1)

        unit_cell = TestUnitCell(lattice_system, lattice_constants)
        assert unit_cell.lattice_system == LatticeSystem.CUBIC

    @staticmethod
    def test_init_invalid_args():
        """
        Test argument checks for `UnitCell.__init__()`.
        """
        # --- Preparations

        # define concrete subclass of UnitCell for testing purposes
        class TestUnitCell(UnitCell):
            def __init__(
                self,
                lattice_system: LatticeSystem,
                lattice_constants: LatticeConstants,
                symmetry: UnitCellSymmetry = PRIMITIVE_UNIT_CELL_SYMMETRY,
            ):
                super().__init__(lattice_system, lattice_constants, symmetry=symmetry)

            def to_julia(self):
                return None

            def __repr__(self):
                return ""

        # --- Tests

        # ------ type(lattice_system) not LatticeSystem or str

        lattice_system = 4
        lattice_constants = CubicLatticeConstants(1)

        with pytest.raises(ValueError) as exception_info:
            TestUnitCell(lattice_system, lattice_constants)

        expected_error = (
            "`lattice_system` must be a `LatticeSystem` or `str`."
            f"(lattice_system={lattice_system})"
        )
        assert expected_error in str(exception_info)

        # ------ lattice_system is not a valid LatticeSystem

        lattice_system = "triangle"
        lattice_constants = CubicLatticeConstants(1)

        with pytest.raises(ValueError) as exception_info:
            TestUnitCell(lattice_system, lattice_constants)

        expected_error = "'triangle' is not a valid LatticeSystem"
        assert expected_error in str(exception_info)

        # ------ lattice_constants incompatible with lattice_system

        # triclinic
        lattice_system = LatticeSystem.TRICLINIC
        lattice_constants = CubicLatticeConstants(1)

        with pytest.raises(ValueError) as exception_info:
            TestUnitCell(lattice_system, lattice_constants)

        expected_error = (
            "`lattice_constants` for triclinic unit cell must be a "
            "TriclinicLatticeConstants object. "
            f"(lattice_constants={lattice_constants})"
        )
        assert expected_error in str(exception_info)

        # monoclinic
        lattice_system = LatticeSystem.MONOCLINIC
        lattice_constants = CubicLatticeConstants(1)

        with pytest.raises(ValueError) as exception_info:
            TestUnitCell(lattice_system, lattice_constants)

        expected_error = (
            "`lattice_constants` for monoclinic unit cell must be a "
            "MonoclinicLatticeConstants object. "
            f"(lattice_constants={lattice_constants})"
        )
        assert expected_error in str(exception_info)

        # orthorhombic
        lattice_system = LatticeSystem.ORTHORHOMBIC
        lattice_constants = CubicLatticeConstants(1)

        with pytest.raises(ValueError) as exception_info:
            TestUnitCell(lattice_system, lattice_constants)

        expected_error = (
            "`lattice_constants` for orthorhombic unit cell must be a "
            "OrthorhombicLatticeConstants object. "
            f"(lattice_constants={lattice_constants})"
        )
        assert expected_error in str(exception_info)

        # hexagonal
        lattice_system = LatticeSystem.HEXAGONAL
        lattice_constants = CubicLatticeConstants(1)

        with pytest.raises(ValueError) as exception_info:
            TestUnitCell(lattice_system, lattice_constants)

        expected_error = (
            "`lattice_constants` for hexagonal unit cell must be a "
            "HexagonalLatticeConstants object. "
            f"(lattice_constants={lattice_constants})"
        )
        assert expected_error in str(exception_info)

        # rhombohedral
        lattice_system = LatticeSystem.RHOMBOHEDRAL
        lattice_constants = CubicLatticeConstants(1)

        with pytest.raises(ValueError) as exception_info:
            TestUnitCell(lattice_system, lattice_constants)

        expected_error = (
            "`lattice_constants` for rhombohedral unit cell must be a "
            "RhombohedralLatticeConstants object. "
            f"(lattice_constants={lattice_constants})"
        )
        assert expected_error in str(exception_info)

        # tetragonal
        lattice_system = LatticeSystem.TETRAGONAL
        lattice_constants = CubicLatticeConstants(1)

        with pytest.raises(ValueError) as exception_info:
            TestUnitCell(lattice_system, lattice_constants)

        expected_error = (
            "`lattice_constants` for tetragonal unit cell must be a "
            "TetragonalLatticeConstants object. "
            f"(lattice_constants={lattice_constants})"
        )
        assert expected_error in str(exception_info)

        # cubic
        lattice_system = LatticeSystem.CUBIC
        lattice_constants = TriclinicLatticeConstants(1, 2, 3, 0.1, 0.2, 0.3)

        with pytest.raises(ValueError) as exception_info:
            TestUnitCell(lattice_system, lattice_constants)

        expected_error = (
            "`lattice_constants` for cubic unit cell must be a "
            "CubicLatticeConstants object. "
            f"(lattice_constants={lattice_constants})"
        )
        assert expected_error in str(exception_info)

        # ------ symmetry is not a UnitCellSymmetry object

        lattice_system = LatticeSystem.CUBIC
        lattice_constants = CubicLatticeConstants(1)
        symmetry = 3

        with pytest.raises(ValueError) as exception_info:
            TestUnitCell(lattice_system, lattice_constants, symmetry=symmetry)

        expected_error = (
            f"`symmetry` must be a `UnitCellSymmetry` object. (symmetry={symmetry})"
        )
        assert expected_error in str(exception_info)

    @staticmethod
    def test_fields_and_properties():
        """
        Test fields and properties.
        """
        # --- Preparations

        a = 1
        b = 2
        c = 3
        alpha = 0.1
        beta = 0.2
        gamma = 0.3

        # --- Tests

        unit_cell = TriclinicUnitCell(a, b, c, alpha, beta, gamma)

        # fields
        assert unit_cell.lattice_constants == TriclinicLatticeConstants(
            a, b, c, alpha, beta, gamma
        )
        assert unit_cell.lattice_system == LatticeSystem.TRICLINIC
        assert unit_cell.symmetry == UnitCellSymmetry(
            centering=Centering.PRIMITIVE, symmetry_elements=set()
        )

        # properties
        assert unit_cell.centering == Centering.PRIMITIVE
        assert unit_cell.symmetry_elements == set()

    @staticmethod
    def test_isclose():
        """
        Test `UnitCell.isclose()` special cases
        """
        # --- Preparations

        # define concrete subclass of UnitCell for testing purposes
        class TestUnitCell(UnitCell):
            def __init__(
                self, lattice_system: LatticeSystem, lattice_constants: LatticeConstants
            ):
                super().__init__(lattice_system, lattice_constants)

            def to_julia(self):
                return None

            def __repr__(self):
                return ""

        # --- Tests

        # ------ lattice systems do not match

        lattice_system_1 = LatticeSystem.CUBIC
        lattice_constants_1 = CubicLatticeConstants(1)
        unit_cell_1 = TestUnitCell(lattice_system_1, lattice_constants_1)

        lattice_system_2 = LatticeSystem.TRICLINIC
        lattice_constants_2 = TriclinicLatticeConstants(1, 2, 3, 0.1, 0.2, 0.3)
        unit_cell_2 = TestUnitCell(lattice_system_2, lattice_constants_2)

        assert not unit_cell_1.isclose(unit_cell_2)

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
        unit_cell_ref = TriclinicUnitCell(1, 2, 3, 0.1, 0.2, 0.3)
        unit_cell_jl = unit_cell_ref.to_julia()
        unit_cell = UnitCell.from_julia(unit_cell_jl)
        assert isinstance(unit_cell, TriclinicUnitCell)
        assert unit_cell == unit_cell_ref
        assert unit_cell is not unit_cell_ref

        # monoclinic
        unit_cell_ref = MonoclinicUnitCell(1, 2, 3, 1.8)
        unit_cell_jl = unit_cell_ref.to_julia()
        unit_cell = UnitCell.from_julia(unit_cell_jl)
        assert isinstance(unit_cell, MonoclinicUnitCell)
        assert unit_cell == unit_cell_ref
        assert unit_cell is not unit_cell_ref

        # orthorhombic
        unit_cell_ref = OrthorhombicUnitCell(1, 2, 3)
        unit_cell_jl = unit_cell_ref.to_julia()
        unit_cell = UnitCell.from_julia(unit_cell_jl)
        assert isinstance(unit_cell, OrthorhombicUnitCell)
        assert unit_cell == unit_cell_ref
        assert unit_cell is not unit_cell_ref

        # tetragonal
        unit_cell_ref = TetragonalUnitCell(1, 2)
        unit_cell_jl = unit_cell_ref.to_julia()
        unit_cell = UnitCell.from_julia(unit_cell_jl)
        assert isinstance(unit_cell, TetragonalUnitCell)
        assert unit_cell == unit_cell_ref
        assert unit_cell is not unit_cell_ref

        # rhombohedral
        unit_cell_ref = RhombohedralUnitCell(1, 1.8)
        unit_cell_jl = unit_cell_ref.to_julia()
        unit_cell = UnitCell.from_julia(unit_cell_jl)
        assert isinstance(unit_cell, RhombohedralUnitCell)
        assert unit_cell == unit_cell_ref
        assert unit_cell is not unit_cell_ref

        # hexagonal
        unit_cell_ref = HexagonalUnitCell(1, 2)
        unit_cell_jl = unit_cell_ref.to_julia()
        unit_cell = UnitCell.from_julia(unit_cell_jl)
        assert isinstance(unit_cell, HexagonalUnitCell)
        assert unit_cell == unit_cell_ref
        assert unit_cell is not unit_cell_ref

        # cubic
        unit_cell_ref = CubicUnitCell(1)
        unit_cell_jl = unit_cell_ref.to_julia()
        unit_cell = UnitCell.from_julia(unit_cell_jl)
        assert isinstance(unit_cell, CubicUnitCell)
        assert unit_cell == unit_cell_ref
        assert unit_cell is not unit_cell_ref

    def test_from_julia_invalid_args(self):
        """
        Test argument checks for `UnitCell.from_julia()`.
        """
        # --- Tests

        # ------ `unit_cell_jl` is not a Julia UnitCell object

        unit_cell_jl_invalid = "not a Julia UnitCell object"
        with pytest.raises(ValueError) as exception_info:
            UnitCell.from_julia(unit_cell_jl_invalid)

        expected_error = (
            "`unit_cell_jl` must be a Julia `UnitCell` object. "
            f"(unit_cell_jl={unit_cell_jl_invalid})."
        )
        assert expected_error in str(exception_info)

        # ------ `unit_cell_jl` is a Julia UnitCell object for a custom lattice system

        self.jl.seval("import Xtallography: lattice_system")
        self.jl.seval("struct UnsupportedLatticeSystem <: LatticeSystem end")
        self.jl.seval("UnsupportedUnitCell = UnitCell{UnsupportedLatticeSystem}")
        self.jl.seval(
            "lattice_system(::UnsupportedUnitCell) = UnsupportedLatticeSystem()"
        )
        lattice_constants_jl = self.jl.seval("(a=1,)")
        unit_cell_jl_invalid = self.jl.UnsupportedUnitCell(lattice_constants_jl)

        with pytest.raises(ValueError) as exception_info:
            UnitCell.from_julia(unit_cell_jl_invalid)

        expected_error = (
            "Unsupported LatticeSystem type. "
            "(lattice_system_jl=UnsupportedLatticeSystem)"
        )
        assert expected_error in str(exception_info)
