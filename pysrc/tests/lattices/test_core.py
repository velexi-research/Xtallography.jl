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
Unit tests for `xtallography.lattices.core` module
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
from xtallography.lattices import LatticeSystem, Centering, Lattice
from xtallography.lattices import UnitCell, CubicUnitCell

# Local packages/modules


# --- Test Suites


class test_xtallography_lattices_core(unittest.TestCase):
    """
    Test suite for the `xtallography.lattices.core` module
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
    def test_LatticeSystem():
        """
        Basic tests for `LatticeSystem` class.
        """
        # --- Tests

        # Get full list of values
        # TODO

    @staticmethod
    def test_Centering():
        """
        Basic tests for `Centering` class.
        """
        # --- Tests

        # Get full list of values
        # TODO

    def test_Centering_to_julia(self):
        """
        Test `Centering.to_julia()` method.
        """
        # --- Tests

        # primitive
        centering = Centering.PRIMITIVE
        centering_jl = centering.to_julia()
        assert self.jl.isa(centering_jl, self.jl.Primitive)

        # base_centered
        centering = Centering.BASE_CENTERED
        centering_jl = centering.to_julia()
        assert self.jl.isa(centering_jl, self.jl.BaseCentered)

        # body_centered
        centering = Centering.BODY_CENTERED
        centering_jl = centering.to_julia()
        assert self.jl.isa(centering_jl, self.jl.BodyCentered)

        # face_centered
        centering = Centering.FACE_CENTERED
        centering_jl = centering.to_julia()
        assert self.jl.isa(centering_jl, self.jl.FaceCentered)

    @staticmethod
    def test_Centering_from_julia():
        """
        Test `Centering.from_julia()` method.
        """
        # --- Tests

        # primitive
        centering_jl = _JL.primitive
        centering = Centering.from_julia(centering_jl)
        assert centering == Centering.PRIMITIVE

        # base_centered
        centering_jl = _JL.base_centered
        centering = Centering.from_julia(centering_jl)
        assert centering == Centering.BASE_CENTERED

        # body_centered
        centering_jl = _JL.body_centered
        centering = Centering.from_julia(centering_jl)
        assert centering == Centering.BODY_CENTERED

        # face_centered
        centering_jl = _JL.face_centered
        centering = Centering.from_julia(centering_jl)
        assert centering == Centering.FACE_CENTERED

    @staticmethod
    def test_Centering_from_julia_invalid_args():
        """
        Test argument checks for `Centering.from_julia()`.
        """
        # --- Tests

        # ------ `centering_jl` is not a Julia Centering object

        centering_jl_invalid = "not a Julia Centering object"
        with pytest.raises(ValueError) as exception_info:
            xtallography.lattices.Centering.from_julia(centering_jl_invalid)

        expected_error = (
            "`centering_jl` must be a Julia `Centering` object. "
            f"(centering_jl={centering_jl_invalid})."
        )
        assert expected_error in str(exception_info)

        # ------ `centering_jl` is an unsupported Julia Centering object

        _JL.seval("struct UnsupportedCentering <: Centering end")
        centering_jl_invalid = _JL.UnsupportedCentering()
        with pytest.raises(ValueError) as exception_info:
            xtallography.lattices.Centering.from_julia(centering_jl_invalid)

        expected_error = (
            f"Unsupported Centering type. (centering_jl={centering_jl_invalid})"
        )
        assert expected_error in str(exception_info)

    @staticmethod
    def test_Lattice():
        """
        Test `Lattice` type.
        """
        # --- Tests

        # ------ Valid field values

        # field values provided as LatticeSystem and Centering types
        try:
            Lattice(LatticeSystem.CUBIC, Centering.PRIMITIVE)
            assert True
        except Exception:
            pytest.fail("Lattice() constructor raised unexpected error")

        # field values provided as strings
        try:
            Lattice("cubic", "primitive")
            assert True
        except Exception:
            pytest.fail("Lattice() constructor raised unexpected error")

        # ------ Invalid field values

        # Nonsense `system` and `centering` values
        #
        # Note: Lattice() constructor doesn't validate values
        try:
            Lattice("abc", "def")
            assert True
        except Exception:
            pytest.fail("Lattice() constructor raised unexpected error")

        # `system` value has wrong type
        #
        # Note: Lattice() constructor doesn't validate value types
        try:
            Lattice(1, "def")
            assert True
        except Exception:
            pytest.fail("Lattice() constructor raised unexpected error")

        # `centering` value has wrong type
        #
        # Note: Lattice() constructor doesn't validate value types
        try:
            Lattice("abc", None)
            assert True
        except Exception:
            pytest.fail("Lattice() constructor raised unexpected error")

    @staticmethod
    def test_create_lattice():
        """
        Test `create_lattice()`.
        """
        # --- Tests

        # Arguments valid and lowercase
        lattice = xtallography.lattices.create_lattice("cubic", "primitive")

        assert lattice.lattice_system == LatticeSystem.CUBIC
        assert lattice.centering == Centering.PRIMITIVE

        # Arguments valid but not lowercase
        lattice = xtallography.lattices.create_lattice("cuBIc", "pRIMitIVe")

        assert lattice.lattice_system == LatticeSystem.CUBIC
        assert lattice.centering == Centering.PRIMITIVE

        # Arguments provided as keyword arguments
        lattice = xtallography.lattices.create_lattice(
            lattice_system="monoclinic",
            centering="body_centered",
        )

        assert lattice.lattice_system == LatticeSystem.MONOCLINIC
        assert lattice.centering == Centering.BODY_CENTERED

        # Arguments provided as LatticeSystem and Centering types
        lattice = xtallography.lattices.create_lattice(
            LatticeSystem.CUBIC, Centering.PRIMITIVE
        )

        assert lattice.lattice_system == LatticeSystem.CUBIC
        assert lattice.centering == Centering.PRIMITIVE

    @staticmethod
    def test_create_lattice_invalid_args():
        """
        Test argument checks for `create_lattice()`.
        """
        # --- Tests

        # default arguments
        with pytest.raises(ValueError) as exception_info:
            xtallography.lattices.create_lattice()

        expected_error = "`lattice_system` must be a string. (lattice_system=None)"
        assert expected_error in str(exception_info)

        # valid `system`, default `centering`
        with pytest.raises(ValueError) as exception_info:
            xtallography.lattices.create_lattice("cubic")

        expected_error = "`centering` must be a string. (centering=None)"
        assert expected_error in str(exception_info)

        # default `system`, valid `centering` (via keyword argument)
        with pytest.raises(ValueError) as exception_info:
            xtallography.lattices.create_lattice(centering="primitive")

        expected_error = "`lattice_system` must be a string. (lattice_system=None)"
        assert expected_error in str(exception_info)

        # `system` has wrong type
        with pytest.raises(ValueError) as exception_info:
            xtallography.lattices.create_lattice(1, "def")

        expected_error = "`lattice_system` must be a string. (lattice_system=1)"
        assert expected_error in str(exception_info)

        # `centering` has wrong type
        with pytest.raises(ValueError) as exception_info:
            xtallography.lattices.create_lattice("abc", None)

        expected_error = "`centering` must be a string. (centering=None)"
        assert expected_error in str(exception_info)

        # Nonsense `system` and `centering` values
        with pytest.raises(ValueError) as exception_info:
            xtallography.lattices.create_lattice("abc", "def")

        expected_error = (
            f"`lattice_system` must be one of "
            f"{xtallography.lattices.LatticeSystem.values()}. "
            "(lattice_system='abc')"
        )
        assert expected_error in str(exception_info)

        # Nonsense `centering` value
        with pytest.raises(ValueError) as exception_info:
            xtallography.lattices.create_lattice("hexagonal", "def")

        expected_error = (
            "`centering` must be one of "
            f"{xtallography.lattices.Centering.values()}. "
            "(centering='def')"
        )
        assert expected_error in str(exception_info)

    @staticmethod
    def test_standardize_lattice():
        """
        Test `standardize_lattice()`.
        """
        # --- Tests

        # Field values are valid and already standardized
        lattice = Lattice("cubic", "primitive")

        lattice = xtallography.lattices.standardize_lattice(lattice)

        assert isinstance(lattice.lattice_system, LatticeSystem)
        assert lattice.lattice_system == LatticeSystem.CUBIC

        assert isinstance(lattice.centering, Centering)
        assert lattice.centering == Centering.PRIMITIVE

        # Field values are valid but not standardized
        lattice = Lattice("cuBIc", "pRIMitIVe")

        lattice = xtallography.lattices.standardize_lattice(lattice)

        assert isinstance(lattice.lattice_system, LatticeSystem)
        assert lattice.lattice_system == LatticeSystem.CUBIC

        assert isinstance(lattice.centering, Centering)
        assert lattice.centering == Centering.PRIMITIVE

    @staticmethod
    def test_standardize_lattice_invalid_args():
        """
        Test argument checks for `standardize_lattice()`.
        """
        # --- Tests

        # `system` has wrong type
        with pytest.raises(ValueError) as exception_info:
            xtallography.lattices.standardize_lattice(Lattice(1, "def"))

        expected_error = (
            "`lattice.lattice_system` must be a string. (lattice.lattice_system=1)"
        )
        assert expected_error in str(exception_info)

        # `centering` has wrong type
        with pytest.raises(ValueError) as exception_info:
            xtallography.lattices.standardize_lattice(Lattice("abc", None))

        expected_error = (
            "`lattice.centering` must be a string. (lattice.centering=None)"
        )
        assert expected_error in str(exception_info)

        # Nonsense `system` and `centering` values
        with pytest.raises(ValueError) as exception_info:
            xtallography.lattices.standardize_lattice(Lattice("abc", "def"))

        expected_error = (
            "`lattice_system` must be one of "
            f"{xtallography.lattices.LatticeSystem.values()}. "
            "(lattice_system='abc')"
        )
        assert expected_error in str(exception_info)

        # Nonsense `centering` value
        with pytest.raises(ValueError) as exception_info:
            xtallography.lattices.standardize_lattice(Lattice("hexagonal", "def"))

        expected_error = (
            "`centering` must be one of "
            f"{xtallography.lattices.Centering.values()}. "
            "(centering='def')"
        )
        assert expected_error in str(exception_info)

    @staticmethod
    def test_is_bravais_lattice():
        """
        Test `is_bravais_lattice()`.
        """
        # --- Tests

        # --------- triclinic

        # (triclinic, primitive)
        lattice = Lattice("triclinic", "primitive")

        assert xtallography.lattices.is_bravais_lattice(lattice)

        # --------- monoclinic

        # (monoclinic, primitive)
        lattice = Lattice("monoclinic", "primitive")

        assert xtallography.lattices.is_bravais_lattice(lattice)

        # (monoclinic, body_centered)
        lattice = Lattice("monoclinic", "body_centered")

        assert xtallography.lattices.is_bravais_lattice(lattice)

        # (monoclinic, base_centered)
        lattice = Lattice("monoclinic", "base_centered")

        assert xtallography.lattices.is_bravais_lattice(lattice)

        # --------- orthorhombic

        # (orthorhombic, primitive)
        lattice = Lattice("orthorhombic", Centering.PRIMITIVE)

        assert xtallography.lattices.is_bravais_lattice(lattice)

        # (orthorhombic, body_centered)
        lattice = Lattice(LatticeSystem.ORTHORHOMBIC, "body_centered")

        assert xtallography.lattices.is_bravais_lattice(lattice)

        # (orthorhombic, base_centered)
        lattice = Lattice(LatticeSystem.ORTHORHOMBIC, Centering.BASE_CENTERED)

        assert xtallography.lattices.is_bravais_lattice(lattice)

        # (orthorhombic, face_centered)
        lattice = Lattice("orthorhombic", "face_centered")

        assert xtallography.lattices.is_bravais_lattice(lattice)

        # --------- tetragonal

        # (tetragonal, primitive)
        lattice = Lattice("tetragonal", Centering.PRIMITIVE)

        assert xtallography.lattices.is_bravais_lattice(lattice)

        # (tetragonal, body_centered)
        lattice = Lattice(LatticeSystem.TETRAGONAL, Centering.BODY_CENTERED)

        assert xtallography.lattices.is_bravais_lattice(lattice)

        # --------- rhomobohedral

        # (rhombohedral, primitive)
        lattice = Lattice("rhombohedral", Centering.PRIMITIVE)

        assert xtallography.lattices.is_bravais_lattice(lattice)

        # --------- hexagonal

        # (hexagonal, primitive)
        lattice = Lattice("hexagonal", "primitive")

        assert xtallography.lattices.is_bravais_lattice(lattice)

        # --------- cubic

        # (cubic, primitive)
        lattice = Lattice("cubic", "primitive")

        assert xtallography.lattices.is_bravais_lattice(lattice)

        # (cubic, body_centered)
        lattice = Lattice(LatticeSystem.CUBIC, "body_centered")

        assert xtallography.lattices.is_bravais_lattice(lattice)

        # (cubic, face_centered)
        lattice = Lattice("cubic", "face_centered")

        assert xtallography.lattices.is_bravais_lattice(lattice)

    @staticmethod
    def test_is_bravais_lattice_invalid_bravais_lattices():
        """
        Test `is_bravais_lattice()`.
        """
        # --- Tests

        # --------- triclinic

        # (triclinic, body_centered)
        lattice = Lattice(LatticeSystem.TRICLINIC, Centering.BODY_CENTERED)
        assert not xtallography.lattices.is_bravais_lattice(lattice)

        # (triclinic, base_centered)
        lattice = Lattice(LatticeSystem.TRICLINIC, Centering.BASE_CENTERED)
        assert not xtallography.lattices.is_bravais_lattice(lattice)

        # (triclinic, face_centered)
        lattice = Lattice(LatticeSystem.TRICLINIC, Centering.FACE_CENTERED)
        assert not xtallography.lattices.is_bravais_lattice(lattice)

        # --------- monoclinic

        # (monoclinic, face_centered)
        lattice = Lattice(LatticeSystem.MONOCLINIC, Centering.FACE_CENTERED)
        assert not xtallography.lattices.is_bravais_lattice(lattice)

        # --------- orthorhombic

        # all centerings are valid

        # --------- tetragonal

        # (tetragonal, base_centered)
        lattice = Lattice(LatticeSystem.TETRAGONAL, Centering.BASE_CENTERED)
        assert not xtallography.lattices.is_bravais_lattice(lattice)

        # (tetragonal, face_centered)
        lattice = Lattice(LatticeSystem.TETRAGONAL, Centering.FACE_CENTERED)
        assert not xtallography.lattices.is_bravais_lattice(lattice)

        # --------- rhomobohedral

        # (rhombohedral, body_centered)
        lattice = Lattice(LatticeSystem.RHOMBOHEDRAL, Centering.BODY_CENTERED)
        assert not xtallography.lattices.is_bravais_lattice(lattice)

        # (rhombohedral, base_centered)
        lattice = Lattice(LatticeSystem.RHOMBOHEDRAL, Centering.BASE_CENTERED)
        assert not xtallography.lattices.is_bravais_lattice(lattice)

        # (rhombohedral, face_centered)
        lattice = Lattice(LatticeSystem.RHOMBOHEDRAL, Centering.FACE_CENTERED)
        assert not xtallography.lattices.is_bravais_lattice(lattice)

        # --------- hexagonal

        # (hexagonal, body_centered)
        lattice = Lattice(LatticeSystem.HEXAGONAL, Centering.BODY_CENTERED)
        assert not xtallography.lattices.is_bravais_lattice(lattice)

        # (hexagonal, base_centered)
        lattice = Lattice(LatticeSystem.HEXAGONAL, Centering.BASE_CENTERED)
        assert not xtallography.lattices.is_bravais_lattice(lattice)

        # (hexagonal, face_centered)
        lattice = Lattice(LatticeSystem.HEXAGONAL, Centering.FACE_CENTERED)
        assert not xtallography.lattices.is_bravais_lattice(lattice)

        # --------- cubic

        # (cubic, base_centered)
        lattice = Lattice(LatticeSystem.CUBIC, Centering.BASE_CENTERED)
        assert not xtallography.lattices.is_bravais_lattice(lattice)

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
