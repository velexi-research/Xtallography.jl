# Copyright (c) 2024 Velexi Corporation
#
# Use of this software is governed by the Business Source License included
# in the LICENSE file and at www.mariadb.com/bsl11.
#
# Change Date: Four years from the date the Licensed Work is published
#
# Change License: Mozilla Public License 2.0
"""
Unit tests for `xtallography.lattices.core` module
"""
# --- Imports

# Standard library
import copy
import math
import unittest

# External packages
import juliacall
import pytest
import xtallography
from xtallography.lattices import LatticeSystem, Centering, BravaisLattice

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
        self.jl = juliacall.newmodule("PyPxrdTest")
        self.jl.seval("using XtallographyUtils")

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
    def test_BravaisLattice():
        """
        Test `BravaisLattice` type.
        """
        # --- Tests

        # ------ Valid field values

        # field values provided as LatticeSystem and Centering types
        try:
            BravaisLattice(LatticeSystem.CUBIC, Centering.PRIMITIVE)
            assert True
        except Exception:
            pytest.fail("BravaisLattice() constructor raised unexpected error")

        # field values provided as strings
        try:
            BravaisLattice("cubic", "primitive")
            assert True
        except Exception:
            pytest.fail("BravaisLattice() constructor raised unexpected error")

        # ------ Invalid field values

        # Nonsense `system` and `centering` values
        #
        # Note: BravaisLattice() constructor doesn't validate values
        try:
            BravaisLattice("abc", "def")
            assert True
        except Exception:
            pytest.fail("BravaisLattice() constructor raised unexpected error")

        # `system` value has wrong type
        #
        # Note: BravaisLattice() constructor doesn't validate value types
        try:
            BravaisLattice(1, "def")
            assert True
        except Exception:
            pytest.fail("BravaisLattice() constructor raised unexpected error")

        # `centering` value has wrong type
        #
        # Note: BravaisLattice() constructor doesn't validate value types
        try:
            BravaisLattice("abc", None)
            assert True
        except Exception:
            pytest.fail("BravaisLattice() constructor raised unexpected error")

    @staticmethod
    def test_create_bravais_lattice():
        """
        Test `create_bravais_lattice()`.
        """
        # --- Tests

        # Arguments valid and lowercase
        lattice = xtallography.lattices.create_bravais_lattice("cubic", "primitive")

        assert lattice.lattice_system == LatticeSystem.CUBIC
        assert lattice.centering == Centering.PRIMITIVE

        # Arguments valid but not lowercase
        lattice = xtallography.lattices.create_bravais_lattice("cuBIc", "pRIMitIVe")

        assert lattice.lattice_system == LatticeSystem.CUBIC
        assert lattice.centering == Centering.PRIMITIVE

        # Arguments provided as keyword arguments
        lattice = xtallography.lattices.create_bravais_lattice(
            lattice_system="monoclinic",
            centering="body_centered",
        )

        assert lattice.lattice_system == LatticeSystem.MONOCLINIC
        assert lattice.centering == Centering.BODY_CENTERED

        # Arguments provided as LatticeSystem and Centering types
        lattice = xtallography.lattices.create_bravais_lattice(
            LatticeSystem.CUBIC, Centering.PRIMITIVE
        )

        assert lattice.lattice_system == LatticeSystem.CUBIC
        assert lattice.centering == Centering.PRIMITIVE

    @staticmethod
    def test_create_bravais_lattice_invalid_args():
        """
        Test argument checks for `create_bravais_lattice()`.
        """
        # --- Tests

        # default arguments
        with pytest.raises(ValueError) as exception_info:
            xtallography.lattices.create_bravais_lattice()

        expected_error = "`lattice_system` must be a string. (lattice_system=None)"
        assert expected_error in str(exception_info)

        # valid `system`, default `centering`
        with pytest.raises(ValueError) as exception_info:
            xtallography.lattices.create_bravais_lattice("cubic")

        expected_error = "`centering` must be a string. (centering=None)"
        assert expected_error in str(exception_info)

        # default `system`, valid `centering` (via keyword argument)
        with pytest.raises(ValueError) as exception_info:
            xtallography.lattices.create_bravais_lattice(centering="primitive")

        expected_error = "`lattice_system` must be a string. (lattice_system=None)"
        assert expected_error in str(exception_info)

        # `system` has wrong type
        with pytest.raises(ValueError) as exception_info:
            xtallography.lattices.create_bravais_lattice(1, "def")

        expected_error = "`lattice_system` must be a string. (lattice_system=1)"
        assert expected_error in str(exception_info)

        # `centering` has wrong type
        with pytest.raises(ValueError) as exception_info:
            xtallography.lattices.create_bravais_lattice("abc", None)

        expected_error = "`centering` must be a string. (centering=None)"
        assert expected_error in str(exception_info)

        # Nonsense `system` and `centering` values
        with pytest.raises(ValueError) as exception_info:
            xtallography.lattices.create_bravais_lattice("abc", "def")

        expected_error = (
            f"`lattice_system` must be one of "
            f"{xtallography.lattices.LatticeSystem.values()}. "
            "(lattice_system='abc')"
        )
        assert expected_error in str(exception_info)

        # Nonsense `centering` value
        with pytest.raises(ValueError) as exception_info:
            xtallography.lattices.create_bravais_lattice("hexagonal", "def")

        expected_error = (
            "`centering` must be one of "
            f"{xtallography.lattices.Centering.values()}. "
            "(centering='def')"
        )
        assert expected_error in str(exception_info)

    @staticmethod
    def test_standardize_bravais_lattice():
        """
        Test `standardize_bravais_lattice()`.
        """
        # --- Tests

        # Field values are valid and already standardized
        lattice = BravaisLattice("cubic", "primitive")

        lattice = xtallography.lattices.standardize_bravais_lattice(lattice)

        assert isinstance(lattice.lattice_system, LatticeSystem)
        assert lattice.lattice_system == LatticeSystem.CUBIC

        assert isinstance(lattice.centering, Centering)
        assert lattice.centering == Centering.PRIMITIVE

        # Field values are valid but not standardized
        lattice = BravaisLattice("cuBIc", "pRIMitIVe")

        lattice = xtallography.lattices.standardize_bravais_lattice(lattice)

        assert isinstance(lattice.lattice_system, LatticeSystem)
        assert lattice.lattice_system == LatticeSystem.CUBIC

        assert isinstance(lattice.centering, Centering)
        assert lattice.centering == Centering.PRIMITIVE

    @staticmethod
    def test_standardize_bravais_lattice_invalid_args():
        """
        Test argument checks for `standardize_bravais_lattice()`.
        """
        # --- Tests

        # `system` has wrong type
        with pytest.raises(ValueError) as exception_info:
            xtallography.lattices.standardize_bravais_lattice(BravaisLattice(1, "def"))

        expected_error = (
            "`lattice.lattice_system` must be a string. (lattice.lattice_system=1)"
        )
        assert expected_error in str(exception_info)

        # `centering` has wrong type
        with pytest.raises(ValueError) as exception_info:
            xtallography.lattices.standardize_bravais_lattice(
                BravaisLattice("abc", None)
            )

        expected_error = (
            "`lattice.centering` must be a string. (lattice.centering=None)"
        )
        assert expected_error in str(exception_info)

        # Nonsense `system` and `centering` values
        with pytest.raises(ValueError) as exception_info:
            xtallography.lattices.standardize_bravais_lattice(
                BravaisLattice("abc", "def")
            )

        expected_error = (
            "`lattice_system` must be one of "
            f"{xtallography.lattices.LatticeSystem.values()}. "
            "(lattice_system='abc')"
        )
        assert expected_error in str(exception_info)

        # Nonsense `centering` value
        with pytest.raises(ValueError) as exception_info:
            xtallography.lattices.standardize_bravais_lattice(
                BravaisLattice("hexagonal", "def")
            )

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
        lattice = BravaisLattice("triclinic", "primitive")

        assert xtallography.lattices.is_bravais_lattice(lattice)

        # --------- monoclinic

        # (monoclinic, primitive)
        lattice = BravaisLattice("monoclinic", "primitive")

        assert xtallography.lattices.is_bravais_lattice(lattice)

        # (monoclinic, body_centered)
        lattice = BravaisLattice("monoclinic", "body_centered")

        assert xtallography.lattices.is_bravais_lattice(lattice)

        # (monoclinic, base_centered)
        lattice = BravaisLattice("monoclinic", "base_centered")

        assert xtallography.lattices.is_bravais_lattice(lattice)

        # --------- orthorhombic

        # (orthorhombic, primitive)
        lattice = BravaisLattice("orthorhombic", Centering.PRIMITIVE)

        assert xtallography.lattices.is_bravais_lattice(lattice)

        # (orthorhombic, body_centered)
        lattice = BravaisLattice(LatticeSystem.ORTHORHOMBIC, "body_centered")

        assert xtallography.lattices.is_bravais_lattice(lattice)

        # (orthorhombic, base_centered)
        lattice = BravaisLattice(LatticeSystem.ORTHORHOMBIC, Centering.BASE_CENTERED)

        assert xtallography.lattices.is_bravais_lattice(lattice)

        # (orthorhombic, face_centered)
        lattice = BravaisLattice("orthorhombic", "face_centered")

        assert xtallography.lattices.is_bravais_lattice(lattice)

        # --------- tetragonal

        # (tetragonal, primitive)
        lattice = BravaisLattice("tetragonal", Centering.PRIMITIVE)

        assert xtallography.lattices.is_bravais_lattice(lattice)

        # (tetragonal, body_centered)
        lattice = BravaisLattice(LatticeSystem.TETRAGONAL, Centering.BODY_CENTERED)

        assert xtallography.lattices.is_bravais_lattice(lattice)

        # --------- rhomobohedral

        # (rhombohedral, primitive)
        lattice = BravaisLattice("rhombohedral", Centering.PRIMITIVE)

        assert xtallography.lattices.is_bravais_lattice(lattice)

        # --------- hexagonal

        # (hexagonal, primitive)
        lattice = BravaisLattice("hexagonal", "primitive")

        assert xtallography.lattices.is_bravais_lattice(lattice)

        # --------- cubic

        # (cubic, primitive)
        lattice = BravaisLattice("cubic", "primitive")

        assert xtallography.lattices.is_bravais_lattice(lattice)

        # (cubic, body_centered)
        lattice = BravaisLattice(LatticeSystem.CUBIC, "body_centered")

        assert xtallography.lattices.is_bravais_lattice(lattice)

        # (cubic, face_centered)
        lattice = BravaisLattice("cubic", "face_centered")

        assert xtallography.lattices.is_bravais_lattice(lattice)

    @staticmethod
    def test_is_bravais_lattice_invalid_bravais_lattices():
        """
        Test `is_bravais_lattice()`.
        """
        # --- Tests

        # --------- triclinic

        # (triclinic, body_centered)
        lattice = BravaisLattice(LatticeSystem.TRICLINIC, Centering.BODY_CENTERED)
        assert not xtallography.lattices.is_bravais_lattice(lattice)

        # (triclinic, base_centered)
        lattice = BravaisLattice(LatticeSystem.TRICLINIC, Centering.BASE_CENTERED)
        assert not xtallography.lattices.is_bravais_lattice(lattice)

        # (triclinic, face_centered)
        lattice = BravaisLattice(LatticeSystem.TRICLINIC, Centering.FACE_CENTERED)
        assert not xtallography.lattices.is_bravais_lattice(lattice)

        # --------- monoclinic

        # (monoclinic, face_centered)
        lattice = BravaisLattice(LatticeSystem.MONOCLINIC, Centering.FACE_CENTERED)
        assert not xtallography.lattices.is_bravais_lattice(lattice)

        # --------- orthorhombic

        # all centerings are valid

        # --------- tetragonal

        # (tetragonal, base_centered)
        lattice = BravaisLattice(LatticeSystem.TETRAGONAL, Centering.BASE_CENTERED)
        assert not xtallography.lattices.is_bravais_lattice(lattice)

        # (tetragonal, face_centered)
        lattice = BravaisLattice(LatticeSystem.TETRAGONAL, Centering.FACE_CENTERED)
        assert not xtallography.lattices.is_bravais_lattice(lattice)

        # --------- rhomobohedral

        # (rhombohedral, body_centered)
        lattice = BravaisLattice(LatticeSystem.RHOMBOHEDRAL, Centering.BODY_CENTERED)
        assert not xtallography.lattices.is_bravais_lattice(lattice)

        # (rhombohedral, base_centered)
        lattice = BravaisLattice(LatticeSystem.RHOMBOHEDRAL, Centering.BASE_CENTERED)
        assert not xtallography.lattices.is_bravais_lattice(lattice)

        # (rhombohedral, face_centered)
        lattice = BravaisLattice(LatticeSystem.RHOMBOHEDRAL, Centering.FACE_CENTERED)
        assert not xtallography.lattices.is_bravais_lattice(lattice)

        # --------- hexagonal

        # (hexagonal, body_centered)
        lattice = BravaisLattice(LatticeSystem.HEXAGONAL, Centering.BODY_CENTERED)
        assert not xtallography.lattices.is_bravais_lattice(lattice)

        # (hexagonal, base_centered)
        lattice = BravaisLattice(LatticeSystem.HEXAGONAL, Centering.BASE_CENTERED)
        assert not xtallography.lattices.is_bravais_lattice(lattice)

        # (hexagonal, face_centered)
        lattice = BravaisLattice(LatticeSystem.HEXAGONAL, Centering.FACE_CENTERED)
        assert not xtallography.lattices.is_bravais_lattice(lattice)

        # --------- cubic

        # (cubic, base_centered)
        lattice = BravaisLattice(LatticeSystem.CUBIC, Centering.BASE_CENTERED)
        assert not xtallography.lattices.is_bravais_lattice(lattice)

    @unittest.skip("BROKEN")
    @staticmethod
    def test_is_valid_unit_cell():
        """
        Test `is_valid_unit_cell()`.
        """
        # --- Tests

        # ------ triclinic

        lattice = xtallography.lattices.create_bravais_lattice("triclinic", "primitive")

        unit_cell = {"a": 1, "b": 2, "c": 3, "alpha": 80, "beta": 80, "gamma": 80}
        assert xtallography.lattices.is_valid_unit_cell(lattice, unit_cell)

        unit_cell = {"a": 1, "b": 2, "c": 3, "beta": 80}
        with pytest.raises(ValueError) as exception_info:
            xtallography.lattices.is_valid_unit_cell(lattice, unit_cell)

        expected_error = (
            "`unit_cell` is not compatible with `lattice`. "
            f"(lattice={lattice},unit_cell.keys={tuple(unit_cell.keys())})"
        )
        assert expected_error in str(exception_info)

        # ------ monoclinic

        lattice = xtallography.lattices.create_bravais_lattice(
            "monoclinic", "primitive"
        )

        unit_cell = {"a": 1, "b": 2, "c": 3, "beta": 80}
        assert xtallography.lattices.is_valid_unit_cell(lattice, unit_cell)

        unit_cell = {"a": 1, "b": 2, "c": 3, "alpha": 80, "beta": 80, "gamma": 80}
        with pytest.raises(ValueError) as exception_info:
            xtallography.lattices.is_valid_unit_cell(lattice, unit_cell)

        expected_error = (
            "`unit_cell` is not compatible with `lattice`. "
            f"(lattice={lattice},unit_cell.keys={tuple(unit_cell.keys())})"
        )
        assert expected_error in str(exception_info)

        # ------ orthorhombic

        lattice = xtallography.lattices.create_bravais_lattice(
            "orthorhombic", "primitive"
        )

        unit_cell = {"a": 1, "b": 2, "c": 3}
        assert xtallography.lattices.is_valid_unit_cell(lattice, unit_cell)

        unit_cell = {"a": 1, "b": 2, "c": 3, "alpha": 80, "beta": 80, "gamma": 80}
        with pytest.raises(ValueError) as exception_info:
            xtallography.lattices.is_valid_unit_cell(lattice, unit_cell)

        expected_error = (
            "`unit_cell` is not compatible with `lattice`. "
            f"(lattice={lattice},unit_cell.keys={tuple(unit_cell.keys())})"
        )
        assert expected_error in str(exception_info)

        # ------ tetragonal

        lattice = xtallography.lattices.create_bravais_lattice(
            "tetragonal", "body_centered"
        )

        unit_cell = {"a": 1, "c": 3}
        assert xtallography.lattices.is_valid_unit_cell(lattice, unit_cell)

        unit_cell = {"a": 1}
        with pytest.raises(ValueError) as exception_info:
            xtallography.lattices.is_valid_unit_cell(lattice, unit_cell)

        expected_error = (
            "`unit_cell` is not compatible with `lattice`. "
            f"(lattice={lattice},unit_cell.keys={tuple(unit_cell.keys())})"
        )
        assert expected_error in str(exception_info)

        # ------ rhombohedral

        lattice = xtallography.lattices.create_bravais_lattice(
            "rhombohedral", "primitive"
        )

        unit_cell = {"a": 1, "alpha": 45}
        assert xtallography.lattices.is_valid_unit_cell(lattice, unit_cell)

        unit_cell = {"a": 1, "c": 3}
        with pytest.raises(ValueError) as exception_info:
            xtallography.lattices.is_valid_unit_cell(lattice, unit_cell)

        expected_error = (
            "`unit_cell` is not compatible with `lattice`. "
            f"(lattice={lattice},unit_cell.keys={tuple(unit_cell.keys())})"
        )
        assert expected_error in str(exception_info)

        # ------ hexagonal

        lattice = xtallography.lattices.create_bravais_lattice("hexagonal", "primitive")

        unit_cell = {"a": 1, "c": 3}
        assert xtallography.lattices.is_valid_unit_cell(lattice, unit_cell)

        unit_cell = {"a": 1, "b": 2, "c": 3, "beta": 80}
        with pytest.raises(ValueError) as exception_info:
            xtallography.lattices.is_valid_unit_cell(lattice, unit_cell)

        expected_error = (
            "`unit_cell` is not compatible with `lattice`. "
            f"(lattice={lattice},unit_cell.keys={tuple(unit_cell.keys())})"
        )
        assert expected_error in str(exception_info)

        # ------ cubic

        lattice = xtallography.lattices.create_bravais_lattice("cubic", "face_centered")

        unit_cell = {"a": 1}
        assert xtallography.lattices.is_valid_unit_cell(lattice, unit_cell)

        unit_cell = {"a": 1, "c": 3}
        with pytest.raises(ValueError) as exception_info:
            xtallography.lattices.is_valid_unit_cell(lattice, unit_cell)

        expected_error = (
            "`unit_cell` is not compatible with `lattice`. "
            f"(lattice={lattice},unit_cell.keys={tuple(unit_cell.keys())})"
        )
        assert expected_error in str(exception_info)

        # ------ some `unit_cell` values are edge cases, allow_edge_case=True

        lattice = xtallography.lattices.create_bravais_lattice("hexagonal", "primitive")
        unit_cell = {"a": 1, "c": 3}

        # some `unit_cell` values are 0
        test_unit_cell = copy.copy(unit_cell)
        test_unit_cell["a"] = 0
        try:
            assert xtallography.lattices.is_valid_unit_cell(
                lattice, test_unit_cell, allow_edge_cases=True
            )
        except Exception:
            pytest.fail("BravaisLattice() constructor raised unexpected error")

        assert expected_error in str(exception_info)

        # some `unit_cell` values are inf
        test_unit_cell = copy.copy(unit_cell)
        test_unit_cell["a"] = math.inf
        try:
            assert xtallography.lattices.is_valid_unit_cell(
                lattice, test_unit_cell, allow_edge_cases=True
            )
        except Exception:
            pytest.fail("BravaisLattice() constructor raised unexpected error")

    @unittest.skip("BROKEN")
    @staticmethod
    def test_is_valid_unit_cell_invalid_args():
        """
        Test argument checks for `is_valid_unit_cell()`.
        """
        # --- Preparations

        # valid arguments
        lattice = xtallography.lattices.create_bravais_lattice(
            "orthorhombic", "body_centered"
        )
        unit_cell = {"a": 3, "b": 4, "c": 5}

        # --- Tests

        # `lattice` is not a `BravaisLattice` type
        invalid_lattice = "not a lattice"
        with pytest.raises(ValueError) as exception_info:
            xtallography.lattices.is_valid_unit_cell(invalid_lattice, unit_cell)

        expected_error = (
            f"`lattice` must be a `BravaisLattice` object. (lattice={invalid_lattice})"
        )
        assert expected_error in str(exception_info)

        # `lattice` is not a valid Bravais lattice
        invalid_lattice = xtallography.lattices.create_bravais_lattice(
            "hexagonal", "body_centered"
        )
        with pytest.raises(ValueError) as exception_info:
            xtallography.lattices.is_valid_unit_cell(invalid_lattice, unit_cell)

        expected_error = (
            f"`lattice` must be a Bravais lattice. (lattice={invalid_lattice})"
        )
        assert expected_error in str(exception_info)

        # `unit_cell` is not a dict
        invalid_unit_cell = "not a dict"
        with pytest.raises(ValueError) as exception_info:
            xtallography.lattices.is_valid_unit_cell(lattice, invalid_unit_cell)

        expected_error = f"`unit_cell` must be a dict. (unit_cell={invalid_unit_cell})"
        assert expected_error in str(exception_info)

        # `unit_cell` values are not numbers
        invalid_unit_cell = copy.copy(unit_cell)
        invalid_unit_cell["a"] = "not a number"
        with pytest.raises(ValueError) as exception_info:
            xtallography.lattices.is_valid_unit_cell(lattice, invalid_unit_cell)

        expected_error = (
            "`unit_cell['a']` must be a number. (unit_cell['a']=not a number)"
        )

        assert expected_error in str(exception_info)

        # some `unit_cell` values are 0, allow_edge_cases = False
        invalid_unit_cell = copy.copy(unit_cell)
        invalid_unit_cell["a"] = 0
        with pytest.raises(ValueError) as exception_info:
            xtallography.lattices.is_valid_unit_cell(
                lattice, invalid_unit_cell, allow_edge_cases=False
            )

        expected_error = "`unit_cell['a']` must be positive. (unit_cell['a']=0)"

        assert expected_error in str(exception_info)

        # some `unit_cell` values are negative, allow_edge_cases = True
        invalid_unit_cell = copy.copy(unit_cell)
        invalid_unit_cell["a"] = -10
        with pytest.raises(ValueError) as exception_info:
            xtallography.lattices.is_valid_unit_cell(
                lattice, invalid_unit_cell, allow_edge_cases=True
            )

        expected_error = "`unit_cell['a']` must be nonnegative. (unit_cell['a']=-10)"

        assert expected_error in str(exception_info)

        # some `unit_cell` values are negative, allow_edge_cases = False
        invalid_unit_cell = copy.copy(unit_cell)
        invalid_unit_cell["a"] = -10
        with pytest.raises(ValueError) as exception_info:
            xtallography.lattices.is_valid_unit_cell(
                lattice, invalid_unit_cell, allow_edge_cases=False
            )

        expected_error = "`unit_cell['a']` must be positive. (unit_cell['a']=-10)"

        assert expected_error in str(exception_info)

        # some `unit_cell` values are 0, allow_edge_cases = False
        invalid_unit_cell = copy.copy(unit_cell)
        invalid_unit_cell["a"] = math.inf
        with pytest.raises(ValueError) as exception_info:
            xtallography.lattices.is_valid_unit_cell(
                lattice, invalid_unit_cell, allow_edge_cases=False
            )

        expected_error = "`unit_cell['a']` must be finite. (unit_cell['a']=inf)"

        assert expected_error in str(exception_info)
