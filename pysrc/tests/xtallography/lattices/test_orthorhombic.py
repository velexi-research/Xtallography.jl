# Copyright (c) 2025 Velexi Corporation
#
# Use of this software is governed by the Business Source License included
# in the LICENSE file and at www.mariadb.com/bsl11.
#
# Change Date: Four years from the date the Licensed Work is published
#
# Change License: Mozilla Public License 2.0
"""
Unit tests for `xtallography.lattices.orthorhombic` module
"""
# --- Imports

# Standard library
import unittest

# External packages
import pytest

# Local packages/modules
from xtallography.lattices import jl
from xtallography.lattices import LatticeSystem, Centering, OrthorhombicUnitCell


# --- Test Suites


class test_xtallography_lattice_orthorhombic(unittest.TestCase):
    """
    Test suite for the `xtallography.lattice.orthorhombic` module
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
        b = 2
        c = 3

        # --- Tests

        # default centering
        unit_cell = OrthorhombicUnitCell(a, b, c)
        assert unit_cell.a == a
        assert unit_cell.b == b
        assert unit_cell.c == c
        assert unit_cell.lattice_system == LatticeSystem.ORTHORHOMBIC
        assert unit_cell.centering == Centering.PRIMITIVE

        # centering = primitive
        unit_cell = OrthorhombicUnitCell(a, b, c, centering=Centering.PRIMITIVE)
        assert unit_cell.a == a
        assert unit_cell.b == b
        assert unit_cell.c == c
        assert unit_cell.lattice_system == LatticeSystem.ORTHORHOMBIC
        assert unit_cell.centering == Centering.PRIMITIVE

        # centering = body-centered
        unit_cell = OrthorhombicUnitCell(a, b, c, centering=Centering.BODY_CENTERED)
        assert unit_cell.a == a
        assert unit_cell.b == b
        assert unit_cell.c == c
        assert unit_cell.lattice_system == LatticeSystem.ORTHORHOMBIC
        assert unit_cell.centering == Centering.BODY_CENTERED

        # centering = base-centered
        unit_cell = OrthorhombicUnitCell(a, b, c, centering=Centering.BASE_CENTERED)
        assert unit_cell.a == a
        assert unit_cell.b == b
        assert unit_cell.c == c
        assert unit_cell.lattice_system == LatticeSystem.ORTHORHOMBIC
        assert unit_cell.centering == Centering.BASE_CENTERED

        # centering = face-centered
        unit_cell = OrthorhombicUnitCell(a, b, c, centering=Centering.FACE_CENTERED)
        assert unit_cell.a == a
        assert unit_cell.b == b
        assert unit_cell.c == c
        assert unit_cell.lattice_system == LatticeSystem.ORTHORHOMBIC
        assert unit_cell.centering == Centering.FACE_CENTERED

        # centering argument is a str
        unit_cell = OrthorhombicUnitCell(a, b, c, centering="primitive")

        assert unit_cell.a == a
        assert unit_cell.b == b
        assert unit_cell.c == c
        assert isinstance(unit_cell.lattice_system, LatticeSystem)
        assert unit_cell.lattice_system == LatticeSystem.ORTHORHOMBIC
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

        # --- Tests

        # ------ Invalid `a`

        # a < 0
        invalid_value = -1
        with pytest.raises(ValueError) as exception_info:
            OrthorhombicUnitCell(invalid_value, b, c)

        expected_error = f"`a` must be positive. (a={invalid_value})"
        assert expected_error in str(exception_info)

        # a = 0
        invalid_value = 0
        with pytest.raises(ValueError) as exception_info:
            OrthorhombicUnitCell(invalid_value, b, c)

        expected_error = f"`a` must be positive. (a={invalid_value})"
        assert expected_error in str(exception_info)

        # ------ Invalid `b`

        # b < 0
        invalid_value = -1
        with pytest.raises(ValueError) as exception_info:
            OrthorhombicUnitCell(a, invalid_value, c)

        expected_error = f"`b` must be positive. (b={invalid_value})"
        assert expected_error in str(exception_info)

        # b = 0
        invalid_value = 0
        with pytest.raises(ValueError) as exception_info:
            OrthorhombicUnitCell(a, invalid_value, c)

        expected_error = f"`b` must be positive. (b={invalid_value})"
        assert expected_error in str(exception_info)

        # ------ Invalid `c`

        # c < 0
        invalid_value = -1
        with pytest.raises(ValueError) as exception_info:
            OrthorhombicUnitCell(a, b, invalid_value)

        expected_error = f"`c` must be positive. (c={invalid_value})"
        assert expected_error in str(exception_info)

        # c = 0
        invalid_value = 0
        with pytest.raises(ValueError) as exception_info:
            OrthorhombicUnitCell(a, b, invalid_value)

        expected_error = f"`c` must be positive. (c={invalid_value})"
        assert expected_error in str(exception_info)

    @staticmethod
    def test_to_julia():
        """
        Test `to_julia()`.
        """
        # --- Preparations

        a = 1
        b = 2
        c = 3

        unit_cell = OrthorhombicUnitCell(a, b, c)

        # --- Tests

        unit_cell_jl = unit_cell.to_julia()
        assert jl.isa(unit_cell_jl, jl.OrthorhombicLatticeConstants)
