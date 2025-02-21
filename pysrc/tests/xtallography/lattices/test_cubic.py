# Copyright (c) 2025 Velexi Corporation
#
# Use of this software is governed by the Business Source License included
# in the LICENSE file and at www.mariadb.com/bsl11.
#
# Change Date: Four years from the date the Licensed Work is published
#
# Change License: Mozilla Public License 2.0
"""
Unit tests for `xtallography.lattices.cubic` module
"""
# --- Imports

# Standard library
import unittest

# External packages
import pytest

# Local packages/modules
from xtallography import jl
from xtallography.lattices import LatticeSystem, Centering, CubicUnitCell


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
        assert jl.isa(unit_cell_jl, jl.CubicLatticeConstants)
