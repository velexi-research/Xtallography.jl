# Copyright (c) 2025 Velexi Corporation
#
# Use of this software is governed by the Business Source License included
# in the LICENSE file and at www.mariadb.com/bsl11.
#
# Change Date: Four years from the date the Licensed Work is published
#
# Change License: Mozilla Public License 2.0
"""
Unit tests for `xtallography.lattices.hexagonal` module
"""
# --- Imports

# Standard library
import unittest

# External packages
import pytest

# Local packages/modules
from xtallography import jl
from xtallography.lattices import LatticeSystem, Centering, HexagonalUnitCell


# --- Test Suites


class test_xtallography_lattice_hexagonal(unittest.TestCase):
    """
    Test suite for the `xtallography.lattice.hexagonal` module
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
        c = 3

        # --- Tests

        unit_cell = HexagonalUnitCell(a, c)
        assert unit_cell.a == a
        assert unit_cell.c == c
        assert unit_cell.lattice_system == LatticeSystem.HEXAGONAL
        assert unit_cell.centering == Centering.PRIMITIVE

    @staticmethod
    def test_init_invalid_args():
        """
        Test `__init__()`: invalid arguments.
        """
        # --- Preparations

        a = 1
        c = 3

        # --- Tests

        # ------ Invalid `a`

        # a < 0
        invalid_value = -1
        with pytest.raises(ValueError) as exception_info:
            HexagonalUnitCell(invalid_value, c)

        expected_error = f"`a` must be positive. (a={invalid_value})"
        assert expected_error in str(exception_info)

        # a = 0
        invalid_value = 0
        with pytest.raises(ValueError) as exception_info:
            HexagonalUnitCell(invalid_value, c)

        expected_error = f"`a` must be positive. (a={invalid_value})"
        assert expected_error in str(exception_info)

        # c < 0
        invalid_value = -1
        with pytest.raises(ValueError) as exception_info:
            HexagonalUnitCell(a, invalid_value)

        expected_error = f"`c` must be positive. (c={invalid_value})"
        assert expected_error in str(exception_info)

        # c = 0
        invalid_value = 0
        with pytest.raises(ValueError) as exception_info:
            HexagonalUnitCell(a, invalid_value)

        expected_error = f"`c` must be positive. (c={invalid_value})"
        assert expected_error in str(exception_info)

    @staticmethod
    def test_to_julia():
        """
        Test `to_julia()`.
        """
        # --- Preparations

        a = 1
        c = 4

        unit_cell = HexagonalUnitCell(a, c)

        # --- Tests

        unit_cell_jl = unit_cell.to_julia()
        assert jl.isa(unit_cell_jl, jl.HexagonalLatticeConstants)
