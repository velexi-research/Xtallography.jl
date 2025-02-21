# Copyright (c) 2025 Velexi Corporation
#
# Use of this software is governed by the Business Source License included
# in the LICENSE file and at www.mariadb.com/bsl11.
#
# Change Date: Four years from the date the Licensed Work is published
#
# Change License: Mozilla Public License 2.0
"""
Unit tests for `xtallography.lattices.monoclinic` module
"""
# --- Imports

# Standard library
import math
import unittest

# External packages
import pytest

# Local packages/modules
from xtallography.lattices import jl
from xtallography.lattices import LatticeSystem, Centering, MonoclinicUnitCell


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

    @staticmethod
    def test_to_julia():
        """
        Test `to_julia()`.
        """
        # --- Preparations

        a = 1
        b = 2
        c = 3
        beta = 0.1

        unit_cell = MonoclinicUnitCell(a, b, c, beta)

        # --- Tests

        unit_cell_jl = unit_cell.to_julia()
        assert jl.isa(unit_cell_jl, jl.MonoclinicLatticeConstants)
