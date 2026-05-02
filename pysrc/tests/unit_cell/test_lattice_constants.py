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
Unit tests for the `xtallography.unit_cell.lattice_constants` module
"""
# --- Imports

# Standard library
import dataclasses
import unittest

# External packages
from xtallography.unit_cell.lattice_constants import (
    LatticeConstants,
    TriclinicLatticeConstants,
    MonoclinicLatticeConstants,
    OrthorhombicLatticeConstants,
    TetragonalLatticeConstants,
    RhombohedralLatticeConstants,
    HexagonalLatticeConstants,
    CubicLatticeConstants,
)

# Local packages/modules


# --- Test Suites


class test_xtallography_unit_cell_LatticeConstants(unittest.TestCase):
    """
    Test suite for the `LatticeConstants` class
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
    def test_LatticeConstants():
        """
        Check subclasses of LatticeConstants
        """
        assert issubclass(TriclinicLatticeConstants, LatticeConstants)
        assert issubclass(MonoclinicLatticeConstants, LatticeConstants)
        assert issubclass(OrthorhombicLatticeConstants, LatticeConstants)
        assert issubclass(HexagonalLatticeConstants, LatticeConstants)
        assert issubclass(RhombohedralLatticeConstants, LatticeConstants)
        assert issubclass(TetragonalLatticeConstants, LatticeConstants)
        assert issubclass(CubicLatticeConstants, LatticeConstants)

    @staticmethod
    def test_TriclinicLatticeConstants():
        """
        Test TriclinicLatticeConstants
        """
        # --- Tests

        a = 1
        b = 2
        c = 3
        alpha = 0.1
        beta = 0.2
        gamma = 0.3

        lattice_constants = TriclinicLatticeConstants(a, b, c, alpha, beta, gamma)

        assert lattice_constants.a == a
        assert lattice_constants.b == b
        assert lattice_constants.c == c
        assert lattice_constants.alpha == alpha
        assert lattice_constants.beta == beta
        assert lattice_constants.gamma == gamma

        field_names = set(
            [field.name for field in dataclasses.fields(lattice_constants)]
        )
        assert field_names == set(["a", "b", "c", "alpha", "beta", "gamma"])

    @staticmethod
    def test_MonoclinicLatticeConstants():
        """
        Test MonoclinicLatticeConstants
        """
        # --- Tests

        a = 1
        b = 2
        c = 3
        beta = 0.2

        lattice_constants = MonoclinicLatticeConstants(a, b, c, beta)

        assert lattice_constants.a == a
        assert lattice_constants.b == b
        assert lattice_constants.c == c
        assert lattice_constants.beta == beta

        field_names = set(
            [field.name for field in dataclasses.fields(lattice_constants)]
        )
        assert field_names == set(["a", "b", "c", "beta"])

    @staticmethod
    def test_OrthorhombicLatticeConstants():
        """
        Test OrthorhombicLatticeConstants
        """
        # --- Tests

        a = 1
        b = 2
        c = 3

        lattice_constants = OrthorhombicLatticeConstants(a, b, c)

        assert lattice_constants.a == a
        assert lattice_constants.b == b
        assert lattice_constants.c == c

        field_names = set(
            [field.name for field in dataclasses.fields(lattice_constants)]
        )
        assert field_names == set(["a", "b", "c"])

    @staticmethod
    def test_HexagonalLatticeConstants():
        """
        Test HexagonalLatticeConstants
        """
        # --- Tests

        a = 1
        c = 3

        lattice_constants = HexagonalLatticeConstants(a, c)

        assert lattice_constants.a == a
        assert lattice_constants.c == c

        field_names = set(
            [field.name for field in dataclasses.fields(lattice_constants)]
        )
        assert field_names == set(["a", "c"])

    @staticmethod
    def test_TetragonalLatticeConstants():
        """
        Test TetragonalLatticeConstants
        """
        # --- Tests

        a = 1
        c = 3

        lattice_constants = TetragonalLatticeConstants(a, c)

        assert lattice_constants.a == a
        assert lattice_constants.c == c

        field_names = set(
            [field.name for field in dataclasses.fields(lattice_constants)]
        )
        assert field_names == set(["a", "c"])

    @staticmethod
    def test_RhombohedralLatticeConstants():
        """
        Test RhombohedralLatticeConstants
        """
        # --- Tests

        a = 1
        alpha = 3

        lattice_constants = RhombohedralLatticeConstants(a, alpha)

        assert lattice_constants.a == a
        assert lattice_constants.alpha == alpha

        field_names = set(
            [field.name for field in dataclasses.fields(lattice_constants)]
        )
        assert field_names == set(["a", "alpha"])

    @staticmethod
    def test_CubicLatticeConstants():
        """
        Test CubicLatticeConstants
        """
        # --- Tests

        a = 1

        lattice_constants = CubicLatticeConstants(a)

        assert lattice_constants.a == a

        field_names = set(
            [field.name for field in dataclasses.fields(lattice_constants)]
        )
        assert field_names == set(["a"])
