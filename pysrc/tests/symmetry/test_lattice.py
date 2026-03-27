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
Unit tests for `xtallography.symmetry.lattice` module
"""
# --- Imports

# Standard library
import unittest

# External packages
import juliacall
import pytest
from xtallography.symmetry import Centering, Lattice, LatticeSystem

# Local packages/modules


# --- Test Suites


class test_xtallography_symmetry_lattice_Lattice(unittest.TestCase):
    """
    Test suite for the `Lattice` class
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
    def test_init():
        """
        Test `__init__()`.
        """
        # --- Tests

        # field values provided as LatticeSystem and Centering objects
        lattice = Lattice(LatticeSystem.CUBIC, Centering.PRIMITIVE)
        print(lattice)
        assert lattice.lattice_system == LatticeSystem.CUBIC
        assert lattice.centering == Centering.PRIMITIVE

        # field values provided as strings
        lattice = Lattice("cubic", "primitive")
        assert lattice.lattice_system == LatticeSystem.CUBIC
        assert lattice.centering == Centering.PRIMITIVE

        # field values provided as strings but not lowercase
        lattice = Lattice("cuBIc", "pRIMitIVe")
        assert lattice.lattice_system == LatticeSystem.CUBIC
        assert lattice.centering == Centering.PRIMITIVE

        # default centering argument
        lattice = Lattice(lattice_system="monoclinic")
        assert lattice.lattice_system == LatticeSystem.MONOCLINIC
        assert lattice.centering == Centering.PRIMITIVE

        # centering set as keyword argument
        lattice = Lattice(lattice_system="monoclinic", centering=Centering.BODY)
        assert lattice.lattice_system == LatticeSystem.MONOCLINIC
        assert lattice.centering == Centering.BODY

    @staticmethod
    def test_init_invalid_arguments():
        """
        Test argument checks for `__init__()`
        """
        # --- Tests

        # ------ lattice_system

        # invalid `lattice_system` type
        with pytest.raises(ValueError) as exception_info:
            Lattice(1)

        expected_error = (
            "`lattice_system` must be a LatticeSystem object or a string. "
            "(lattice_system=1)"
        )
        assert expected_error in str(exception_info)

        # invalid `lattice_system` value
        with pytest.raises(ValueError) as exception_info:
            Lattice("abc")

        expected_error = (
            f"`lattice_system` must be one of {LatticeSystem.values()}. "
            "(lattice_system='abc')"
        )
        assert expected_error in str(exception_info)

        # ------ centering

        # invalid `centering` type
        with pytest.raises(ValueError) as exception_info:
            Lattice(LatticeSystem.TETRAGONAL, centering=1)

        expected_error = (
            "`centering` must be a Centering object or a string. (centering=1)"
        )
        assert expected_error in str(exception_info)

        # invalid `centering` value
        # invalid `lattice_system` value
        with pytest.raises(ValueError) as exception_info:
            Lattice("rhombohedral", "def")

        expected_error = (
            "`centering` must be one of {Centering.values()}. (centering='def')"
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
        assert lattice.is_bravais_lattice()

        # (triclinic, body)
        lattice = Lattice(LatticeSystem.TRICLINIC, Centering.BODY)
        assert not lattice.is_bravais_lattice()

        # (triclinic, base)
        lattice = Lattice(LatticeSystem.TRICLINIC, Centering.BASE)
        assert not lattice.is_bravais_lattice()

        # (triclinic, face)
        lattice = Lattice(LatticeSystem.TRICLINIC, Centering.FACE)
        assert not lattice.is_bravais_lattice()

        # --------- monoclinic

        # (monoclinic, primitive)
        lattice = Lattice("monoclinic", "primitive")
        assert lattice.is_bravais_lattice()

        # (monoclinic, body)
        lattice = Lattice("monoclinic", "body")
        assert lattice.is_bravais_lattice()

        # (monoclinic, base)
        lattice = Lattice("monoclinic", "base")
        assert lattice.is_bravais_lattice()

        # (monoclinic, face)
        lattice = Lattice(LatticeSystem.MONOCLINIC, Centering.FACE)
        assert not lattice.is_bravais_lattice()

        # --------- orthorhombic

        # (orthorhombic, primitive)
        lattice = Lattice("orthorhombic", Centering.PRIMITIVE)
        assert lattice.is_bravais_lattice()

        # (orthorhombic, body)
        lattice = Lattice(LatticeSystem.ORTHORHOMBIC, "body")
        assert lattice.is_bravais_lattice()

        # (orthorhombic, base)
        lattice = Lattice(LatticeSystem.ORTHORHOMBIC, Centering.BASE)
        assert lattice.is_bravais_lattice()

        # (orthorhombic, face)
        lattice = Lattice("orthorhombic", "face")
        assert lattice.is_bravais_lattice()

        # --------- tetragonal

        # (tetragonal, primitive)
        lattice = Lattice("tetragonal", Centering.PRIMITIVE)
        assert lattice.is_bravais_lattice()

        # (tetragonal, body)
        lattice = Lattice(LatticeSystem.TETRAGONAL, Centering.BODY)
        assert lattice.is_bravais_lattice()

        # (tetragonal, base)
        lattice = Lattice(LatticeSystem.TETRAGONAL, Centering.BASE)
        assert not lattice.is_bravais_lattice()

        # (tetragonal, face)
        lattice = Lattice(LatticeSystem.TETRAGONAL, Centering.FACE)
        assert not lattice.is_bravais_lattice()

        # --------- rhomobohedral

        # (rhombohedral, primitive)
        lattice = Lattice("rhombohedral", Centering.PRIMITIVE)
        assert lattice.is_bravais_lattice()

        # (rhombohedral, body)
        lattice = Lattice(LatticeSystem.RHOMBOHEDRAL, Centering.BODY)
        assert not lattice.is_bravais_lattice()

        # (rhombohedral, base)
        lattice = Lattice(LatticeSystem.RHOMBOHEDRAL, Centering.BASE)
        assert not lattice.is_bravais_lattice()

        # (rhombohedral, face)
        lattice = Lattice(LatticeSystem.RHOMBOHEDRAL, Centering.FACE)
        assert not lattice.is_bravais_lattice()

        # --------- hexagonal

        # (hexagonal, primitive)
        lattice = Lattice("hexagonal", "primitive")
        assert lattice.is_bravais_lattice()

        # (hexagonal, body)
        lattice = Lattice(LatticeSystem.HEXAGONAL, Centering.BODY)
        assert not lattice.is_bravais_lattice()

        # (hexagonal, base)
        lattice = Lattice(LatticeSystem.HEXAGONAL, Centering.BASE)
        assert not lattice.is_bravais_lattice()

        # (hexagonal, face)
        lattice = Lattice(LatticeSystem.HEXAGONAL, Centering.FACE)
        assert not lattice.is_bravais_lattice()

        # --------- cubic

        # (cubic, primitive)
        lattice = Lattice("cubic", "primitive")
        assert lattice.is_bravais_lattice()

        # (cubic, body)
        lattice = Lattice(LatticeSystem.CUBIC, "body")
        assert lattice.is_bravais_lattice()

        # (cubic, base)
        lattice = Lattice(LatticeSystem.CUBIC, Centering.BASE)
        assert not lattice.is_bravais_lattice()

        # (cubic, face)
        lattice = Lattice("cubic", "face")
        assert lattice.is_bravais_lattice()
