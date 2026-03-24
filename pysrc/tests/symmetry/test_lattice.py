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
from xtallography.symmetry import create_lattice, standardize_lattice
from xtallography.symmetry import is_bravais_lattice

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


class test_xtallography_symmetry_lattice_functions(unittest.TestCase):
    """
    Test suite for lattice functions
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
    def test_create_lattice():
        """
        Test `create_lattice()`.
        """
        # --- Tests

        # Arguments valid and lowercase
        lattice = create_lattice("cubic", "primitive")

        assert lattice.lattice_system == LatticeSystem.CUBIC
        assert lattice.centering == Centering.PRIMITIVE

        # Arguments valid but not lowercase
        lattice = create_lattice("cuBIc", "pRIMitIVe")

        assert lattice.lattice_system == LatticeSystem.CUBIC
        assert lattice.centering == Centering.PRIMITIVE

        # Arguments provided as keyword arguments
        lattice = create_lattice(
            lattice_system="monoclinic",
            centering="body",
        )

        assert lattice.lattice_system == LatticeSystem.MONOCLINIC
        assert lattice.centering == Centering.BODY

        # Arguments provided as LatticeSystem and Centering types
        lattice = create_lattice(LatticeSystem.CUBIC, Centering.PRIMITIVE)

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
            create_lattice()

        expected_error = "`lattice_system` must be a string. (lattice_system=None)"
        assert expected_error in str(exception_info)

        # valid `system`, default `centering`
        with pytest.raises(ValueError) as exception_info:
            create_lattice("cubic")

        expected_error = "`centering` must be a string. (centering=None)"
        assert expected_error in str(exception_info)

        # default `system`, valid `centering` (via keyword argument)
        with pytest.raises(ValueError) as exception_info:
            create_lattice(centering="primitive")

        expected_error = "`lattice_system` must be a string. (lattice_system=None)"
        assert expected_error in str(exception_info)

        # `system` has wrong type
        with pytest.raises(ValueError) as exception_info:
            create_lattice(1, "def")

        expected_error = "`lattice_system` must be a string. (lattice_system=1)"
        assert expected_error in str(exception_info)

        # `centering` has wrong type
        with pytest.raises(ValueError) as exception_info:
            create_lattice("abc", None)

        expected_error = "`centering` must be a string. (centering=None)"
        assert expected_error in str(exception_info)

        # Nonsense `system` and `centering` values
        with pytest.raises(ValueError) as exception_info:
            create_lattice("abc", "def")

        expected_error = (
            f"`lattice_system` must be one of "
            f"{LatticeSystem.values()}. "
            "(lattice_system='abc')"
        )
        assert expected_error in str(exception_info)

        # Nonsense `centering` value
        with pytest.raises(ValueError) as exception_info:
            create_lattice("hexagonal", "def")

        expected_error = (
            "`centering` must be one of " f"{Centering.values()}. " "(centering='def')"
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

        lattice = standardize_lattice(lattice)

        assert isinstance(lattice.lattice_system, LatticeSystem)
        assert lattice.lattice_system == LatticeSystem.CUBIC

        assert isinstance(lattice.centering, Centering)
        assert lattice.centering == Centering.PRIMITIVE

        # Field values are valid but not standardized
        lattice = Lattice("cuBIc", "pRIMitIVe")

        lattice = standardize_lattice(lattice)

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
            standardize_lattice(Lattice(1, "def"))

        expected_error = (
            "`lattice.lattice_system` must be a string. (lattice.lattice_system=1)"
        )
        assert expected_error in str(exception_info)

        # `centering` has wrong type
        with pytest.raises(ValueError) as exception_info:
            standardize_lattice(Lattice("abc", None))

        expected_error = (
            "`lattice.centering` must be a string. (lattice.centering=None)"
        )
        assert expected_error in str(exception_info)

        # Nonsense `system` and `centering` values
        with pytest.raises(ValueError) as exception_info:
            standardize_lattice(Lattice("abc", "def"))

        expected_error = (
            "`lattice_system` must be one of "
            f"{LatticeSystem.values()}. "
            "(lattice_system='abc')"
        )
        assert expected_error in str(exception_info)

        # Nonsense `centering` value
        with pytest.raises(ValueError) as exception_info:
            standardize_lattice(Lattice("hexagonal", "def"))

        expected_error = (
            "`centering` must be one of " f"{Centering.values()}. " "(centering='def')"
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

        assert is_bravais_lattice(lattice)

        # --------- monoclinic

        # (monoclinic, primitive)
        lattice = Lattice("monoclinic", "primitive")

        assert is_bravais_lattice(lattice)

        # (monoclinic, body)
        lattice = Lattice("monoclinic", "body")

        assert is_bravais_lattice(lattice)

        # (monoclinic, base)
        lattice = Lattice("monoclinic", "base")

        assert is_bravais_lattice(lattice)

        # --------- orthorhombic

        # (orthorhombic, primitive)
        lattice = Lattice("orthorhombic", Centering.PRIMITIVE)

        assert is_bravais_lattice(lattice)

        # (orthorhombic, body)
        lattice = Lattice(LatticeSystem.ORTHORHOMBIC, "body")

        assert is_bravais_lattice(lattice)

        # (orthorhombic, base)
        lattice = Lattice(LatticeSystem.ORTHORHOMBIC, Centering.BASE)

        assert is_bravais_lattice(lattice)

        # (orthorhombic, face)
        lattice = Lattice("orthorhombic", "face")

        assert is_bravais_lattice(lattice)

        # --------- tetragonal

        # (tetragonal, primitive)
        lattice = Lattice("tetragonal", Centering.PRIMITIVE)

        assert is_bravais_lattice(lattice)

        # (tetragonal, body)
        lattice = Lattice(LatticeSystem.TETRAGONAL, Centering.BODY)

        assert is_bravais_lattice(lattice)

        # --------- rhomobohedral

        # (rhombohedral, primitive)
        lattice = Lattice("rhombohedral", Centering.PRIMITIVE)

        assert is_bravais_lattice(lattice)

        # --------- hexagonal

        # (hexagonal, primitive)
        lattice = Lattice("hexagonal", "primitive")

        assert is_bravais_lattice(lattice)

        # --------- cubic

        # (cubic, primitive)
        lattice = Lattice("cubic", "primitive")

        assert is_bravais_lattice(lattice)

        # (cubic, body)
        lattice = Lattice(LatticeSystem.CUBIC, "body")

        assert is_bravais_lattice(lattice)

        # (cubic, face)
        lattice = Lattice("cubic", "face")

        assert is_bravais_lattice(lattice)

    @staticmethod
    def test_is_bravais_lattice_invalid_lattice():
        """
        Test `is_bravais_lattice()`.
        """
        # --- Tests

        # --------- triclinic

        # (triclinic, body)
        lattice = Lattice(LatticeSystem.TRICLINIC, Centering.BODY)
        assert not is_bravais_lattice(lattice)

        # (triclinic, base)
        lattice = Lattice(LatticeSystem.TRICLINIC, Centering.BASE)
        assert not is_bravais_lattice(lattice)

        # (triclinic, face)
        lattice = Lattice(LatticeSystem.TRICLINIC, Centering.FACE)
        assert not is_bravais_lattice(lattice)

        # --------- monoclinic

        # (monoclinic, face)
        lattice = Lattice(LatticeSystem.MONOCLINIC, Centering.FACE)
        assert not is_bravais_lattice(lattice)

        # --------- orthorhombic

        # all centerings are valid

        # --------- tetragonal

        # (tetragonal, base)
        lattice = Lattice(LatticeSystem.TETRAGONAL, Centering.BASE)
        assert not is_bravais_lattice(lattice)

        # (tetragonal, face)
        lattice = Lattice(LatticeSystem.TETRAGONAL, Centering.FACE)
        assert not is_bravais_lattice(lattice)

        # --------- rhomobohedral

        # (rhombohedral, body)
        lattice = Lattice(LatticeSystem.RHOMBOHEDRAL, Centering.BODY)
        assert not is_bravais_lattice(lattice)

        # (rhombohedral, base)
        lattice = Lattice(LatticeSystem.RHOMBOHEDRAL, Centering.BASE)
        assert not is_bravais_lattice(lattice)

        # (rhombohedral, face)
        lattice = Lattice(LatticeSystem.RHOMBOHEDRAL, Centering.FACE)
        assert not is_bravais_lattice(lattice)

        # --------- hexagonal

        # (hexagonal, body)
        lattice = Lattice(LatticeSystem.HEXAGONAL, Centering.BODY)
        assert not is_bravais_lattice(lattice)

        # (hexagonal, base)
        lattice = Lattice(LatticeSystem.HEXAGONAL, Centering.BASE)
        assert not is_bravais_lattice(lattice)

        # (hexagonal, face)
        lattice = Lattice(LatticeSystem.HEXAGONAL, Centering.FACE)
        assert not is_bravais_lattice(lattice)

        # --------- cubic

        # (cubic, base)
        lattice = Lattice(LatticeSystem.CUBIC, Centering.BASE)
        assert not is_bravais_lattice(lattice)
