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
Unit tests for `xtallography.lattices.core` module (excluding UnitCell tests)
"""
# --- Imports

# Standard library
import unittest

# External packages
import juliacall
import pytest
import xtallography
from xtallography.lattices import LatticeSystem, Centering, Lattice

# Local packages/modules


# --- Test Suites


class test_xtallography_lattices_core_LatticeSystem(unittest.TestCase):
    """
    Test suite for the `xtallography.lattices.core.LatticeSystem` class
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

    def test_to_julia(self):
        """
        Test `to_julia()` method.
        """
        # --- Tests

        # triclinic
        lattice_system = LatticeSystem.TRICLINIC
        lattice_system_jl = lattice_system.to_julia()
        assert self.jl.isa(lattice_system_jl, self.jl.Triclinic)

        # monoclinic
        lattice_system = LatticeSystem.MONOCLINIC
        lattice_system_jl = lattice_system.to_julia()
        assert self.jl.isa(lattice_system_jl, self.jl.Monoclinic)

        # orthorhombic
        lattice_system = LatticeSystem.ORTHORHOMBIC
        lattice_system_jl = lattice_system.to_julia()
        assert self.jl.isa(lattice_system_jl, self.jl.Orthorhombic)

        # tetragonal
        lattice_system = LatticeSystem.TETRAGONAL
        lattice_system_jl = lattice_system.to_julia()
        assert self.jl.isa(lattice_system_jl, self.jl.Tetragonal)

        # rhombohedral
        lattice_system = LatticeSystem.RHOMBOHEDRAL
        lattice_system_jl = lattice_system.to_julia()
        assert self.jl.isa(lattice_system_jl, self.jl.Rhombohedral)

        # hexagonal
        lattice_system = LatticeSystem.HEXAGONAL
        lattice_system_jl = lattice_system.to_julia()
        assert self.jl.isa(lattice_system_jl, self.jl.Hexagonal)

        # cubic
        lattice_system = LatticeSystem.CUBIC
        lattice_system_jl = lattice_system.to_julia()
        assert self.jl.isa(lattice_system_jl, self.jl.Cubic)

    def test_from_julia(self):
        """
        Test `from_julia()` method.
        """
        # --- Tests

        # triclinic
        lattice_system_jl = self.jl.triclinic
        lattice_system = LatticeSystem.from_julia(lattice_system_jl)
        assert lattice_system == LatticeSystem.TRICLINIC

        # monoclinic
        lattice_system_jl = self.jl.monoclinic
        lattice_system = LatticeSystem.from_julia(lattice_system_jl)
        assert lattice_system == LatticeSystem.MONOCLINIC

        # orthorhombic
        lattice_system_jl = self.jl.orthorhombic
        lattice_system = LatticeSystem.from_julia(lattice_system_jl)
        assert lattice_system == LatticeSystem.ORTHORHOMBIC

        # tetragonal
        lattice_system_jl = self.jl.tetragonal
        lattice_system = LatticeSystem.from_julia(lattice_system_jl)
        assert lattice_system == LatticeSystem.TETRAGONAL

        # rhombohedral
        lattice_system_jl = self.jl.rhombohedral
        lattice_system = LatticeSystem.from_julia(lattice_system_jl)
        assert lattice_system == LatticeSystem.RHOMBOHEDRAL

        # hexagonal
        lattice_system_jl = self.jl.hexagonal
        lattice_system = LatticeSystem.from_julia(lattice_system_jl)
        assert lattice_system == LatticeSystem.HEXAGONAL

        # cubic
        lattice_system_jl = self.jl.cubic
        lattice_system = LatticeSystem.from_julia(lattice_system_jl)
        assert lattice_system == LatticeSystem.CUBIC

    def test_from_julia_invalid_args(self):
        """
        Test argument checks for `from_julia()`.
        """
        # --- Tests

        # ------ `lattice_system_jl` is not a Julia LatticeSystem object

        lattice_system_jl_invalid = "not a Julia LatticeSystem object"
        with pytest.raises(ValueError) as exception_info:
            xtallography.lattices.LatticeSystem.from_julia(lattice_system_jl_invalid)

        expected_error = (
            "`lattice_system_jl` must be a Julia `LatticeSystem` object. "
            f"(lattice_system_jl={lattice_system_jl_invalid})."
        )
        assert expected_error in str(exception_info)

        # ------ `lattice_system_jl` is an unsupported Julia LatticeSystem object

        self.jl.seval("struct UnsupportedLatticeSystem <: LatticeSystem end")
        lattice_system_jl_invalid = self.jl.UnsupportedLatticeSystem()
        with pytest.raises(ValueError) as exception_info:
            xtallography.lattices.LatticeSystem.from_julia(lattice_system_jl_invalid)

        expected_error = (
            "Unsupported LatticeSystem type. "
            "(lattice_system_jl=UnsupportedLatticeSystem)"
        )
        assert expected_error in str(exception_info)

    @staticmethod
    def test_values():
        """
        Test values() method.
        """
        # --- Tests

        # Get full list of values
        values = set(LatticeSystem.values())
        assert values == {
            "triclinic",
            "monoclinic",
            "orthorhombic",
            "tetragonal",
            "rhombohedral",
            "hexagonal",
            "cubic",
        }


class test_xtallography_lattices_core_Centering(unittest.TestCase):
    """
    Test suite for the `xtallography.lattices.core.Centering` class
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

    def test_to_julia(self):
        """
        Test `to_julia()` method.
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

    def test_from_julia(self):
        """
        Test `from_julia()` method.
        """
        # --- Tests

        # primitive
        centering_jl = self.jl.primitive
        centering = Centering.from_julia(centering_jl)
        assert centering == Centering.PRIMITIVE

        # base_centered
        centering_jl = self.jl.base_centered
        centering = Centering.from_julia(centering_jl)
        assert centering == Centering.BASE_CENTERED

        # body_centered
        centering_jl = self.jl.body_centered
        centering = Centering.from_julia(centering_jl)
        assert centering == Centering.BODY_CENTERED

        # face_centered
        centering_jl = self.jl.face_centered
        centering = Centering.from_julia(centering_jl)
        assert centering == Centering.FACE_CENTERED

    def test_from_julia_invalid_args(self):
        """
        Test argument checks for `from_julia()`.
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

        self.jl.seval("struct UnsupportedCentering <: Centering end")
        centering_jl_invalid = self.jl.UnsupportedCentering()
        with pytest.raises(ValueError) as exception_info:
            xtallography.lattices.Centering.from_julia(centering_jl_invalid)

        expected_error = (
            "Unsupported Centering type. (centering_jl=UnsupportedCentering)"
        )
        assert expected_error in str(exception_info)

    @staticmethod
    def test_values():
        """
        Test values() method.
        """
        # --- Tests

        # Get full list of values
        values = set(Centering.values())
        assert values == {
            "primitive",
            "base_centered",
            "body_centered",
            "face_centered",
        }


class test_xtallography_lattices_core_Lattice(unittest.TestCase):
    """
    Test suite for the `xtallography.lattices.core.Lattice` class
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


class test_xtallography_lattices_core_functions(unittest.TestCase):
    """
    Test suite for the `xtallography.lattices.core` module functions (excluding UnitCell
    tests)
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
