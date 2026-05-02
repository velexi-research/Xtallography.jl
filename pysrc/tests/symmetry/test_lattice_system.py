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
Unit tests for `xtallography.symmetry.lattice_system` module
"""
# --- Imports

# Standard library
import unittest

# External packages
import juliacall
import pytest
from xtallography.symmetry import LatticeSystem

# Local packages/modules


# --- Test Suites


class test_xtallography_symmetry_lattice_system_LatticeSystem(unittest.TestCase):
    """
    Test suite for the `LatticeSystem` class
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
            LatticeSystem.from_julia(lattice_system_jl_invalid)

        expected_error = (
            "`lattice_system_jl` must be a Julia `LatticeSystem` object. "
            f"(lattice_system_jl={lattice_system_jl_invalid})."
        )
        assert expected_error in str(exception_info)

        # ------ `lattice_system_jl` is an unsupported Julia LatticeSystem object

        self.jl.seval("struct UnsupportedLatticeSystem <: LatticeSystem end")
        lattice_system_jl_invalid = self.jl.UnsupportedLatticeSystem()
        with pytest.raises(ValueError) as exception_info:
            LatticeSystem.from_julia(lattice_system_jl_invalid)

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
