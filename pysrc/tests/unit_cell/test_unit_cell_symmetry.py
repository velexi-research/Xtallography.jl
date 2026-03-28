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
Unit tests for the `xtallography.unit_cell_symmetry.unit_cell_symmetry` module
"""
# --- Imports

# Standard library
from dataclasses import FrozenInstanceError
import unittest

# External packages
import juliacall
import pytest
from xtallography.symmetry import Centering
from xtallography.symmetry import GlidePlane, ScrewAxis
from xtallography.unit_cell import UnitCellSymmetry

# Local packages/modules


# --- Test Suites


class test_xtallography_unit_cell_symmetry_UnitCellSymmetry(unittest.TestCase):
    """
    Test suite for the `UnitCellSymmetry` class
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
        Tests `__init__()`.
        """
        # --- Tests

        # default keyword arguments
        unit_cell_symmetry = UnitCellSymmetry()

        assert unit_cell_symmetry.centering == Centering.PRIMITIVE
        assert unit_cell_symmetry.symmetry_elements == set()

        # all keyword arguments specified
        centering = Centering.FACE
        symmetry_elements = [GlidePlane("1,0,0", "0,1,0"), ScrewAxis("1,0,0", 3, 2)]
        unit_cell_symmetry = UnitCellSymmetry(
            centering=centering, symmetry_elements=symmetry_elements
        )

        assert unit_cell_symmetry.centering == centering
        assert unit_cell_symmetry.symmetry_elements == set(symmetry_elements)

        # keyword arguments: only `centering` specified
        centering = Centering.FACE
        unit_cell_symmetry = UnitCellSymmetry(centering=centering)

        assert unit_cell_symmetry.centering == centering
        assert unit_cell_symmetry.symmetry_elements == set()

        # keyword arguments: only `symmetry_elements` specified as a list
        symmetry_elements = [GlidePlane("1,0,0", "0,1,0"), ScrewAxis("1,0,0", 3, 2)]
        unit_cell_symmetry = UnitCellSymmetry(symmetry_elements=symmetry_elements)

        assert unit_cell_symmetry.centering == Centering.PRIMITIVE
        assert unit_cell_symmetry.symmetry_elements == set(symmetry_elements)

        # keyword arguments: only `symmetry_elements` specified as a set
        symmetry_elements = set(
            [GlidePlane("1,0,0", "0,1,0"), ScrewAxis("1,0,0", 3, 2)]
        )
        unit_cell_symmetry = UnitCellSymmetry(symmetry_elements=symmetry_elements)

        assert unit_cell_symmetry.centering == Centering.PRIMITIVE
        assert unit_cell_symmetry.symmetry_elements == symmetry_elements

    @staticmethod
    def test_init_invalid_arguments():
        """
        Test argument checks for `__init__()`
        """
        # --- Tests

        # `symmetry_elements` contains an element that is not a SymmetryElement object
        symmetry_elements = [
            GlidePlane("1,0,0", "0,1,0"),
            ScrewAxis("1,0,0", 3, 2),
            "not a SymmetryElement",
        ]
        with pytest.raises(ValueError) as exception_info:
            UnitCellSymmetry(Centering.FACE, symmetry_elements)

        expected_error = (
            "`symmetry_elements` contains elements that are not SymmetryElements "
            f"objects. (symmetry_elements={symmetry_elements})"
        )
        assert expected_error in str(exception_info)

    @staticmethod
    def test_frozen():
        """
        Test that UnitCellSymmetry objects are immutable.
        """
        # --- Preparations

        unit_cell_symmetry = UnitCellSymmetry()

        # --- Tests

        # attempt to change `centering` field
        with pytest.raises(FrozenInstanceError) as exception_info:
            unit_cell_symmetry.centering = Centering.FACE

        expected_error = "cannot assign to field 'centering'"
        assert expected_error in str(exception_info)

        # attempt to change `symmetry_elements` field
        with pytest.raises(FrozenInstanceError) as exception_info:
            unit_cell_symmetry.symmetry_elements = set(
                [GlidePlane("1,0,0", "0,1,0"), ScrewAxis("1,0,0", 3, 2)]
            )

        expected_error = "cannot assign to field 'symmetry_elements'"
        assert expected_error in str(exception_info)

    def test_to_julia(self):
        """
        Test `to_julia()`.
        """
        # --- Tests

        # unit cell symmetry: primitive centering, no symmetry elements
        unit_cell_symmetry_jl = UnitCellSymmetry().to_julia()

        assert self.jl.isa(unit_cell_symmetry_jl, self.jl.UnitCellSymmetry)
        assert unit_cell_symmetry_jl.centering == Centering.PRIMITIVE.to_julia()
        assert unit_cell_symmetry_jl.symmetry_elements == self.jl.Set()

        # unit cell symmetry: non-primitive centering, no symmetry elements
        centering = Centering.FACE
        unit_cell_symmetry_jl = UnitCellSymmetry(centering).to_julia()

        assert self.jl.isa(unit_cell_symmetry_jl, self.jl.UnitCellSymmetry)
        assert unit_cell_symmetry_jl.centering == Centering.FACE.to_julia()
        assert unit_cell_symmetry_jl.symmetry_elements == self.jl.Set()

        # unit cell symmetry: non-primitive centering, non-trivial symmetry elements
        centering = Centering.FACE
        symmetry_elements = {GlidePlane("1,0,0", "0,1,0"), ScrewAxis("1,0,0", 3, 2)}
        unit_cell_symmetry_jl = UnitCellSymmetry(
            centering, symmetry_elements
        ).to_julia()

        assert self.jl.isa(unit_cell_symmetry_jl, self.jl.UnitCellSymmetry)
        assert unit_cell_symmetry_jl.centering == Centering.FACE.to_julia()
        assert unit_cell_symmetry_jl.symmetry_elements == self.jl.Set(
            [element.to_julia() for element in symmetry_elements]
        )

    @staticmethod
    def test_from_julia():
        """
        Test `from_julia()`.
        """
        # --- Tests

        # unit cell symmetry: primitive centering, no symmetry elements
        unit_cell_symmetry_jl = UnitCellSymmetry().to_julia()
        unit_cell_symmetry = UnitCellSymmetry.from_julia(unit_cell_symmetry_jl)

        assert isinstance(unit_cell_symmetry, UnitCellSymmetry)
        assert unit_cell_symmetry.centering == Centering.PRIMITIVE
        assert unit_cell_symmetry.symmetry_elements == set()

        # unit cell symmetry: non-primitive centering, no symmetry elements
        centering = Centering.FACE
        unit_cell_symmetry_jl = UnitCellSymmetry(centering).to_julia()
        unit_cell_symmetry = UnitCellSymmetry.from_julia(unit_cell_symmetry_jl)

        assert isinstance(unit_cell_symmetry, UnitCellSymmetry)
        assert unit_cell_symmetry.centering == centering
        assert unit_cell_symmetry.symmetry_elements == set()

        # unit cell symmetry: non-primitive centering, non-trivial symmetry elements
        centering = Centering.FACE
        symmetry_elements = {GlidePlane("1,0,0", "0,1,0"), ScrewAxis("1,0,0", 3, 2)}
        unit_cell_symmetry_jl = UnitCellSymmetry(
            centering, symmetry_elements
        ).to_julia()
        unit_cell_symmetry = UnitCellSymmetry.from_julia(unit_cell_symmetry_jl)

        assert isinstance(unit_cell_symmetry, UnitCellSymmetry)
        assert unit_cell_symmetry.centering == centering
        assert unit_cell_symmetry.symmetry_elements == symmetry_elements

    def test_from_julia_invalid_args(self):
        """
        Test argument checks for `UnitCellSymmetry.from_julia()`.
        """
        # --- Tests

        # ------ `unit_cell_symmetry_jl` is not a Julia UnitCellSymmetry object

        unit_cell_symmetry_jl_invalid = "not a Julia UnitCellSymmetry object"
        with pytest.raises(ValueError) as exception_info:
            UnitCellSymmetry.from_julia(unit_cell_symmetry_jl_invalid)

        expected_error = (
            "`unit_cell_symmetry_jl` must be a Julia `UnitCellSymmetry` object. "
            f"(unit_cell_symmetry_jl={unit_cell_symmetry_jl_invalid})."
        )
        assert expected_error in str(exception_info)

    @staticmethod
    def test__eq__():
        """
        Test `__eq__()`.
        """
        # --- Tests

        # equal UnitCellSymmetry objects
        centering = Centering.FACE
        symmetry_elements = [GlidePlane("1,0,0", "0,1,0"), ScrewAxis("1,0,0", 3, 2)]
        unit_cell_symmetry = UnitCellSymmetry(centering, symmetry_elements)

        other = UnitCellSymmetry(centering, symmetry_elements)

        assert unit_cell_symmetry == other

        # equal UnitCellSymmetry objects
        centering = Centering.FACE
        symmetry_elements = [GlidePlane("1,0,0", "0,1,0"), ScrewAxis("1,0,0", 3, 2)]
        unit_cell_symmetry = UnitCellSymmetry(centering, symmetry_elements)
        other = UnitCellSymmetry(centering)

        assert unit_cell_symmetry != other
