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
Unit tests for `xtallography.unit_cell.rhombohedral` module
"""
# --- Imports

# Standard library
import math
import unittest

# External packages
import juliacall
import pytest

# Local packages/modules
from xtallography.symmetry import LatticeSystem, Centering
from xtallography.symmetry import GlidePlane, ScrewAxis, SymmetryElement
from xtallography.unit_cell import (
    RhombohedralLatticeConstants,
    RhombohedralUnitCell,
    UnitCellSymmetry,
    TetragonalUnitCell,
)


# --- Test Suites


class test_xtallography_unit_cell_rhombohedral(unittest.TestCase):
    """
    Test suite for the `RhombohedralUnitCell` class
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

        Also tests:
        * `centering` and `symmetry_elements` properties
        * `a`, `alpha` properties
        """
        # --- Preparations

        a = 1
        alpha = 0.1

        # --- Tests

        # ------ default keyword arguments

        unit_cell = RhombohedralUnitCell(a, alpha)

        assert unit_cell.lattice_constants == RhombohedralLatticeConstants(a, alpha)
        assert unit_cell.lattice_system == LatticeSystem.RHOMBOHEDRAL
        assert unit_cell.symmetry == UnitCellSymmetry(
            centering=Centering.PRIMITIVE, symmetry_elements=set()
        )

        assert unit_cell.a == a
        assert unit_cell.alpha == alpha
        assert unit_cell.centering == Centering.PRIMITIVE
        assert unit_cell.symmetry_elements == set()

        # ------ non-default keyword arguments

        # centering = primitive
        unit_cell = RhombohedralUnitCell(a, alpha, centering=Centering.PRIMITIVE)

        assert unit_cell.lattice_constants == RhombohedralLatticeConstants(a, alpha)
        assert unit_cell.lattice_system == LatticeSystem.RHOMBOHEDRAL
        assert unit_cell.symmetry == UnitCellSymmetry(
            centering=Centering.PRIMITIVE, symmetry_elements=set()
        )

        assert unit_cell.a == a
        assert unit_cell.alpha == alpha
        assert unit_cell.centering == Centering.PRIMITIVE
        assert unit_cell.symmetry_elements == set()

        # centering = base
        unit_cell = RhombohedralUnitCell(a, alpha, centering=Centering.BASE)

        assert unit_cell.lattice_constants == RhombohedralLatticeConstants(a, alpha)
        assert unit_cell.lattice_system == LatticeSystem.RHOMBOHEDRAL
        assert unit_cell.symmetry == UnitCellSymmetry(
            centering=Centering.BASE, symmetry_elements=set()
        )

        assert unit_cell.a == a
        assert unit_cell.alpha == alpha
        assert unit_cell.centering == Centering.BASE
        assert unit_cell.symmetry_elements == set()

        # centering = body
        unit_cell = RhombohedralUnitCell(a, alpha, centering=Centering.BODY)

        assert unit_cell.lattice_constants == RhombohedralLatticeConstants(a, alpha)
        assert unit_cell.lattice_system == LatticeSystem.RHOMBOHEDRAL
        assert unit_cell.symmetry == UnitCellSymmetry(
            centering=Centering.BODY, symmetry_elements=set()
        )

        assert unit_cell.a == a
        assert unit_cell.alpha == alpha
        assert unit_cell.centering == Centering.BODY
        assert unit_cell.symmetry_elements == set()

        # centering = face
        unit_cell = RhombohedralUnitCell(a, alpha, centering=Centering.FACE)

        assert unit_cell.lattice_constants == RhombohedralLatticeConstants(a, alpha)
        assert unit_cell.lattice_system == LatticeSystem.RHOMBOHEDRAL
        assert unit_cell.symmetry == UnitCellSymmetry(
            centering=Centering.FACE, symmetry_elements=set()
        )

        assert unit_cell.a == a
        assert unit_cell.alpha == alpha
        assert unit_cell.centering == Centering.FACE
        assert unit_cell.symmetry_elements == set()

        # symmetry_elements non-empty
        symmetry_elements = [GlidePlane("1,0,0", "0,1,0"), ScrewAxis("1,0,0", 3, 2)]
        unit_cell = RhombohedralUnitCell(a, alpha, symmetry_elements=symmetry_elements)

        assert unit_cell.lattice_constants == RhombohedralLatticeConstants(a, alpha)
        assert unit_cell.lattice_system == LatticeSystem.RHOMBOHEDRAL
        assert unit_cell.symmetry == UnitCellSymmetry(
            centering=Centering.PRIMITIVE, symmetry_elements=symmetry_elements
        )

        assert unit_cell.a == a
        assert unit_cell.alpha == alpha
        assert unit_cell.centering == Centering.PRIMITIVE
        assert unit_cell.symmetry_elements == set(symmetry_elements)

    @staticmethod
    def test_init_invalid_args():
        """
        Test `__init__()`: invalid arguments.
        """
        # --- Preparations

        a = 1
        alpha = 0.1

        # --- Tests

        # ------ Invalid `a`

        # a < 0
        invalid_value = -1
        with pytest.raises(ValueError) as exception_info:
            RhombohedralUnitCell(invalid_value, alpha)

        expected_error = f"`a` must be positive. (a={invalid_value})"
        assert expected_error in str(exception_info)

        # a = 0
        invalid_value = 0
        with pytest.raises(ValueError) as exception_info:
            RhombohedralUnitCell(invalid_value, alpha)

        expected_error = f"`a` must be positive. (a={invalid_value})"
        assert expected_error in str(exception_info)

        # ------ Invalid `alpha`

        # alpha < 0
        invalid_value = -1
        with pytest.raises(ValueError) as exception_info:
            RhombohedralUnitCell(a, invalid_value)

        expected_error = (
            f"`alpha` must lie in the interval (0, 2 pi). (alpha={invalid_value})"
        )
        assert expected_error in str(exception_info)

        # alpha = 0
        invalid_value = 0
        with pytest.raises(ValueError) as exception_info:
            RhombohedralUnitCell(a, invalid_value)

        expected_error = (
            f"`alpha` must lie in the interval (0, 2 pi). (alpha={invalid_value})"
        )
        assert expected_error in str(exception_info)

        # alpha = 2 pi
        invalid_value = 2 * math.pi
        with pytest.raises(ValueError) as exception_info:
            RhombohedralUnitCell(a, invalid_value)

        expected_error = (
            f"`alpha` must lie in the interval (0, 2 pi). (alpha={invalid_value})"
        )
        assert expected_error in str(exception_info)

        # alpha > 2 pi
        invalid_value = 7
        with pytest.raises(ValueError) as exception_info:
            RhombohedralUnitCell(a, invalid_value)

        expected_error = (
            f"`alpha` must lie in the interval (0, 2 pi). (alpha={invalid_value})"
        )
        assert expected_error in str(exception_info)

    def test_to_julia(self):
        """
        Test `to_julia()`.
        """
        # --- Preparations

        a = 1
        alpha = 0.1

        # ------ default centering and symmetry elements

        unit_cell = RhombohedralUnitCell(a, alpha)
        unit_cell_jl = unit_cell.to_julia()

        assert self.jl.isa(unit_cell_jl, self.jl.UnitCell)
        assert self.jl.isa(unit_cell_jl, self.jl.RhombohedralUnitCell)
        assert unit_cell_jl.lattice_constants.a == a
        assert unit_cell_jl.lattice_constants.α == alpha
        assert unit_cell_jl.symmetry.centering == self.jl.primitive_centering
        assert (
            set(
                [
                    SymmetryElement.from_julia(element)
                    for element in unit_cell_jl.symmetry.symmetry_elements
                ]
            )
            == set()
        )

        # ------ non-default centering

        # centering = primitive
        unit_cell = RhombohedralUnitCell(a, alpha, centering=Centering.PRIMITIVE)
        unit_cell_jl = unit_cell.to_julia()

        assert self.jl.isa(unit_cell_jl, self.jl.UnitCell)
        assert self.jl.isa(unit_cell_jl, self.jl.RhombohedralUnitCell)
        assert unit_cell_jl.lattice_constants.a == a
        assert unit_cell_jl.lattice_constants.α == alpha
        assert unit_cell_jl.symmetry.centering == self.jl.primitive_centering
        assert (
            set(
                [
                    SymmetryElement.from_julia(element)
                    for element in unit_cell_jl.symmetry.symmetry_elements
                ]
            )
            == set()
        )

        # centering = base
        unit_cell = RhombohedralUnitCell(a, alpha, centering=Centering.BASE)
        unit_cell_jl = unit_cell.to_julia()

        assert self.jl.isa(unit_cell_jl, self.jl.UnitCell)
        assert self.jl.isa(unit_cell_jl, self.jl.RhombohedralUnitCell)
        assert unit_cell_jl.lattice_constants.a == a
        assert unit_cell_jl.lattice_constants.α == alpha
        assert unit_cell_jl.symmetry.centering == self.jl.base_centering
        assert (
            set(
                [
                    SymmetryElement.from_julia(element)
                    for element in unit_cell_jl.symmetry.symmetry_elements
                ]
            )
            == set()
        )

        # centering = body
        unit_cell = RhombohedralUnitCell(a, alpha, centering=Centering.BODY)
        unit_cell_jl = unit_cell.to_julia()

        assert self.jl.isa(unit_cell_jl, self.jl.UnitCell)
        assert self.jl.isa(unit_cell_jl, self.jl.RhombohedralUnitCell)
        assert unit_cell_jl.lattice_constants.a == a
        assert unit_cell_jl.lattice_constants.α == alpha
        assert unit_cell_jl.symmetry.centering == self.jl.body_centering
        assert (
            set(
                [
                    SymmetryElement.from_julia(element)
                    for element in unit_cell_jl.symmetry.symmetry_elements
                ]
            )
            == set()
        )

        # centering = face
        unit_cell = RhombohedralUnitCell(a, alpha, centering=Centering.FACE)
        unit_cell_jl = unit_cell.to_julia()

        assert self.jl.isa(unit_cell_jl, self.jl.UnitCell)
        assert self.jl.isa(unit_cell_jl, self.jl.RhombohedralUnitCell)
        assert unit_cell_jl.lattice_constants.a == a
        assert unit_cell_jl.lattice_constants.α == alpha
        assert unit_cell_jl.symmetry.centering == self.jl.face_centering
        assert (
            set(
                [
                    SymmetryElement.from_julia(element)
                    for element in unit_cell_jl.symmetry.symmetry_elements
                ]
            )
            == set()
        )

        # ------ non-default symmetry elements

        symmetry_elements = [GlidePlane("1,0,0", "0,1,0"), ScrewAxis("1,0,0", 3, 2)]
        unit_cell = RhombohedralUnitCell(a, alpha, symmetry_elements=symmetry_elements)
        unit_cell_jl = unit_cell.to_julia()

        assert self.jl.isa(unit_cell_jl, self.jl.UnitCell)
        assert self.jl.isa(unit_cell_jl, self.jl.RhombohedralUnitCell)
        assert unit_cell_jl.lattice_constants.a == a
        assert unit_cell_jl.lattice_constants.α == alpha
        assert unit_cell_jl.symmetry.centering == self.jl.primitive_centering
        assert set(
            [
                SymmetryElement.from_julia(element)
                for element in unit_cell_jl.symmetry.symmetry_elements
            ]
        ) == set(symmetry_elements)

    def test_from_julia(self):
        """
        Test `from_julia()`.
        """
        # --- Preparations

        # lattice constants
        a = 1
        alpha = 0.1

        # --- Tests

        # ------ default centering and symmetry elements

        unit_cell_jl = self.jl.RhombohedralUnitCell(a, alpha)

        unit_cell = RhombohedralUnitCell.from_julia(unit_cell_jl)
        assert unit_cell == RhombohedralUnitCell(a, alpha)

        # ------ non-default centering

        # centering = primitive
        unit_cell_jl = self.jl.RhombohedralUnitCell(
            a, alpha, centering=self.jl.primitive_centering
        )
        unit_cell = RhombohedralUnitCell.from_julia(unit_cell_jl)
        assert unit_cell == RhombohedralUnitCell(
            a, alpha, centering=Centering.PRIMITIVE
        )

        # centering = base
        unit_cell_jl = self.jl.RhombohedralUnitCell(
            a, alpha, centering=self.jl.base_centering
        )
        unit_cell = RhombohedralUnitCell.from_julia(unit_cell_jl)
        assert unit_cell == RhombohedralUnitCell(a, alpha, centering=Centering.BASE)

        # centering = body
        unit_cell_jl = self.jl.RhombohedralUnitCell(
            a, alpha, centering=self.jl.body_centering
        )
        unit_cell = RhombohedralUnitCell.from_julia(unit_cell_jl)
        assert unit_cell == RhombohedralUnitCell(a, alpha, centering=Centering.BODY)

        # centering = face
        unit_cell_jl = self.jl.RhombohedralUnitCell(
            a, alpha, centering=self.jl.face_centering
        )
        unit_cell = RhombohedralUnitCell.from_julia(unit_cell_jl)
        assert unit_cell == RhombohedralUnitCell(a, alpha, centering=Centering.FACE)

        # ------ non-default symmetry elements

        symmetry_elements = [GlidePlane("1,0,0", "0,1,0"), ScrewAxis("1,0,0", 3, 2)]
        symmetry_elements_jl = self.jl.Vector(
            [element.to_julia() for element in symmetry_elements]
        )
        unit_cell_jl = self.jl.RhombohedralUnitCell(
            a,
            alpha,
            symmetry_elements=symmetry_elements_jl,
        )

        unit_cell = RhombohedralUnitCell.from_julia(unit_cell_jl)
        assert unit_cell == RhombohedralUnitCell(
            a, alpha, symmetry_elements=symmetry_elements
        )

    def test_from_julia_invalid_arguments(self):
        """
        Test `from_julia()`.
        """
        # --- Tests

        # unit_cell_jl is not a RhombohedralUnitCell
        unit_cell_jl_invalid = self.jl.CubicUnitCell(1)

        with pytest.raises(ValueError) as exception_info:
            RhombohedralUnitCell.from_julia(unit_cell_jl_invalid)

        expected_error = (
            "`unit_cell_jl` must be a Julia `RhombohedralUnitCell` object. "
            f"(unit_cell_jl={unit_cell_jl_invalid})."
        )
        assert expected_error in str(exception_info)

    @staticmethod
    def test_repr():
        """
        Test `__repr__()`.
        """
        # --- Preparations

        # lattice constants
        a = 1
        alpha = 0.1

        # --- Tests

        # ------ default centering and symmetry

        unit_cell = RhombohedralUnitCell(a, alpha)
        assert str(unit_cell) == (
            f"RhombohedralUnitCell(a={a},alpha={alpha},"
            "centering=primitive,"
            "symmetry_elements=[])"
        )

        # ------ non-default centering

        # centering = primitive
        unit_cell = RhombohedralUnitCell(a, alpha, centering=Centering.PRIMITIVE)
        assert str(unit_cell) == (
            f"RhombohedralUnitCell(a={a},alpha={alpha},"
            "centering=primitive,"
            "symmetry_elements=[])"
        )

        # centering = base
        unit_cell = RhombohedralUnitCell(a, alpha, centering=Centering.BASE)
        assert str(unit_cell) == (
            f"RhombohedralUnitCell(a={a},alpha={alpha},"
            "centering=base,"
            "symmetry_elements=[])"
        )

        # centering = body
        unit_cell = RhombohedralUnitCell(a, alpha, centering=Centering.BODY)
        assert str(unit_cell) == (
            f"RhombohedralUnitCell(a={a},alpha={alpha},"
            "centering=body,"
            "symmetry_elements=[])"
        )

        # centering = face
        unit_cell = RhombohedralUnitCell(a, alpha, centering=Centering.FACE)
        assert str(unit_cell) == (
            f"RhombohedralUnitCell(a={a},alpha={alpha},"
            "centering=face,"
            "symmetry_elements=[])"
        )

        # ------ non-default symmetry elements

        symmetry_elements = [GlidePlane("1,0,0", "0,1,0"), ScrewAxis("1,0,0", 3, 2)]
        unit_cell = RhombohedralUnitCell(a, alpha, symmetry_elements=symmetry_elements)
        assert str(unit_cell) == (
            f"RhombohedralUnitCell(a={a},alpha={alpha},"
            "centering=primitive,"
            "symmetry_elements=["
            "GlidePlane(translation='1,0,0',reflection_plane='0,1,0'),"
            "ScrewAxis(axis='1,0,0',n=3,m=2)"
            "])"
        )

    @staticmethod
    def test_eq():
        """
        Test `__eq__()`.
        """
        # --- Preparations

        # lattice constants
        a = 1
        alpha = 0.1

        # --- Tests

        # ------ types differ

        unit_cell_1 = RhombohedralUnitCell(a, alpha)
        unit_cell_2 = TetragonalUnitCell(a, a + 2)
        assert unit_cell_1 != unit_cell_2

        # ------ lattice constants and symmmetry are the same

        unit_cell_1 = RhombohedralUnitCell(a, alpha)
        unit_cell_2 = RhombohedralUnitCell(a, alpha)
        assert unit_cell_1 == unit_cell_2

        # ------ lattice constants differ

        # `a` values differ
        unit_cell_1 = RhombohedralUnitCell(a, alpha)
        unit_cell_2 = RhombohedralUnitCell(a + 1, alpha)
        assert unit_cell_1 != unit_cell_2

        # `alpha` values differ
        unit_cell_1 = RhombohedralUnitCell(a, alpha)
        unit_cell_2 = RhombohedralUnitCell(a, alpha + 0.1)
        assert unit_cell_1 != unit_cell_2

        # ------ symmetry differs

        # centerings differ
        unit_cell_1 = RhombohedralUnitCell(a, alpha, centering=Centering.FACE)
        unit_cell_2 = RhombohedralUnitCell(a, alpha, centering=Centering.BODY)
        assert unit_cell_1 != unit_cell_2

        # symmetry elements differ
        symmetry_elements = [GlidePlane("1,0,0", "0,1,0"), ScrewAxis("1,0,0", 3, 2)]
        unit_cell_1 = RhombohedralUnitCell(
            a, alpha, symmetry_elements=symmetry_elements
        )
        unit_cell_2 = RhombohedralUnitCell(
            a, alpha, symmetry_elements=symmetry_elements[0:-1]
        )
        assert unit_cell_1 != unit_cell_2

    @staticmethod
    def test_isclose():
        """
        Test `isclose()`.
        """
        # --- Preparations

        # lattice constants
        a = 1
        alpha = 0.1

        # --- Tests

        # ------ types differ

        unit_cell_1 = RhombohedralUnitCell(a, alpha)
        unit_cell_2 = TetragonalUnitCell(a, a + 1)
        assert not unit_cell_1.isclose(unit_cell_2)

        # ------ `a`

        # `a` values the equal to within tolerance
        unit_cell_1 = RhombohedralUnitCell(a, alpha)
        unit_cell_2 = RhombohedralUnitCell(a + 0.1, alpha)
        assert unit_cell_1.isclose(unit_cell_2, atol=0.2)

        # `a` values differ by more than tolerance
        unit_cell_1 = RhombohedralUnitCell(a, alpha)
        unit_cell_2 = RhombohedralUnitCell(a + 1, alpha)
        assert not unit_cell_1.isclose(unit_cell_2, atol=0.2)

        # ------ `alpha`

        # `alpha` values the equal to within tolerance
        unit_cell_1 = RhombohedralUnitCell(a, alpha)
        unit_cell_2 = RhombohedralUnitCell(a, alpha + 0.1)
        assert unit_cell_1.isclose(unit_cell_2, atol=0.2)

        # `alpha` values differ by more than tolerance
        unit_cell_1 = RhombohedralUnitCell(a, alpha)
        unit_cell_2 = RhombohedralUnitCell(a, alpha + 0.3)
        assert not unit_cell_1.isclose(unit_cell_2, atol=0.2)

        # ------ symmetry differs

        # centerings differ
        unit_cell_1 = RhombohedralUnitCell(a, alpha, centering=Centering.FACE)
        unit_cell_2 = RhombohedralUnitCell(a, alpha, centering=Centering.BODY)
        assert not unit_cell_1.isclose(unit_cell_2)

        # symmetry elements differ
        symmetry_elements = [GlidePlane("1,0,0", "0,1,0"), ScrewAxis("1,0,0", 3, 2)]
        unit_cell_1 = RhombohedralUnitCell(
            a, alpha, symmetry_elements=symmetry_elements
        )
        unit_cell_2 = RhombohedralUnitCell(
            a, alpha, symmetry_elements=symmetry_elements[0:-1]
        )
        assert not unit_cell_1.isclose(unit_cell_2)
