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
Unit tests for `xtallography.unit_cell.monoclinic` module
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
    MonoclinicLatticeConstants,
    MonoclinicUnitCell,
    UnitCellSymmetry,
    TetragonalUnitCell,
)


# --- Test Suites


class test_xtallography_unit_cell_monoclinic(unittest.TestCase):
    """
    Test suite for the `MonoclinicUnitCell` class
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
        * `a`, `b`, `c`, `beta` properties
        """
        # --- Preparations

        a = 1
        b = 2
        c = 3
        beta = 0.2

        # --- Tests

        # ------ default keyword arguments

        unit_cell = MonoclinicUnitCell(a, b, c, beta)

        assert unit_cell.lattice_constants == MonoclinicLatticeConstants(a, b, c, beta)
        assert unit_cell.lattice_system == LatticeSystem.MONOCLINIC
        assert unit_cell.symmetry == UnitCellSymmetry(
            centering=Centering.PRIMITIVE, symmetry_elements=set()
        )

        assert unit_cell.a == a
        assert unit_cell.b == b
        assert unit_cell.c == c
        assert unit_cell.beta == beta
        assert unit_cell.centering == Centering.PRIMITIVE
        assert unit_cell.symmetry_elements == set()

        # ------ non-default keyword arguments

        # centering = primitive
        unit_cell = MonoclinicUnitCell(a, b, c, beta, centering=Centering.PRIMITIVE)

        assert unit_cell.lattice_constants == MonoclinicLatticeConstants(a, b, c, beta)
        assert unit_cell.lattice_system == LatticeSystem.MONOCLINIC
        assert unit_cell.symmetry == UnitCellSymmetry(
            centering=Centering.PRIMITIVE, symmetry_elements=set()
        )

        assert unit_cell.a == a
        assert unit_cell.b == b
        assert unit_cell.c == c
        assert unit_cell.beta == beta
        assert unit_cell.centering == Centering.PRIMITIVE
        assert unit_cell.symmetry_elements == set()

        # centering = base
        unit_cell = MonoclinicUnitCell(a, b, c, beta, centering=Centering.BASE)

        assert unit_cell.lattice_constants == MonoclinicLatticeConstants(a, b, c, beta)
        assert unit_cell.lattice_system == LatticeSystem.MONOCLINIC
        assert unit_cell.symmetry == UnitCellSymmetry(
            centering=Centering.BASE, symmetry_elements=set()
        )

        assert unit_cell.a == a
        assert unit_cell.b == b
        assert unit_cell.c == c
        assert unit_cell.beta == beta
        assert unit_cell.centering == Centering.BASE
        assert unit_cell.symmetry_elements == set()

        # centering = body
        unit_cell = MonoclinicUnitCell(a, b, c, beta, centering=Centering.BODY)

        assert unit_cell.lattice_constants == MonoclinicLatticeConstants(a, b, c, beta)
        assert unit_cell.lattice_system == LatticeSystem.MONOCLINIC
        assert unit_cell.symmetry == UnitCellSymmetry(
            centering=Centering.BODY, symmetry_elements=set()
        )

        assert unit_cell.a == a
        assert unit_cell.b == b
        assert unit_cell.c == c
        assert unit_cell.beta == beta
        assert unit_cell.centering == Centering.BODY
        assert unit_cell.symmetry_elements == set()

        # centering = face
        unit_cell = MonoclinicUnitCell(a, b, c, beta, centering=Centering.FACE)

        assert unit_cell.lattice_constants == MonoclinicLatticeConstants(a, b, c, beta)
        assert unit_cell.lattice_system == LatticeSystem.MONOCLINIC
        assert unit_cell.symmetry == UnitCellSymmetry(
            centering=Centering.FACE, symmetry_elements=set()
        )

        assert unit_cell.a == a
        assert unit_cell.b == b
        assert unit_cell.c == c
        assert unit_cell.beta == beta
        assert unit_cell.centering == Centering.FACE
        assert unit_cell.symmetry_elements == set()

        # symmetry_elements non-empty
        symmetry_elements = [GlidePlane("1,0,0", "0,1,0"), ScrewAxis("1,0,0", 3, 2)]
        unit_cell = MonoclinicUnitCell(
            a, b, c, beta, symmetry_elements=symmetry_elements
        )

        assert unit_cell.lattice_constants == MonoclinicLatticeConstants(a, b, c, beta)
        assert unit_cell.lattice_system == LatticeSystem.MONOCLINIC
        assert unit_cell.symmetry == UnitCellSymmetry(
            centering=Centering.PRIMITIVE, symmetry_elements=symmetry_elements
        )

        assert unit_cell.a == a
        assert unit_cell.b == b
        assert unit_cell.c == c
        assert unit_cell.beta == beta
        assert unit_cell.centering == Centering.PRIMITIVE
        assert unit_cell.symmetry_elements == set(symmetry_elements)

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

    def test_to_julia(self):
        """
        Test `to_julia()`.
        """
        # --- Preparations

        a = 1
        b = 2
        c = 3
        beta = 0.1

        # ------ default centering and symmetry elements

        unit_cell = MonoclinicUnitCell(a, b, c, beta)
        unit_cell_jl = unit_cell.to_julia()

        assert self.jl.isa(unit_cell_jl, self.jl.UnitCell)
        assert self.jl.isa(unit_cell_jl, self.jl.MonoclinicUnitCell)
        assert unit_cell_jl.lattice_constants.a == a
        assert unit_cell_jl.lattice_constants.b == b
        assert unit_cell_jl.lattice_constants.c == c
        assert unit_cell_jl.lattice_constants.β == beta
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
        unit_cell = MonoclinicUnitCell(a, b, c, beta, centering=Centering.PRIMITIVE)
        unit_cell_jl = unit_cell.to_julia()

        assert self.jl.isa(unit_cell_jl, self.jl.UnitCell)
        assert self.jl.isa(unit_cell_jl, self.jl.MonoclinicUnitCell)
        assert unit_cell_jl.lattice_constants.a == a
        assert unit_cell_jl.lattice_constants.b == b
        assert unit_cell_jl.lattice_constants.c == c
        assert unit_cell_jl.lattice_constants.β == beta
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
        unit_cell = MonoclinicUnitCell(a, b, c, beta, centering=Centering.BASE)
        unit_cell_jl = unit_cell.to_julia()

        assert self.jl.isa(unit_cell_jl, self.jl.UnitCell)
        assert self.jl.isa(unit_cell_jl, self.jl.MonoclinicUnitCell)
        assert unit_cell_jl.lattice_constants.a == a
        assert unit_cell_jl.lattice_constants.b == b
        assert unit_cell_jl.lattice_constants.c == c
        assert unit_cell_jl.lattice_constants.β == beta
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
        unit_cell = MonoclinicUnitCell(a, b, c, beta, centering=Centering.BODY)
        unit_cell_jl = unit_cell.to_julia()

        assert self.jl.isa(unit_cell_jl, self.jl.UnitCell)
        assert self.jl.isa(unit_cell_jl, self.jl.MonoclinicUnitCell)
        assert unit_cell_jl.lattice_constants.a == a
        assert unit_cell_jl.lattice_constants.b == b
        assert unit_cell_jl.lattice_constants.c == c
        assert unit_cell_jl.lattice_constants.β == beta
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
        unit_cell = MonoclinicUnitCell(a, b, c, beta, centering=Centering.FACE)
        unit_cell_jl = unit_cell.to_julia()

        assert self.jl.isa(unit_cell_jl, self.jl.UnitCell)
        assert self.jl.isa(unit_cell_jl, self.jl.MonoclinicUnitCell)
        assert unit_cell_jl.lattice_constants.a == a
        assert unit_cell_jl.lattice_constants.b == b
        assert unit_cell_jl.lattice_constants.c == c
        assert unit_cell_jl.lattice_constants.β == beta
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
        unit_cell = MonoclinicUnitCell(
            a, b, c, beta, symmetry_elements=symmetry_elements
        )
        unit_cell_jl = unit_cell.to_julia()

        assert self.jl.isa(unit_cell_jl, self.jl.UnitCell)
        assert self.jl.isa(unit_cell_jl, self.jl.MonoclinicUnitCell)
        assert unit_cell_jl.lattice_constants.a == a
        assert unit_cell_jl.lattice_constants.b == b
        assert unit_cell_jl.lattice_constants.c == c
        assert unit_cell_jl.lattice_constants.β == beta
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
        b = 2
        c = 3
        beta = 0.1

        # --- Tests

        # ------ default centering and symmetry elements

        unit_cell_jl = self.jl.MonoclinicUnitCell(a, b, c, beta)

        unit_cell = MonoclinicUnitCell.from_julia(unit_cell_jl)
        assert unit_cell == MonoclinicUnitCell(a, b, c, beta)

        # ------ non-default centering

        # centering = primitive
        unit_cell_jl = self.jl.MonoclinicUnitCell(
            a, b, c, beta, centering=self.jl.primitive_centering
        )
        unit_cell = MonoclinicUnitCell.from_julia(unit_cell_jl)
        assert unit_cell == MonoclinicUnitCell(
            a, b, c, beta, centering=Centering.PRIMITIVE
        )

        # centering = base
        unit_cell_jl = self.jl.MonoclinicUnitCell(
            a, b, c, beta, centering=self.jl.base_centering
        )
        unit_cell = MonoclinicUnitCell.from_julia(unit_cell_jl)
        assert unit_cell == MonoclinicUnitCell(a, b, c, beta, centering=Centering.BASE)

        # centering = body
        unit_cell_jl = self.jl.MonoclinicUnitCell(
            a, b, c, beta, centering=self.jl.body_centering
        )
        unit_cell = MonoclinicUnitCell.from_julia(unit_cell_jl)
        assert unit_cell == MonoclinicUnitCell(a, b, c, beta, centering=Centering.BODY)

        # centering = face
        unit_cell_jl = self.jl.MonoclinicUnitCell(
            a, b, c, beta, centering=self.jl.face_centering
        )
        unit_cell = MonoclinicUnitCell.from_julia(unit_cell_jl)
        assert unit_cell == MonoclinicUnitCell(a, b, c, beta, centering=Centering.FACE)

        # ------ non-default symmetry elements

        symmetry_elements = [GlidePlane("1,0,0", "0,1,0"), ScrewAxis("1,0,0", 3, 2)]
        symmetry_elements_jl = self.jl.Vector(
            [element.to_julia() for element in symmetry_elements]
        )
        unit_cell_jl = self.jl.MonoclinicUnitCell(
            a,
            b,
            c,
            beta,
            symmetry_elements=symmetry_elements_jl,
        )

        unit_cell = MonoclinicUnitCell.from_julia(unit_cell_jl)
        assert unit_cell == MonoclinicUnitCell(
            a, b, c, beta, symmetry_elements=symmetry_elements
        )

    def test_from_julia_invalid_arguments(self):
        """
        Test `from_julia()`.
        """
        # --- Tests

        # unit_cell_jl is not a MonoclinicUnitCell
        unit_cell_jl_invalid = self.jl.CubicUnitCell(1)

        with pytest.raises(ValueError) as exception_info:
            MonoclinicUnitCell.from_julia(unit_cell_jl_invalid)

        expected_error = (
            "`unit_cell_jl` must be a Julia `MonoclinicUnitCell` object. "
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
        b = 2
        c = 3
        beta = 0.1

        # --- Tests

        # ------ default centering and symmetry

        unit_cell = MonoclinicUnitCell(a, b, c, beta)
        assert str(unit_cell) == (
            f"MonoclinicUnitCell(a={a},b={b},c={c},beta={beta},"
            "centering=primitive,"
            "symmetry_elements=[])"
        )

        # ------ non-default centering

        # centering = primitive
        unit_cell = MonoclinicUnitCell(a, b, c, beta, centering=Centering.PRIMITIVE)
        assert str(unit_cell) == (
            f"MonoclinicUnitCell(a={a},b={b},c={c},beta={beta},"
            "centering=primitive,"
            "symmetry_elements=[])"
        )

        # centering = base
        unit_cell = MonoclinicUnitCell(a, b, c, beta, centering=Centering.BASE)
        assert str(unit_cell) == (
            f"MonoclinicUnitCell(a={a},b={b},c={c},beta={beta},"
            "centering=base,"
            "symmetry_elements=[])"
        )

        # centering = body
        unit_cell = MonoclinicUnitCell(a, b, c, beta, centering=Centering.BODY)
        assert str(unit_cell) == (
            f"MonoclinicUnitCell(a={a},b={b},c={c},beta={beta},"
            "centering=body,"
            "symmetry_elements=[])"
        )

        # centering = face
        unit_cell = MonoclinicUnitCell(a, b, c, beta, centering=Centering.FACE)
        assert str(unit_cell) == (
            f"MonoclinicUnitCell(a={a},b={b},c={c},beta={beta},"
            "centering=face,"
            "symmetry_elements=[])"
        )

        # ------ non-default symmetry elements

        symmetry_elements = [GlidePlane("1,0,0", "0,1,0"), ScrewAxis("1,0,0", 3, 2)]
        unit_cell = MonoclinicUnitCell(
            a, b, c, beta, symmetry_elements=symmetry_elements
        )
        assert str(unit_cell) == (
            f"MonoclinicUnitCell(a={a},b={b},c={c},beta={beta},"
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
        b = 2
        c = 3
        beta = 0.1

        # --- Tests

        # ------ types differ

        unit_cell_1 = MonoclinicUnitCell(a, b, c, beta)
        unit_cell_2 = TetragonalUnitCell(a, a + 1)
        assert unit_cell_1 != unit_cell_2

        # ------ lattice constants and symmetry are the same

        unit_cell_1 = MonoclinicUnitCell(a, b, c, beta)
        unit_cell_2 = MonoclinicUnitCell(a, b, c, beta)
        assert unit_cell_1 == unit_cell_2

        # ------ lattice constants differ

        # `a` values differ
        unit_cell_1 = MonoclinicUnitCell(a, b, c, beta)
        unit_cell_2 = MonoclinicUnitCell(a + 1, b, c, beta)
        assert unit_cell_1 != unit_cell_2

        # `b` values differ
        unit_cell_1 = MonoclinicUnitCell(a, b, c, beta)
        unit_cell_2 = MonoclinicUnitCell(a, b + 2, c, beta)
        assert unit_cell_1 != unit_cell_2

        # `c` values differ
        unit_cell_1 = MonoclinicUnitCell(a, b, c, beta)
        unit_cell_2 = MonoclinicUnitCell(a, b, c + 3, beta)
        assert unit_cell_1 != unit_cell_2

        # `beta` values differ
        unit_cell_1 = MonoclinicUnitCell(a, b, c, beta)
        unit_cell_2 = MonoclinicUnitCell(a, b, c, beta + 0.1)
        assert unit_cell_1 != unit_cell_2

        # ------ symmetry differs

        # centerings differ
        unit_cell_1 = MonoclinicUnitCell(a, b, c, beta, centering="primitive")
        unit_cell_2 = MonoclinicUnitCell(a, b, c, beta, centering="body")
        assert unit_cell_1 != unit_cell_2

        # symmetry elements differ
        symmetry_elements = [GlidePlane("1,0,0", "0,1,0"), ScrewAxis("1,0,0", 3, 2)]
        unit_cell_1 = MonoclinicUnitCell(
            a, b, c, beta, symmetry_elements=symmetry_elements
        )
        unit_cell_2 = MonoclinicUnitCell(
            a, b, c, beta, symmetry_elements=symmetry_elements[0:-1]
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
        b = 2
        c = 3
        beta = 0.1

        # --- Tests

        # ------ types differ

        unit_cell_1 = MonoclinicUnitCell(a, b, c, beta)
        unit_cell_2 = TetragonalUnitCell(a, a + 1)
        assert not unit_cell_1.isclose(unit_cell_2)

        # ------ `a`

        # `a` values differ by less than tolerance
        unit_cell_1 = MonoclinicUnitCell(a, b, c, beta)
        unit_cell_2 = MonoclinicUnitCell(a + 0.1, b, c, beta)
        assert unit_cell_1.isclose(unit_cell_2, atol=0.2)

        # `a` values differ by more than tolerance
        unit_cell_1 = MonoclinicUnitCell(a, b, c, beta)
        unit_cell_2 = MonoclinicUnitCell(a + 1, b, c, beta)
        assert not unit_cell_1.isclose(unit_cell_2, atol=0.2)

        # ------ `b`

        # `b` values differ by less than tolerance
        unit_cell_1 = MonoclinicUnitCell(a, b, c, beta)
        unit_cell_2 = MonoclinicUnitCell(a, b + 0.1, c, beta)
        assert unit_cell_1.isclose(unit_cell_2, atol=0.2)

        # `b` values differ by more than tolerance
        unit_cell_1 = MonoclinicUnitCell(a, b, c, beta)
        unit_cell_2 = MonoclinicUnitCell(a, b + 2, c, beta)
        assert not unit_cell_1.isclose(unit_cell_2, atol=0.2)

        # ------ `c`

        # `c` values differ by less than tolerance
        unit_cell_1 = MonoclinicUnitCell(a, b, c, beta)
        unit_cell_2 = MonoclinicUnitCell(a, b, c - 0.05, beta)
        assert unit_cell_1.isclose(unit_cell_2, atol=0.2)

        # `c` values differ by more than tolerance
        unit_cell_1 = MonoclinicUnitCell(a, b, c, beta)
        unit_cell_2 = MonoclinicUnitCell(a, b, c + 3, beta)
        assert not unit_cell_1.isclose(unit_cell_2, atol=0.2)

        # ------ `beta`

        # `beta` values differ by less than tolerance
        unit_cell_1 = MonoclinicUnitCell(a, b, c, beta)
        unit_cell_2 = MonoclinicUnitCell(a, b, c, beta + 0.05)
        assert unit_cell_1.isclose(unit_cell_2, atol=0.2)

        # `beta` values differ by more than tolerance
        unit_cell_1 = MonoclinicUnitCell(a, b, c, beta)
        unit_cell_2 = MonoclinicUnitCell(a, b, c, beta + 0.25)
        assert not unit_cell_1.isclose(unit_cell_2, atol=0.2)

        # ------ symmetry differs

        # centerings differ
        unit_cell_1 = MonoclinicUnitCell(a, b, c, beta, centering=Centering.FACE)
        unit_cell_2 = MonoclinicUnitCell(a, b, c, beta, centering=Centering.BODY)
        assert not unit_cell_1.isclose(unit_cell_2)

        # symmetry elements differ
        symmetry_elements = [GlidePlane("1,0,0", "0,1,0"), ScrewAxis("1,0,0", 3, 2)]
        unit_cell_1 = MonoclinicUnitCell(
            a, b, c, beta, symmetry_elements=symmetry_elements
        )
        unit_cell_2 = MonoclinicUnitCell(
            a, b, c, beta, symmetry_elements=symmetry_elements[0:-1]
        )
        assert not unit_cell_1.isclose(unit_cell_2)
