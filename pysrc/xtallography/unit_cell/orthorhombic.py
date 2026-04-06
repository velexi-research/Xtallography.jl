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
The orthorhombic module defines classes and methods specific to orthorhombic lattices.
"""

# --- Imports

# Standard library
import copy
from typing import Optional

# Local packages/modules
from xtallography.symmetry import Centering, LatticeSystem

from .. import _JL
from .lattice_constants import OrthorhombicLatticeConstants
from .unit_cell import UnitCell
from .unit_cell_symmetry import UnitCellSymmetry


# --- Classes


class OrthorhombicUnitCell(UnitCell):
    """
    Class representing a orthorhombic unit cell
    """

    # --- Initializer

    def __init__(
        self,
        a: float,
        b: float,
        c: float,
        centering: Optional[Centering] = Centering.PRIMITIVE,
        symmetry_elements: Optional[set] = copy.deepcopy(set()),
    ):
        """
        Initialize OrthorhombicUnitCell object.

        Parameters
        ----------
        `a`, `b`, `c`: lattice constants

        `centering`: centering of unit cell

        `symmetry_elements`: symmetry elements that unit cell is invariant under
        """
        # --- Check arguments

        if a <= 0:
            raise ValueError(f"`a` must be positive. (a={a})")

        if b <= 0:
            raise ValueError(f"`b` must be positive. (b={b})")

        if c <= 0:
            raise ValueError(f"`c` must be positive. (c={c})")

        # --- Initialize parent class

        lattice_constants = OrthorhombicLatticeConstants(a, b, c)
        symmetry = UnitCellSymmetry(centering, symmetry_elements)
        super().__init__(LatticeSystem.ORTHORHOMBIC, lattice_constants, symmetry)

    # --- Properties

    @property
    def a(self):
        """
        Return `a`.
        """
        return self.lattice_constants.a

    @property
    def b(self):
        """
        Return `b`.
        """
        return self.lattice_constants.b

    @property
    def c(self):
        """
        Return `c`.
        """
        return self.lattice_constants.c

    # --- Methods

    def to_julia(self):
        """
        Convert OrthorhombicUnitCell object to a Julia UnitCell object.
        """
        return _JL.OrthorhombicUnitCell(
            self.a,
            self.b,
            self.c,
            centering=self.centering.to_julia(),
            symmetry_elements=_JL.Vector(
                [element.to_julia() for element in self.symmetry_elements]
            ),
        )

    @classmethod
    def from_julia(cls, unit_cell_jl: _JL.UnitCell):
        """
        Convert a Julia UnitCell object to a OrthorhombicUnitCell object.
        """
        # --- Check arguments

        if not _JL.isa(unit_cell_jl, _JL.OrthorhombicUnitCell):
            raise ValueError(
                "`unit_cell_jl` must be a Julia `OrthorhombicUnitCell` object. "
                f"(unit_cell_jl={unit_cell_jl})."
            )

        # --- Convert unit_cell_jl to a OrthorhombicUnitCell object

        unit_cell_symmetry = UnitCellSymmetry.from_julia(unit_cell_jl.symmetry)
        unit_cell = OrthorhombicUnitCell(
            unit_cell_jl.lattice_constants.a,
            unit_cell_jl.lattice_constants.b,
            unit_cell_jl.lattice_constants.c,
            centering=unit_cell_symmetry.centering,
            symmetry_elements=unit_cell_symmetry.symmetry_elements,
        )

        return unit_cell

    def __repr__(self):
        """
        Return string representation of OrthorhombicUnitCell.
        """
        return (
            f"OrthorhombicUnitCell(a={self.a},b={self.b},c={self.c},"
            f"centering={self.centering},"
            "symmetry_elements=["
            f"{','.join(sorted([str(item) for item in self.symmetry_elements]))}"
            "])"
        )
