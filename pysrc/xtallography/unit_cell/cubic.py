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
The cubic module defines classes and methods specific to cubic lattices.
"""

# --- Imports

# Standard library
import copy
from typing import Optional

# Local packages/modules
from xtallography.symmetry import Centering, LatticeSystem

from .. import _JL
from .lattice_constants import CubicLatticeConstants
from .unit_cell import UnitCell
from .unit_cell_symmetry import UnitCellSymmetry


# --- Classes


class CubicUnitCell(UnitCell):
    """
    Class representing a cubic unit cell
    """

    # --- Initializer

    def __init__(
        self,
        a: float,
        centering: Optional[Centering] = Centering.PRIMITIVE,
        symmetry_elements: Optional[set] = copy.deepcopy(set()),
    ):
        """
        Initialize CubicUnitCell object.

        Parameters
        ----------
        `a`: lattice constant

        `centering`: centering of unit cell

        `symmetry_elements`: symmetry elements that unit cell is invariant under
        """
        # --- Check arguments

        if a <= 0:
            raise ValueError(f"`a` must be positive. (a={a})")

        # --- Initialize parent class

        lattice_constants = CubicLatticeConstants(a)
        symmetry = UnitCellSymmetry(centering, symmetry_elements)
        super().__init__(LatticeSystem.CUBIC, lattice_constants, symmetry)

    # --- Properties

    @property
    def a(self):
        """
        Return `a`.
        """
        return self.lattice_constants.a

    # --- Methods

    def to_julia(self):
        """
        Convert CubicUnitCell object to a Julia UnitCell object.
        """
        return _JL.CubicUnitCell(
            self.a,
            centering=self.centering.to_julia(),
            symmetry_elements=_JL.Vector(
                [element.to_julia() for element in self.symmetry_elements]
            ),
        )

    @classmethod
    def from_julia(cls, unit_cell_jl: _JL.UnitCell):
        """
        Convert a Julia UnitCell object to a CubicUnitCell object.
        """
        # --- Check arguments

        if not _JL.isa(unit_cell_jl, _JL.CubicUnitCell):
            raise ValueError(
                "`unit_cell_jl` must be a Julia `CubicUnitCell` object. "
                f"(unit_cell_jl={unit_cell_jl})."
            )

        # --- Convert unit_cell_jl to a CubicUnitCell object

        unit_cell_symmetry = UnitCellSymmetry.from_julia(unit_cell_jl.symmetry)
        unit_cell = CubicUnitCell(
            unit_cell_jl.lattice_constants.a,
            centering=unit_cell_symmetry.centering,
            symmetry_elements=unit_cell_symmetry.symmetry_elements,
        )

        return unit_cell

    def __repr__(self):
        """
        Return string representation of CubicUnitCell.
        """
        return (
            f"CubicUnitCell(a={self.a},"
            f"centering={self.centering},"
            "symmetry_elements=["
            f"{','.join(sorted([str(item) for item in self.symmetry_elements]))}"
            "])"
        )
