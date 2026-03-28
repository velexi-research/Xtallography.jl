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
UnitCellSymmetry class
"""

# --- Imports

# Standard library
from dataclasses import dataclass
import copy
from typing import Optional

# Local packages/modules
from .. import _JL

from xtallography.symmetry import Centering, SymmetryElement


# --- Classes


@dataclass(frozen=True)
class UnitCellSymmetry:
    """
    Class representing the symmetry of a unit cell

    Fields
    ------
    * `centering` (Centering): centering

    * `symmetry_elements` (set): set of symmetry elements
    """

    # --- Fields

    centering: Centering
    symmetry_elements: set

    # --- Initializer

    def __init__(
        self,
        centering: Optional[Centering] = Centering.PRIMITIVE,
        symmetry_elements: Optional[set, list] = copy.deepcopy([]),
    ):
        """
        Initialize UnitCellSymmetry object.

        Parameters
        ----------
        `centering`: centering

        `symmetry_elements`: set of symmetry elements
        """
        # --- Check arguments

        # Check that `symmetry_elements` contains only SymmetryElement objects
        if not all(
            [isinstance(element, SymmetryElement) for element in symmetry_elements]
        ):
            raise ValueError(
                "`symmetry_elements` contains elements that are not SymmetryElements "
                f"objects. (symmetry_elements={symmetry_elements})"
            )

        # --- Initialize field values

        # for frozen DataClasses, field values cannot be set directly
        object.__setattr__(self, "centering", centering)
        object.__setattr__(
            self, "symmetry_elements", frozenset(copy.copy(symmetry_elements))
        )

    # --- Methods

    def to_julia(self):
        """
        Convert Python UnitCellSymmetry object to a Julia UnitCellSymmetry object.
        """

        return _JL.UnitCellSymmetry(
            self.centering.to_julia(),
            symmetry_elements=_JL.Set(
                [element.to_julia() for element in self.symmetry_elements]
            ),
        )

    @classmethod
    def from_julia(cls, unit_cell_symmetry_jl: _JL.UnitCellSymmetry):
        """
        Convert a Julia UnitCellSymmetry object to a Python UnitCellSymmetry object.

        Parameters
        ----------
        unit_cell_symmetry_jl: Julia UnitCellSymmetry object

        Return value
        ------------
        unit_cell_symmetry: Python UnitCellSymmetry object
        """
        # --- Check arguments

        if not _JL.isa(unit_cell_symmetry_jl, _JL.UnitCellSymmetry):
            raise ValueError(
                "`unit_cell_symmetry_jl` must be a Julia `UnitCellSymmetry` object. "
                f"(unit_cell_symmetry_jl={unit_cell_symmetry_jl})."
            )

        # --- Convert unit_cell_symmetry_jl to a Julia UnitCellSymmetry object

        centering = Centering.from_julia(unit_cell_symmetry_jl.centering)
        symmetry_elements = set(
            [
                SymmetryElement.from_julia(element)
                for element in unit_cell_symmetry_jl.symmetry_elements
            ]
        )

        return UnitCellSymmetry(
            centering=centering, symmetry_elements=symmetry_elements
        )


"""
PRIMITIVE_UNIT_CELL_SYMMETRY

`UnitCellSymmetry` object representing a primitive unit cell with no additional symmetry
elements.
"""
PRIMITIVE_UNIT_CELL_SYMMETRY = UnitCellSymmetry()
