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
The rhombohedral module defines classes and methods specific to rhombohedral lattices.
"""

# --- Imports

# Standard library
import math

# Local packages/modules
from .. import _JL
from .core import LatticeSystem, UnitCell


# --- Classes


class RhombohedralUnitCell(UnitCell):
    """
    Lattice constants for a rhombohedral unit cell
    """

    # ---Initializer

    def __init__(self, a: float, alpha: float):
        """
        Initialize RhombohedralUnitCell object.

        Parameters
        ----------
        `a`, `alpha`: lattice constants
        """
        # Check arguments
        if a <= 0:
            raise ValueError(f"`a` must be positive. (a={a})")

        if alpha <= 0 or alpha >= 2 * math.pi:
            raise ValueError(
                f"`alpha` must lie in the interval (0, 2 pi). (alpha={alpha})"
            )

        # Initialize parent class
        super().__init__(LatticeSystem.RHOMBOHEDRAL)

        # Initialize lattice constants
        self._a = a
        self._alpha = alpha

    # --- Properties

    @property
    def a(self):
        """
        Return `a`.
        """
        return self._a

    @property
    def alpha(self):
        """
        Return `alpha`.
        """
        return self._alpha

    # --- Methods

    def to_julia(self):
        """
        Convert RhombohedralUnitCell object to a Julia UnitCell object.
        """
        return _JL.UnitCell(
            _JL.RhombohedralLatticeConstants(self.a, self.alpha), _JL.primitive
        )

    @classmethod
    def from_julia(cls, unit_cell_jl: _JL.UnitCell):
        """
        Convert a Julia UnitCell object to a RhombohedralUnitCell object.
        """
        # --- Check arguments

        if not _JL.isa(unit_cell_jl, _JL.UnitCell):
            raise ValueError(
                "`unit_cell_jl` must be a Julia `UnitCell` object. "
                f"(unit_cell_jl={unit_cell_jl})."
            )

        if not _JL.isa(
            unit_cell_jl.lattice_constants, _JL.RhombohedralLatticeConstants
        ):
            raise ValueError(
                "`unit_cell_jl` must be a Julia `UnitCell` object for rhombohedral "
                f"unit cell. (unit_cell_jl={unit_cell_jl})."
            )

        # --- Convert unit_cell_jl to a RhombohedralUnitCell object

        unit_cell = RhombohedralUnitCell(
            unit_cell_jl.lattice_constants.a,
            unit_cell_jl.lattice_constants.α,
        )

        return unit_cell

    def __repr__(self):
        """
        Return string representation of RhombohedralUnitCell.
        """
        return f"RhombohedralUnitCell(a={self.a},alpha={self.alpha})"
