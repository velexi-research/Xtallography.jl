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
The tetragonal module defines classes and methods specific to tetragonal lattices.
"""

# --- Imports

# Local packages/modules
from .. import _JL
from .core import LatticeSystem, Centering, UnitCell


# --- Classes


class TetragonalUnitCell(UnitCell):
    """
    Lattice constants for a tetragonal unit cell
    """

    # --- Initializer

    def __init__(self, a: float, c: float, centering: Centering = Centering.PRIMITIVE):
        """
        Initialize TetragonalUnitCell object.

        Parameters
        ----------
        `a`, `c`: lattice constants
        """
        # Check arguments
        if a <= 0:
            raise ValueError(f"`a` must be positive. (a={a})")

        if c <= 0:
            raise ValueError(f"`c` must be positive. (c={c})")

        # Initialize parent class
        super().__init__(LatticeSystem.TETRAGONAL, centering=centering)

        # Initialize lattice constants
        self._a = a
        self._c = c

    # --- Properties

    @property
    def a(self):
        """
        Return `a`.
        """
        return self._a

    @property
    def c(self):
        """
        Return `c`.
        """
        return self._c

    # --- Methods

    def to_julia(self):
        """
        Convert TetragonalUnitCell object to a Julia UnitCell object.
        """
        return _JL.UnitCell(
            _JL.TetragonalLatticeConstants(self.a, self.c), self.centering.to_julia()
        )

    @classmethod
    def from_julia(cls, unit_cell_jl: _JL.UnitCell):
        """
        Convert a Julia UnitCell object to a TetragonalUnitCell object.
        """
        # --- Check arguments

        if not _JL.isa(unit_cell_jl, _JL.UnitCell):
            raise ValueError(
                "`unit_cell_jl` must be a Julia `UnitCell` object. "
                f"(unit_cell_jl={unit_cell_jl})."
            )

        if not _JL.isa(unit_cell_jl.lattice_constants, _JL.TetragonalLatticeConstants):
            raise ValueError(
                "`unit_cell_jl` must be a Julia `UnitCell` object for tetragonal "
                f"unit cell. (unit_cell_jl={unit_cell_jl})."
            )

        # --- Convert unit_cell_jl to a TetragonalUnitCell object

        unit_cell = TetragonalUnitCell(
            unit_cell_jl.lattice_constants.a,
            unit_cell_jl.lattice_constants.c,
            centering=Centering.from_julia(unit_cell_jl.centering),
        )

        return unit_cell

    def __repr__(self):
        """
        Return string representation of TetragonalUnitCell.
        """
        return f"TetragonalUnitCell(a={self.a},c={self.c},centering='{self.centering}')"
