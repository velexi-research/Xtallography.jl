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
The hexagonal module defines classes and methods specific to hexagonal lattices.
"""

# --- Imports

# Local packages/modules
from .. import _JL
from .core import LatticeSystem, UnitCell


# --- Classes


class HexagonalUnitCell(UnitCell):
    """
    Lattice constants for a hexagonal unit cell
    """

    # --- Initializer

    def __init__(self, a: float, c: float):
        """
        Initialize HexagonalUnitCell object.

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
        super().__init__(LatticeSystem.HEXAGONAL)

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
        Convert HexagonalUnitCell object to a Julia struct.
        """
        return _JL.UnitCell(
            _JL.HexagonalLatticeConstants(self.a, self.c), _JL.primitive
        )

    @classmethod
    def from_julia(cls, unit_cell_jl: _JL.UnitCell):
        """
        Convert a Julia UnitCell struct to a HexagonalUnitCell object.
        """
        # Check arguments
        if not _JL.isa(unit_cell_jl, _JL.UnitCell):
            raise ValueError(
                "`unit_cell_jl` must be a Julia `UnitCell` struct. "
                f"(unit_cell_jl={unit_cell_jl})."
            )

        # Convert unit_cell_jl to a HexagonalUnitCell object
        unit_cell = HexagonalUnitCell(
            unit_cell_jl.lattice_constants.a,
            unit_cell_jl.lattice_constants.c,
        )

        return unit_cell

    def __repr__(self):
        """
        Return string representation of UnitCell.
        """
        return f"HexagonalUnitCell(a={self.a},c={self.c})"
