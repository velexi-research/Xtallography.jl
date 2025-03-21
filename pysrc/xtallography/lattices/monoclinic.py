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
The monoclinic module defines classes and methods specific to monoclinic lattices.
"""

# --- Imports

# Standard library
import math

# Local packages/modules
from .. import _JL
from .core import LatticeSystem, Centering, UnitCell


# --- Contants

MONOCLINIC_MIN_ANGLE = 0
MONOCLINIC_MAX_ANGLE = math.pi


# --- Classes


class MonoclinicUnitCell(UnitCell):
    """
    Lattice constants for a monoclinic unit cell
    """

    # --- Initializer

    def __init__(
        self,
        a: float,
        b: float,
        c: float,
        beta: float,
        centering: Centering = Centering.PRIMITIVE,
    ):
        """
        Initialize MonoclinicUnitCell object.

        Parameters
        ----------
        `a`, `b`, `c`, `beta`: lattice constants
        """
        # Check arguments
        if a <= 0:
            raise ValueError(f"`a` must be positive. (a={a})")

        if b <= 0:
            raise ValueError(f"`b` must be positive. (b={b})")

        if c <= 0:
            raise ValueError(f"`c` must be positive. (c={c})")

        if beta <= MONOCLINIC_MIN_ANGLE or beta >= MONOCLINIC_MAX_ANGLE:
            raise ValueError(f"`beta` must lie in the interval (0, pi). (beta={beta})")

        # Initialize parent class
        super().__init__(LatticeSystem.MONOCLINIC, centering=centering)

        # Initialize lattice constants
        self._a = a
        self._b = b
        self._c = c
        self._beta = beta

    # --- Properties

    @property
    def a(self):
        """
        Return `a`.
        """
        return self._a

    @property
    def b(self):
        """
        Return `b`.
        """
        return self._b

    @property
    def c(self):
        """
        Return `c`.
        """
        return self._c

    @property
    def beta(self):
        """
        Return `beta`.
        """
        return self._beta

    # --- Methods

    def to_julia(self):
        """
        Convert MonoclinicUnitCell object to a Julia UnitCell object.
        """
        return _JL.UnitCell(
            _JL.MonoclinicLatticeConstants(self.a, self.b, self.c, self.beta),
            self.centering.to_julia(),
        )

    @classmethod
    def from_julia(cls, unit_cell_jl: _JL.UnitCell):
        """
        Convert a Julia UnitCell object to a MonoclinicUnitCell object.
        """
        # --- Check arguments

        if not _JL.isa(unit_cell_jl, _JL.UnitCell):
            raise ValueError(
                "`unit_cell_jl` must be a Julia `UnitCell` object. "
                f"(unit_cell_jl={unit_cell_jl})."
            )

        if not _JL.isa(unit_cell_jl.lattice_constants, _JL.MonoclinicLatticeConstants):
            raise ValueError(
                "`unit_cell_jl` must be a Julia `UnitCell` object for monoclinic "
                f"unit cell. (unit_cell_jl={unit_cell_jl})."
            )

        # --- Convert unit_cell_jl to a MonoclinicUnitCell object

        unit_cell = MonoclinicUnitCell(
            unit_cell_jl.lattice_constants.a,
            unit_cell_jl.lattice_constants.b,
            unit_cell_jl.lattice_constants.c,
            unit_cell_jl.lattice_constants.β,
            centering=Centering.from_julia(unit_cell_jl.centering),
        )

        return unit_cell

    def __repr__(self):
        """
        Return string representation of MonoclinicUnitCell.
        """
        return (
            f"MonoclinicUnitCell(a={self.a},b={self.b},c={self.c},beta={self.beta},"
            f"centering='{self.centering}')"
        )
