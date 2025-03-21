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
The triclinic module defines classes and methods specific to triclinic lattices.
"""

# --- Imports

# Standard library
import math

# Local packages/modules
from .. import _JL
from .core import LatticeSystem, UnitCell


# --- Classes


class TriclinicUnitCell(UnitCell):
    """
    Lattice constants for a triclinic unit cell
    """

    # --- Initializer

    def __init__(
        self, a: float, b: float, c: float, alpha: float, beta: float, gamma: float
    ):
        """
        Initialize TriclinicUnitCell object.

        Parameters
        ----------
        `a`, `b`, `c`, `alpha`, `beta`, `gamma`: lattice constants
        """
        # Check arguments
        if a <= 0:
            raise ValueError(f"`a` must be positive. (a={a})")

        if b <= 0:
            raise ValueError(f"`b` must be positive. (b={b})")

        if c <= 0:
            raise ValueError(f"`c` must be positive. (c={c})")

        if alpha <= 0 or alpha >= 2 * math.pi:
            raise ValueError(
                f"`alpha` must lie in the interval (0, 2 pi). (alpha={alpha})"
            )

        if beta <= 0 or beta >= 2 * math.pi:
            raise ValueError(
                f"`beta` must lie in the interval (0, 2 pi). (beta={beta})"
            )

        if gamma <= 0 or gamma >= 2 * math.pi:
            raise ValueError(
                f"`gamma` must lie in the interval (0, 2 pi). (gamma={gamma})"
            )

        # TODO: add valid angle check

        # Initialize parent class
        super().__init__(LatticeSystem.TRICLINIC)

        # Initialize lattice constants
        self._a = a
        self._b = b
        self._c = c
        self._alpha = alpha
        self._beta = beta
        self._gamma = gamma

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
    def alpha(self):
        """
        Return `alpha`.
        """
        return self._alpha

    @property
    def beta(self):
        """
        Return `beta`.
        """
        return self._beta

    @property
    def gamma(self):
        """
        Return `gamma`.
        """
        return self._gamma

    # --- Methods

    def to_julia(self):
        """
        Convert TriclinicUnitCell object to a Julia UnitCell object.
        """
        return _JL.UnitCell(
            _JL.TriclinicLatticeConstants(
                self.a, self.b, self.c, self.alpha, self.beta, self.gamma
            ),
            _JL.primitive,
        )

    @classmethod
    def from_julia(cls, unit_cell_jl: _JL.UnitCell):
        """
        Convert a Julia UnitCell object to a TriclinicUnitCell object.
        """
        # --- Check arguments

        if not _JL.isa(unit_cell_jl, _JL.UnitCell):
            raise ValueError(
                "`unit_cell_jl` must be a Julia `UnitCell` object. "
                f"(unit_cell_jl={unit_cell_jl})."
            )

        if not _JL.isa(unit_cell_jl.lattice_constants, _JL.TriclinicLatticeConstants):
            raise ValueError(
                "`unit_cell_jl` must be a Julia `UnitCell` object for triclinic "
                f"unit cell. (unit_cell_jl={unit_cell_jl})."
            )

        # --- Convert unit_cell_jl to a TriclinicUnitCell object

        unit_cell = TriclinicUnitCell(
            unit_cell_jl.lattice_constants.a,
            unit_cell_jl.lattice_constants.b,
            unit_cell_jl.lattice_constants.c,
            unit_cell_jl.lattice_constants.α,
            unit_cell_jl.lattice_constants.β,
            unit_cell_jl.lattice_constants.γ,
        )

        return unit_cell

    def __repr__(self):
        """
        Return string representation of TriclinicUnitCell.
        """
        return (
            f"TriclinicUnitCell(a={self.a},b={self.b},c={self.c},"
            f"alpha={self.alpha},beta={self.beta},gamma={self.gamma})"
        )
