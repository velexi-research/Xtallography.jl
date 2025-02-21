# Copyright (c) 2024 Velexi Corporation
#
# Use of this software is governed by the Business Source License included
# in the LICENSE file and at www.mariadb.com/bsl11.
#
# Change Date: Four years from the date the Licensed Work is published
#
# Change License: Mozilla Public License 2.0
"""
The triclinic module defines classes and methods specific to triclinic lattices.
"""

# --- Imports

# Standard library
import math

# Local packages/modules
from .. import jl
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
        Convert TriclinicUnitCell object to a Julia struct.
        """
        return jl.TriclinicLatticeConstants(
            self.a, self.b, self.c, self.alpha, self.beta, self.gamma
        )
