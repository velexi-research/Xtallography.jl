# Copyright (c) 2024 Velexi Corporation
#
# Use of this software is governed by the Business Source License included
# in the LICENSE file and at www.mariadb.com/bsl11.
#
# Change Date: Four years from the date the Licensed Work is published
#
# Change License: Mozilla Public License 2.0
"""
The monoclinic module defines classes and methods specific to monoclinic lattices.
"""

# --- Imports

# Standard library
import math

# Local packages/modules
from .. import jl
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
        Convert MonoclinicUnitCell object to a Julia struct.
        """
        return jl.MonoclinicLatticeConstants(self.a, self.b, self.c, self.beta)
