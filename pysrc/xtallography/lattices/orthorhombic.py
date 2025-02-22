# Copyright (c) 2024 Velexi Corporation
#
# Use of this software is governed by the Business Source License included
# in the LICENSE file and at www.mariadb.com/bsl11.
#
# Change Date: Four years from the date the Licensed Work is published
#
# Change License: Mozilla Public License 2.0
"""
The orthorhombic module defines classes and methods specific to orthorhombic lattices.
"""

# --- Imports

# Local packages/modules
from .. import _JL
from .core import LatticeSystem, Centering, UnitCell


# --- Classes


class OrthorhombicUnitCell(UnitCell):
    """
    Lattice constants for a orthorhombic unit cell
    """

    # --- Initializer

    def __init__(
        self, a: float, b: float, c: float, centering: Centering = Centering.PRIMITIVE
    ):
        """
        Initialize OrthorhombicUnitCell object.

        Parameters
        ----------
        `a`, `b`, `c`: lattice constants
        """
        # Check arguments
        if a <= 0:
            raise ValueError(f"`a` must be positive. (a={a})")

        if b <= 0:
            raise ValueError(f"`b` must be positive. (b={b})")

        if c <= 0:
            raise ValueError(f"`c` must be positive. (c={c})")

        # Initialize parent class
        super().__init__(LatticeSystem.ORTHORHOMBIC, centering=centering)

        # Initialize lattice constants
        self._a = a
        self._b = b
        self._c = c

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

    # --- Methods

    def to_julia(self):
        """
        Convert OrthorhombicUnitCell object to a Julia struct.
        """
        return _JL.OrthorhombicLatticeConstants(self.a, self.b, self.c)
