# Copyright (c) 2024 Velexi Corporation
#
# Use of this software is governed by the Business Source License included
# in the LICENSE file and at www.mariadb.com/bsl11.
#
# Change Date: Four years from the date the Licensed Work is published
#
# Change License: Mozilla Public License 2.0
"""
The rhombohedral module defines classes and methods specific to rhombohedral lattices.
"""

# --- Imports

# Standard library
import math

# Local packages/modules
from . import jl
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
        Convert RhombohedralUnitCell object to a Julia struct.
        """
        return jl.RhombohedralLatticeConstants(self.a, self.alpha)
