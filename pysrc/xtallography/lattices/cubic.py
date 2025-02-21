# Copyright (c) 2024 Velexi Corporation
#
# Use of this software is governed by the Business Source License included
# in the LICENSE file and at www.mariadb.com/bsl11.
#
# Change Date: Four years from the date the Licensed Work is published
#
# Change License: Mozilla Public License 2.0
"""
The cubic module defines classes and methods specific to cubic lattices.
"""

# --- Imports

# Local packages/modules
from .. import jl
from .core import LatticeSystem, Centering, UnitCell


# --- Classes


class CubicUnitCell(UnitCell):
    """
    Lattice constants for a cubic unit cell
    """

    # --- Initializer

    def __init__(self, a: float, centering: Centering = Centering.PRIMITIVE):
        """
        Initialize CubicUnitCell object.

        Parameters
        ----------
        `a`: lattice constant
        """
        # Check arguments
        if a <= 0:
            raise ValueError(f"`a` must be positive. (a={a})")

        # Initialize parent class
        super().__init__(LatticeSystem.CUBIC, centering=centering)

        # Initialize lattice constants
        self._a = a

    # --- Properties

    @property
    def a(self):
        """
        Return `a`.
        """
        return self._a

    # --- Methods

    def to_julia(self):
        """
        Convert CubicUnitCell object to a Julia struct.
        """
        return jl.CubicLatticeConstants(self.a)
