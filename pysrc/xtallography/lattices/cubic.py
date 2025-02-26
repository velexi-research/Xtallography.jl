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
The cubic module defines classes and methods specific to cubic lattices.
"""

# --- Imports

# Local packages/modules
from .. import _JL
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
        return _JL.CubicLatticeConstants(self.a)
