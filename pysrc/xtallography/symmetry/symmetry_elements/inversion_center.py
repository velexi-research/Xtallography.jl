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
InversionCenter class
"""

# --- Imports

# Standard library
from dataclasses import dataclass
from fractions import Fraction
from typing import Optional

# Local packages/modules
from .symmetry_element import SymmetryElement
from ... import _JL


# --- Classes


@dataclass(frozen=True)
class InversionCenter(SymmetryElement):
    """
    Class representing a inversion center

    Fields
    ------
    * `center` (tuple): location of inversion center in unit cell fractional coordinates
    """

    # --- Fields

    center: tuple

    # --- Initializer

    def __init__(self, center: Optional[tuple] = None):
        """
        Initialize InversionCenter object.

        Parameters
        ----------
        `center`: location of inversion center in unit cell fractional coordinates
        """
        # --- Check arguments

        # set default value for `center`
        if center is None:
            center = (Fraction(0), Fraction(0), Fraction(0))

        else:
            # check that `center` is a 3-tuple
            if len(center) != 3:
                raise ValueError(f"`center` must be a 3-tuple (center={center})")

            # convert `center` to Fraction objects
            try:
                center = tuple(Fraction(x) for x in center)
            except ValueError as error:
                raise ValueError(
                    "all elements of `center` must convertible to Fraction objects "
                    f"(center={center}). "
                    f"[caused by {type(error)}({error})]"
                )

        # --- Initialize field values

        # for frozen DataClasses, field values cannot be set directly
        object.__setattr__(self, "center", center)

    # --- Methods

    def to_julia(self):
        """
        Convert Python InversionCenter object to a Julia InversionCenter object.
        """
        return _JL.InversionCenter(self.center)

    @classmethod
    def from_julia(cls, inversion_center_jl: _JL.InversionCenter):
        """
        Convert a Julia InversionCenter object to a Python InversionCenter object.

        Parameters
        ----------
        inversion_center_jl: Julia InversionCenter object

        Return value
        ------------
        Python InversionCenter object
        """
        # Check arguments
        if not _JL.isa(inversion_center_jl, _JL.InversionCenter):
            raise ValueError(
                "`inversion_center_jl` must be a Julia `InversionCenter` object. "
                f"(inversion_center_jl={inversion_center_jl})."
            )

        # Convert inversion_center_jl to a InversionCenter object
        return InversionCenter(inversion_center_jl.center)

    def __repr__(self):
        """
        Return string representation of InversionCenter.
        """
        return f"InversionCenter(center={self.center})"

    def __eq__(self, other):
        """
        Return True if `self` and `other`are identical symmetry elements; otherwise,
        return False.

        Parameters
        ----------
        other:  SymmetryElement object to compare against
        """
        if not isinstance(other, type(self)):
            return False

        # Check that centers are the same
        return all(self.center[i] == other.center[i] for i in range(3))
