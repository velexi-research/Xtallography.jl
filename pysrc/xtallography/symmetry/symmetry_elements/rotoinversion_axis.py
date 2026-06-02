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
RotoinversionAxis class
"""

# --- Imports

# Standard library
from dataclasses import dataclass
from fractions import Fraction
from typing import Optional

# External packages
import numpy

# Local packages/modules
from .symmetry_element import SymmetryElement
from ... import _JL


# --- Classes


@dataclass(frozen=True)
class RotoinversionAxis(SymmetryElement):
    """
    Class representing a rotoinversion axis

    Fields
    ------
    * `n` (int): rotoinversion order

    * `direction` (tuple): direction of rotoinversion axis

    * `center` (tuple): location of inversion center in unit cell fractional coordinates
      coordinates
    """

    # --- Fields

    n: int
    direction: tuple
    center: tuple

    # --- Initializer

    def __init__(self, n: int, direction: tuple, center: Optional[tuple] = None):
        """
        Initialize RotoinversionAxis object.

        Parameters
        ----------
        `n`: rotoinversion order

        `direction`: direction of rotoinversion axis

        `center`: location of inversion center in unit cell fractional coordinates
                  coordinates
        """
        # --- Check arguments

        # check that n is positive
        if n <= 0:
            raise ValueError(f"`n` must be positive (n={n})")

        # check that `direction` is a 3-tuple
        if len(direction) != 3:
            raise ValueError(f"`direction` must be a 3-tuple (direction={direction})")

        # convert `direction` to Fraction objects
        try:
            direction = tuple(Fraction(x) for x in direction)
        except ValueError as error:
            raise ValueError(
                "all elements of `direction` must convertible to Fraction objects "
                f"(direction={direction}). "
                f"[caused by {type(error)}({error})]"
            )

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
        object.__setattr__(self, "n", n)
        object.__setattr__(self, "direction", direction)
        object.__setattr__(self, "center", center)

    # --- Methods

    def to_julia(self):
        """
        Convert Python RotoinversionAxis object to a Julia RotoinversionAxis object.
        """
        return _JL.RotoinversionAxis(self.n, self.direction, self.center)

    @classmethod
    def from_julia(cls, rotoinversion_axis_jl: _JL.RotoinversionAxis):
        """
        Convert a Julia RotoinversionAxis object to a Python RotoinversionAxis object.

        Parameters
        ----------
        rotoinversion_axis_jl: Julia RotoinversionAxis object

        Return value
        ------------
        Python RotoinversionAxis object
        """
        # Check arguments
        if not _JL.isa(rotoinversion_axis_jl, _JL.RotoinversionAxis):
            raise ValueError(
                "`rotoinversion_axis_jl` must be a Julia `RotoinversionAxis` object. "
                f"(rotoinversion_axis_jl={rotoinversion_axis_jl})."
            )

        # Convert rotoinversion_axis_jl to a RotoinversionAxis object
        return RotoinversionAxis(
            rotoinversion_axis_jl.n,
            rotoinversion_axis_jl.direction,
            rotoinversion_axis_jl.center,
        )

    def __repr__(self):
        """
        Return string representation of RotoinversionAxis.
        """
        return (
            f"RotoinversionAxis(n={self.n},direction={self.direction},"
            f"center={self.center})"
        )

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

        # Check that orders are equal
        if self.n != other.n:
            return False

        # Check centers are the same
        if not all(self.center[i] == other.center[i] for i in range(3)):
            return False

        # Check that directions are the same
        return numpy.dot(self.direction, other.direction) ** 2 == numpy.dot(
            self.direction, self.direction
        ) * numpy.dot(other.direction, other.direction)
