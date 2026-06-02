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
ScrewAxis class
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
class ScrewAxis(SymmetryElement):
    """
    Class representing a screw axis

    Fields
    ------
    * `n` (int): rotation order

    * `m` (int): number of translation steps of size 1/n following rotation by 2π/n

    * `direction` (tuple): direction of rotation axis

    * `location` (tuple): a point on the rotation axis in unit cell fractional coordinates
    """

    # --- Fields

    n: int
    m: int
    direction: tuple
    location: tuple

    # --- Initializer

    def __init__(
        self, n: int, m: int, direction: tuple, location: Optional[tuple] = None
    ):
        """
        Initialize ScrewAxis object.

        Parameters
        ----------
        `n`: rotation order

        `m`: number of translation steps of size 1/n following rotation by 2π/n

        `direction`: direction of rotation axis

        `location`: a point on the rotation axis in unit cell fractional coordinates
        """
        # --- Check arguments

        if n <= 0:
            raise ValueError(f"`n` must be positive (n={n})")

        if m <= 0:
            raise ValueError(f"`m` must be positive (m={m})")

        if m >= n:
            raise ValueError(f"`m` must be less than `n` (n={n},m={m})")

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

        # set default value for `location`
        if location is None:
            location = (Fraction(0), Fraction(0), Fraction(0))

        else:
            # check that `location` is a 3-tuple
            if len(location) != 3:
                raise ValueError(f"`location` must be a 3-tuple (location={location})")

            # convert `location` to Fraction objects
            try:
                location = tuple(Fraction(x) for x in location)
            except ValueError as error:
                raise ValueError(
                    "all elements of `location` must convertible to Fraction objects "
                    f"(location={location}). "
                    f"[caused by {type(error)}({error})]"
                )

        # --- Initialize field values

        # for frozen DataClasses, field values cannot be set directly
        object.__setattr__(self, "n", n)
        object.__setattr__(self, "m", m)
        object.__setattr__(self, "direction", direction)
        object.__setattr__(self, "location", location)

    # --- Methods

    def to_julia(self):
        """
        Convert Python ScrewAxis object to a Julia ScrewAxis object.
        """
        return _JL.ScrewAxis(self.n, self.m, self.direction, self.location)

    @classmethod
    def from_julia(cls, screw_axis_jl: _JL.ScrewAxis):
        """
        Convert a Julia ScrewAxis object to a Python ScrewAxis object.

        Parameters
        ----------
        screw_axis_jl: Julia ScrewAxis object

        Return value
        ------------
        Python ScrewAxis object
        """
        # Check arguments
        if not _JL.isa(screw_axis_jl, _JL.ScrewAxis):
            raise ValueError(
                "`screw_axis_jl` must be a Julia `ScrewAxis` object. "
                f"(screw_axis_jl={screw_axis_jl})."
            )

        # Convert screw_axis_jl to a ScrewAxis object
        return ScrewAxis(
            screw_axis_jl.n,
            screw_axis_jl.m,
            screw_axis_jl.direction,
            screw_axis_jl.location,
        )

    def __repr__(self):
        """
        Return string representation of ScrewAxis.
        """
        return (
            f"ScrewAxis(n={self.n},m={self.m},direction={self.direction},"
            f"location={self.location})"
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

        # Check that the number of translation steps are equal
        if self.m != other.m:
            return False

        # Check that directions are the same
        if numpy.dot(self.direction, other.direction) ** 2 != numpy.dot(
            self.direction, self.direction
        ) * numpy.dot(other.direction, other.direction):
            return False

        # Check that line through both locations is in the same direction as both lines
        delta = tuple(self.location[i] - other.location[i] for i in range(3))
        return numpy.dot(delta, self.direction) ** 2 == numpy.dot(
            delta, delta
        ) * numpy.dot(self.direction, self.direction)
