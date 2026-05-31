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
MirrorPlane class
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
class MirrorPlane(SymmetryElement):
    """
    Class representing a mirror plane

    Fields
    ------
    * `normal` (tuple): normal to mirror plane

    * `location` (tuple): a point on the mirror plane in unit cell fractional coordinates
    """

    # --- Fields

    normal: tuple
    location: tuple

    # --- Initializer

    def __init__(self, normal: tuple, location: Optional[tuple] = None):
        """
        Initialize MirrorPlane object.

        Parameters
        ----------
        `normal` (tuple): normal to mirror plane

        `location` (tuple): a point on the mirror plane in unit cell fractional coordinates
        """
        # --- Check arguments

        # check that `normal` is a 3-tuple
        if len(normal) != 3:
            raise ValueError(f"`normal` must be a 3-tuple (normal={normal})")

        # convert `normal` to Fraction objects
        try:
            normal = tuple(Fraction(x) for x in normal)
        except ValueError as error:
            raise ValueError(
                "all elements of `normal` must convertible to Fraction objects "
                f"(normal={normal}). "
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
        object.__setattr__(self, "normal", normal)
        object.__setattr__(self, "location", location)

    # --- Methods

    def to_julia(self):
        """
        Convert Python MirrorPlane object to a Julia MirrorPlane object.
        """
        return _JL.MirrorPlane(self.normal, self.location)

    @classmethod
    def from_julia(cls, mirror_plane_jl: _JL.MirrorPlane):
        """
        Convert a Julia MirrorPlane object to a Python MirrorPlane object.

        Parameters
        ----------
        mirror_plane_jl: Julia MirrorPlane object

        Return value
        ------------
        Python MirrorPlane object
        """
        # Check arguments
        if not _JL.isa(mirror_plane_jl, _JL.MirrorPlane):
            raise ValueError(
                "`mirror_plane_jl` must be a Julia `MirrorPlane` object. "
                f"(mirror_plane_jl={mirror_plane_jl})."
            )

        # Convert mirror_plane_jl to a MirrorPlane object
        return MirrorPlane(mirror_plane_jl.normal, mirror_plane_jl.location)

    def __repr__(self):
        """
        Return string representation of MirrorPlane.
        """
        return f"MirrorPlane(normal={self.normal},location={self.location})"

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

        # Check that normals are the same
        if numpy.dot(self.normal, other.normal) ** 2 != numpy.dot(
            self.normal, self.normal
        ) * numpy.dot(other.normal, other.normal):
            return False

        # Check that line through both locations lies a plane orthogonal to the plane
        # normal vectors
        delta = tuple(self.location[i] - other.location[i] for i in range(3))
        return numpy.dot(delta, self.normal) == 0
