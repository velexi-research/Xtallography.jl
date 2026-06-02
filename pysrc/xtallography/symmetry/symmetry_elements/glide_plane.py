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
GlidePlane class
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
class GlidePlane(SymmetryElement):
    """
    Class representing a glide plane

    Fields
    ------
    * `glide` (tuple): glide direction

    * `normal` (tuple): normal to glide plane

    * `location` (tuple): a point on the glide plane in unit cell fractional coordinates
    """

    # --- Fields

    glide: tuple
    normal: tuple
    location: tuple

    # --- Initializer

    def __init__(self, glide: tuple, normal: tuple, location: Optional[tuple] = None):
        """
        Initialize GlidePlane object.

        Parameters
        ----------
        `glide`: glide direction

        `normal`: normal to glide plane

        `location`: a point on the glide plane in unit cell fractional coordinates
        """
        # --- Check arguments

        # check that `glide` is a 3-tuple
        if len(glide) != 3:
            raise ValueError(f"`glide` must be a 3-tuple (glide={glide})")

        # check that `glide` is not the zero vector
        if all(component == 0 for component in glide):
            raise ValueError(f"`glide` must be a nonzero vector (glide={glide})")

        # convert `glide` to Fraction objects
        try:
            glide = tuple(Fraction(x) for x in glide)
        except ValueError as error:
            raise ValueError(
                "all elements of `glide` must convertible to Fraction objects "
                f"(glide={glide}). "
                f"[caused by {type(error)}({error})]"
            )

        # check that `normal` is a 3-tuple
        if len(normal) != 3:
            raise ValueError(f"`normal` must be a 3-tuple (normal={normal})")

        # check that `normal` is not the zero vector
        if all(component == 0 for component in normal):
            raise ValueError(f"`normal` must be a nonzero vector (normal={normal})")

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

        # Check that glide direction and normal to glide plane are orthogonal
        if numpy.dot(glide, normal) != 0:
            raise ValueError(
                f"`glide` must be orthogonal to `normal` (glide={glide},normal={normal})"
            )

        # --- Initialize field values

        # for frozen DataClasses, field values cannot be set directly
        object.__setattr__(self, "glide", glide)
        object.__setattr__(self, "normal", normal)
        object.__setattr__(self, "location", location)

    # --- Methods

    def to_julia(self):
        """
        Convert Python GlidePlane object to a Julia GlidePlane object.
        """
        return _JL.GlidePlane(self.glide, self.normal, self.location)

    @classmethod
    def from_julia(cls, glide_plane_jl: _JL.GlidePlane):
        """
        Convert a Julia GlidePlane object to a Python GlidePlane object.

        Parameters
        ----------
        glide_plane_jl: Julia GlidePlane object

        Return value
        ------------
        Python GlidePlane object
        """
        # Check arguments
        if not _JL.isa(glide_plane_jl, _JL.GlidePlane):
            raise ValueError(
                "`glide_plane_jl` must be a Julia `GlidePlane` object. "
                f"(glide_plane_jl={glide_plane_jl})."
            )

        # Convert glide_plane_jl to a GlidePlane object
        return GlidePlane(
            glide_plane_jl.glide, glide_plane_jl.normal, glide_plane_jl.location
        )

    def __repr__(self):
        """
        Return string representation of GlidePlane.
        """
        return (
            f"GlidePlane(glide={self.glide},normal={self.normal},"
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

        # Check that glides are the same

        if not all(self.glide[i] == other.glide[i] for i in range(3)):
            return False

        # Check that normals are in the same direction
        if numpy.dot(self.normal, other.normal) ** 2 != numpy.dot(
            self.normal, self.normal
        ) * numpy.dot(other.normal, other.normal):
            return False

        # Check that line through both locations lies a plane orthogonal to the plane
        # normal vectors
        delta = tuple(self.location[i] - other.location[i] for i in range(3))
        return numpy.dot(delta, self.normal) == 0
