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
Core types and functions that support lattice and unit cell computations
"""

# --- Imports

# Standard library
from abc import abstractmethod, ABC
from dataclasses import dataclass

# Local packages/modules
from .. import _JL


# --- Classes


class SymmetryElement(ABC):
    """
    Abstract base class for symmetry element classes
    """

    # --- Methods

    @abstractmethod
    def to_julia(self):
        """
        Convert Python SymmetryElement object to a Julia SymmetryElement object.
        """

    @classmethod
    def from_julia(cls, symmetry_element_jl: _JL.SymmetryElement):
        """
        Convert a Julia SymmetryElement object to a Python SymmetryElement object.

        Parameters
        ----------
        symmetry_element_jl: Julia SymmetryElement object

        Return value
        ------------
        Python SymmetryElement object
        """
        if _JL.isa(symmetry_element_jl, _JL.GlidePlane):
            return GlidePlane.from_julia(symmetry_element_jl)
        elif _JL.isa(symmetry_element_jl, _JL.ScrewAxis):
            return ScrewAxis.from_julia(symmetry_element_jl)

        raise ValueError(
            "`symmetry_element_jl` must be a Julia `GlidePlane` or `ScrewAxis` object. "
            f"(symmetry_element_jl={symmetry_element_jl})."
        )

    @abstractmethod
    def __repr__(self):
        """
        Return string representation of SymmetryElement.
        """


@dataclass(frozen=True)
class GlidePlane(SymmetryElement):
    """
    Class representing a glide plane

    Fields
    ------
    * `translation` (str): glide translation

    * `reflection_plane` (str): reflection plane
    """

    # --- Fields

    translation: str
    reflection_plane: str

    # --- Methods

    def to_julia(self):
        """
        Convert Python GlidePlane object to a Julia GlidePlane object.
        """
        return _JL.GlidePlane(self.translation, self.reflection_plane)

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
        return GlidePlane(glide_plane_jl.translation, glide_plane_jl.reflection_plane)

    def __repr__(self):
        """
        Return string representation of GlidePlane.
        """
        return (
            f"GlidePlane(translation='{self.translation}',"
            f"reflection_plane='{self.reflection_plane}')"
        )


@dataclass(frozen=True)
class ScrewAxis(SymmetryElement):
    """
    Class representing a screw axis

    Fields
    ------
    * `axis` (str): direction of rotation axis

    * `n` (int): rotation order

    * `m` (int): number of translation steps of size 1/n following rotation by 2π/n
    """

    # --- Fields

    axis: str
    n: int
    m: int

    # --- Initializer

    def __init__(self, axis: str, n: int, m: int):
        """
        Initialize ScrewAxis object.

        Parameters
        ----------
        `axis`: direction of rotation axis

        `n`: rotation order

        `m`: number of translation steps of size 1/n following rotation by 2π/n
        """
        # --- Check arguments

        if n <= 0:
            raise ValueError(f"`n` must be positive (n={n})")

        if m <= 0:
            raise ValueError(f"`m` must be positive (m={m})")

        if m >= n:
            raise ValueError(f"`m` must be less than `n` (n={n},m={m})")

        # --- Initialize field values

        # for frozen DataClasses, field values cannot be set directly
        object.__setattr__(self, "axis", axis)
        object.__setattr__(self, "n", n)
        object.__setattr__(self, "m", m)

    # --- Methods

    def to_julia(self):
        """
        Convert Python ScrewAxis object to a Julia ScrewAxis object.
        """
        return _JL.ScrewAxis(self.axis, self.n, self.m)

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
        return ScrewAxis(screw_axis_jl.axis, screw_axis_jl.n, screw_axis_jl.m)

    def __repr__(self):
        """
        Return string representation of ScrewAxis.
        """
        return f"ScrewAxis(axis='{self.axis}',n={self.n},m={self.m})"
