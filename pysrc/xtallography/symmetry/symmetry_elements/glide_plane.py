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
