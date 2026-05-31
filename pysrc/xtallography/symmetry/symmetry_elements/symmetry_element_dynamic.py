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
Dynamic additions to SymmetryElement class
"""

# --- Imports

# Local packages/modules
from xtallography import _JL
from .rotation_axis import RotationAxis
from .glide_plane import GlidePlane
from .screw_axis import ScrewAxis

# --- Methods


def symmetry_element_from_julia(cls, symmetry_element_jl: _JL.SymmetryElement):
    """
    Convert a Julia SymmetryElement object to a Python SymmetryElement object.

    Parameters
    ----------
    symmetry_element_jl: Julia SymmetryElement object

    Return value
    ------------
    Python SymmetryElement object
    """
    # --- Check arguments

    if not _JL.isa(symmetry_element_jl, _JL.SymmetryElement):
        raise ValueError(
            "`symmetry_element_jl` must be a Julia `SymmetryElement` object. "
            f"(symmetry_element_jl={symmetry_element_jl})."
        )

    # --- Call from_julia() from the appropriate subclass of SymmetryElement

    # TODO. Add other symmetry elements

    if _JL.isa(symmetry_element_jl, _JL.RotationAxis):
        symmetry_element = RotationAxis.from_julia(symmetry_element_jl)

    elif _JL.isa(symmetry_element_jl, _JL.GlidePlane):
        symmetry_element = GlidePlane.from_julia(symmetry_element_jl)

    elif _JL.isa(symmetry_element_jl, _JL.ScrewAxis):
        symmetry_element = ScrewAxis.from_julia(symmetry_element_jl)

    return symmetry_element
