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
Centering class
"""
# --- Imports

# Standard library
from enum import auto, StrEnum

# Local packages/modules
from .. import _JL


# --- Classes


class Centering(StrEnum):
    PRIMITIVE = auto()
    BASE = auto()
    BODY = auto()
    FACE = auto()

    def to_julia(self):
        """
        Convert a Python Centering object to a Julia Centering object.
        """
        if self.value == "primitive":
            return _JL.primitive_centering

        elif self.value == "base":
            return _JL.base_centering

        elif self.value == "body":
            return _JL.body_centering

        elif self.value == "face":
            return _JL.face_centering

    @classmethod
    def from_julia(cls, centering_jl: _JL.Centering):
        """
        Convert a Julia Centering object to a Python Centering object.
        """
        # Check arguments
        if not _JL.isa(centering_jl, _JL.Centering):
            raise ValueError(
                "`centering_jl` must be a Julia `Centering` object. "
                f"(centering_jl={centering_jl})."
            )

        # Convert centering_jl to a Centering object
        if centering_jl == _JL.primitive_centering:
            centering = Centering.PRIMITIVE
        elif centering_jl == _JL.base_centering:
            centering = Centering.BASE
        elif centering_jl == _JL.body_centering:
            centering = Centering.BODY
        elif centering_jl == _JL.face_centering:
            centering = Centering.FACE
        else:
            centering_jl_type = _JL.nameof(_JL.typeof(centering_jl))
            raise ValueError(
                f"Unsupported Centering type. (centering_jl={centering_jl_type})"
            )

        return centering

    @classmethod
    def values(cls):
        """
        Return full list of Centering values.
        """
        return list(map(lambda c: c.value, cls))
