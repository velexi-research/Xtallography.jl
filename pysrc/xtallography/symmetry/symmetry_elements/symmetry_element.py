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
Abstract base class for symmetry element classes
"""

# --- Imports

# Standard library
from abc import abstractmethod, ABC

# Local packages/modules
from ... import _JL


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
        # Dynamically defined in symmetry_elements/__init__.py

    @abstractmethod
    def __repr__(self):
        """
        Return string representation of SymmetryElement.
        """
