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
