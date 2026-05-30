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
Lattice system class
"""
# --- Imports

# Standard library
from enum import auto, StrEnum

# Local packages/modules
from .. import _JL


# --- Classes


class LatticeSystem(StrEnum):
    TRICLINIC = auto()
    MONOCLINIC = auto()
    ORTHORHOMBIC = auto()
    TETRAGONAL = auto()
    RHOMBOHEDRAL = auto()
    HEXAGONAL = auto()
    CUBIC = auto()

    def to_julia(self):
        """
        Convert a Python LatticeSystem object to a Julia LatticeSystem object.
        """
        if self.value == "triclinic":
            return _JL.triclinic

        elif self.value == "monoclinic":
            return _JL.monoclinic

        elif self.value == "orthorhombic":
            return _JL.orthorhombic

        elif self.value == "tetragonal":
            return _JL.tetragonal

        elif self.value == "rhombohedral":
            return _JL.rhombohedral

        elif self.value == "hexagonal":
            return _JL.hexagonal

        elif self.value == "cubic":
            return _JL.cubic

    @classmethod
    def from_julia(cls, lattice_system_jl: _JL.LatticeSystem):
        """
        Convert a Julia LatticeSystem object to a Python LatticeSystem object.
        """
        # Check arguments
        if not _JL.isa(lattice_system_jl, _JL.LatticeSystem):
            raise ValueError(
                "`lattice_system_jl` must be a Julia `LatticeSystem` object. "
                f"(lattice_system_jl={lattice_system_jl})."
            )

        # Convert lattice_system_jl to a LatticeSystem object
        if lattice_system_jl == _JL.triclinic:
            lattice_system = LatticeSystem.TRICLINIC
        elif lattice_system_jl == _JL.monoclinic:
            lattice_system = LatticeSystem.MONOCLINIC
        elif lattice_system_jl == _JL.orthorhombic:
            lattice_system = LatticeSystem.ORTHORHOMBIC
        elif lattice_system_jl == _JL.tetragonal:
            lattice_system = LatticeSystem.TETRAGONAL
        elif lattice_system_jl == _JL.rhombohedral:
            lattice_system = LatticeSystem.RHOMBOHEDRAL
        elif lattice_system_jl == _JL.hexagonal:
            lattice_system = LatticeSystem.HEXAGONAL
        elif lattice_system_jl == _JL.cubic:
            lattice_system = LatticeSystem.CUBIC
        else:
            lattice_system_jl_type = _JL.nameof(_JL.typeof(lattice_system_jl))
            raise ValueError(
                "Unsupported LatticeSystem type. "
                f"(lattice_system_jl={lattice_system_jl_type})"
            )

        return lattice_system

    @classmethod
    def values(cls):
        """
        Return full list of LatticeSystem values.
        """
        return list(map(lambda c: c.value, cls))
