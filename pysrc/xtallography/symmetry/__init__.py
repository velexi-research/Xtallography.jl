#   Copyright 2026 Velexi Corporation
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
The `xtallography.symmetry` package defines classes and methods related to crystal symmetry.
"""
# Package-level types

from .lattice import Lattice  # noqa
from .lattice import BRAVAIS_LATTICES  # noqa

from .centering import Centering  # noqa

from .lattice_system import LatticeSystem  # noqa

from .symmetry_elements import SymmetryElement  # noqa
from .symmetry_elements import RotationAxis  # noqa
from .symmetry_elements import GlidePlane  # noqa
from .symmetry_elements import ScrewAxis  # noqa


# --- Auto-doc

BRAVAIS_LATTICES = BRAVAIS_LATTICES
"""
List of Bravais lattices
"""

# --- Exports

__all__ = [
    "LatticeSystem",
    "Centering",
    "Lattice",
    "BRAVAIS_LATTICES",
    "SymmetryElement",
    "RotationAxis",
    "GlidePlane",
    "ScrewAxis",
]
