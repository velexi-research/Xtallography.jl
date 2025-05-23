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
The `xtallography.lattices` module provides support for integration with `Xtallography.jl`
lattice/unit cell types and methods.
"""
# --- Imports

# Local packages/modules
from .core import (
    LatticeSystem,
    Centering,
    Lattice,
    UnitCell,
    BRAVAIS_LATTICES,
    create_lattice,
    standardize_lattice,
    is_bravais_lattice,
)
from .triclinic import TriclinicUnitCell
from .monoclinic import MonoclinicUnitCell
from .orthorhombic import OrthorhombicUnitCell
from .tetragonal import TetragonalUnitCell
from .rhombohedral import RhombohedralUnitCell
from .hexagonal import HexagonalUnitCell
from .cubic import CubicUnitCell
from .unit_cell_dynamic import unit_cell_from_julia


# --- Dynamically update UnitCell class


setattr(UnitCell, "from_julia", classmethod(unit_cell_from_julia))


# --- Exports

__all__ = [
    "LatticeSystem",
    "Centering",
    "Lattice",
    "UnitCell",
    "BRAVAIS_LATTICES",
    "create_lattice",
    "standardize_lattice",
    "is_bravais_lattice",
    "TriclinicUnitCell",
    "MonoclinicUnitCell",
    "OrthorhombicUnitCell",
    "TetragonalUnitCell",
    "RhombohedralUnitCell",
    "HexagonalUnitCell",
    "CubicUnitCell",
]
