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
Dynamic additions to UnitCell class
"""

# --- Imports

# External packages
from xtallography import _JL

# Local packages/modules
from .core import LatticeSystem
from .triclinic import TriclinicUnitCell
from .monoclinic import MonoclinicUnitCell
from .orthorhombic import OrthorhombicUnitCell
from .tetragonal import TetragonalUnitCell
from .rhombohedral import RhombohedralUnitCell
from .hexagonal import HexagonalUnitCell
from .cubic import CubicUnitCell


# --- Methods


def unit_cell_from_julia(cls, unit_cell_jl: _JL.UnitCell):
    """
    Convert a Julia UnitCell object to a Python UnitCell object.

    Parameters
    ----------
    unit_cell_jl: Julia UnitCell object

    Return value
    ------------
    unit_cell: Python UnitCell object
    """
    # --- Check arguments

    if not _JL.isa(unit_cell_jl, _JL.UnitCell):
        raise ValueError(
            "`unit_cell_jl` must be a Julia `UnitCell` object. "
            f"(unit_cell_jl={unit_cell_jl})."
        )

    # --- Call from_julia() from the appropriate subclass of UnitCell

    lattice_system = LatticeSystem.from_julia(_JL.lattice_system(unit_cell_jl))

    if lattice_system == LatticeSystem.TRICLINIC:
        unit_cell = TriclinicUnitCell.from_julia(unit_cell_jl)

    elif lattice_system == LatticeSystem.MONOCLINIC:
        unit_cell = MonoclinicUnitCell.from_julia(unit_cell_jl)

    elif lattice_system == LatticeSystem.ORTHORHOMBIC:
        unit_cell = OrthorhombicUnitCell.from_julia(unit_cell_jl)

    elif lattice_system == LatticeSystem.TETRAGONAL:
        unit_cell = TetragonalUnitCell.from_julia(unit_cell_jl)

    elif lattice_system == LatticeSystem.RHOMBOHEDRAL:
        unit_cell = RhombohedralUnitCell.from_julia(unit_cell_jl)

    elif lattice_system == LatticeSystem.HEXAGONAL:
        unit_cell = HexagonalUnitCell.from_julia(unit_cell_jl)

    elif lattice_system == LatticeSystem.CUBIC:
        unit_cell = CubicUnitCell.from_julia(unit_cell_jl)

    return unit_cell
