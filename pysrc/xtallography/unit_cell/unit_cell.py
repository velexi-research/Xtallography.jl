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
UnitCell class
"""

# --- Imports

# Standard library
from abc import abstractmethod, ABC
from dataclasses import dataclass, fields
import math
import sys
from typing import Optional, Union

# External packages
import numpy

# Local packages/modules
from .. import _JL

from xtallography.symmetry import LatticeSystem

from .lattice_constants import (
    LatticeConstants,
    TriclinicLatticeConstants,
    MonoclinicLatticeConstants,
    OrthorhombicLatticeConstants,
    HexagonalLatticeConstants,
    RhombohedralLatticeConstants,
    TetragonalLatticeConstants,
    CubicLatticeConstants,
)
from .unit_cell_symmetry import UnitCellSymmetry
from .unit_cell_symmetry import PRIMITIVE_UNIT_CELL_SYMMETRY


# --- Classes


@dataclass(frozen=True)
class UnitCell(ABC):
    """
    Abstract base class for unit cell classes

    Fields
    ------
    * `lattice_system` (LatticeSystem): lattice system of unit cell

    * `lattice_constants` (LatticeConstants): lattice constants of unit cell

    * `symmetry` (UnitCellSymmetry): symmetry of unit cell

    * `unit_cell_jl` (_JL.UnitCell): Julia object representing unit cell

    * `volume` (float): volume of unit cell

    * `surface_area` (float): surface area of unit cell
    """

    lattice_system: LatticeSystem
    lattice_constants: LatticeConstants
    symmetry: UnitCellSymmetry

    # attributes that cannot be set until object is fully constructed
    _unit_cell_jl: _JL.UnitCell
    _basis: list
    _volume: float
    _surface_area: float

    # --- Initializer

    def __init__(
        self,
        lattice_system: Union[LatticeSystem, str],
        lattice_constants: LatticeConstants,
        symmetry: Optional[UnitCellSymmetry] = PRIMITIVE_UNIT_CELL_SYMMETRY,
    ):
        """
        Initialize UnitCell object.

        Parameters
        ----------
        `lattice_system`: lattice system

        `lattice_constants`: lattice constants of unit cell

        `symmetry`: symmetry of unit cell
        """
        # --- Check arguments

        # ------ lattice_system

        # Enforce that lattice_system is a LatticeSystem object
        if isinstance(lattice_system, str):
            lattice_system = LatticeSystem(lattice_system.lower())

        elif not isinstance(lattice_system, LatticeSystem):
            raise ValueError(
                "`lattice_system` must be a `LatticeSystem` or `str`."
                f"(lattice_system={lattice_system})"
            )

        # ------ lattice_constants

        # Check type of lattice_constants
        if lattice_system == LatticeSystem.TRICLINIC and not isinstance(
            lattice_constants, TriclinicLatticeConstants
        ):
            raise ValueError(
                "`lattice_constants` for triclinic unit cell must be a "
                "TriclinicLatticeConstants object. "
                f"(lattice_constants={lattice_constants})"
            )

        elif lattice_system == LatticeSystem.MONOCLINIC and not isinstance(
            lattice_constants, MonoclinicLatticeConstants
        ):
            raise ValueError(
                "`lattice_constants` for monoclinic unit cell must be a "
                "MonoclinicLatticeConstants object. "
                f"(lattice_constants={lattice_constants})"
            )

        elif lattice_system == LatticeSystem.ORTHORHOMBIC and not isinstance(
            lattice_constants, OrthorhombicLatticeConstants
        ):
            raise ValueError(
                "`lattice_constants` for orthorhombic unit cell must be a "
                "OrthorhombicLatticeConstants object. "
                f"(lattice_constants={lattice_constants})"
            )

        elif lattice_system == LatticeSystem.HEXAGONAL and not isinstance(
            lattice_constants, HexagonalLatticeConstants
        ):
            raise ValueError(
                "`lattice_constants` for hexagonal unit cell must be a "
                "HexagonalLatticeConstants object. "
                f"(lattice_constants={lattice_constants})"
            )

        elif lattice_system == LatticeSystem.RHOMBOHEDRAL and not isinstance(
            lattice_constants, RhombohedralLatticeConstants
        ):
            raise ValueError(
                "`lattice_constants` for rhombohedral unit cell must be a "
                "RhombohedralLatticeConstants object. "
                f"(lattice_constants={lattice_constants})"
            )

        elif lattice_system == LatticeSystem.TETRAGONAL and not isinstance(
            lattice_constants, TetragonalLatticeConstants
        ):
            raise ValueError(
                "`lattice_constants` for tetragonal unit cell must be a "
                "TetragonalLatticeConstants object. "
                f"(lattice_constants={lattice_constants})"
            )

        elif lattice_system == LatticeSystem.CUBIC and not isinstance(
            lattice_constants, CubicLatticeConstants
        ):
            raise ValueError(
                "`lattice_constants` for cubic unit cell must be a "
                "CubicLatticeConstants object. "
                f"(lattice_constants={lattice_constants})"
            )

        # ------ symmetry

        # Check that symmetry is a UnitCellSymmetry object
        if not isinstance(symmetry, UnitCellSymmetry):
            raise ValueError(
                f"`symmetry` must be a `UnitCellSymmetry` object. (symmetry={symmetry})"
            )

        # --- Initialize field values

        # for frozen DataClasses, field values cannot be set directly
        object.__setattr__(self, "lattice_system", lattice_system)
        object.__setattr__(self, "lattice_constants", lattice_constants)
        object.__setattr__(self, "symmetry", symmetry)
        object.__setattr__(self, "_unit_cell_jl", None)
        object.__setattr__(self, "_basis", None)
        object.__setattr__(self, "_volume", None)
        object.__setattr__(self, "_surface_area", None)

    # --- Properties

    @property
    def centering(self):
        """
        Return the centering of the unit cell.
        """
        return self.symmetry.centering

    @property
    def symmetry_elements(self):
        """
        Return the symmetry elements of the unit cell.
        """
        return self.symmetry.symmetry_elements

    @property
    def unit_cell_jl(self):
        """
        Julia object representing unit cell.
        """
        if self._unit_cell_jl is None:
            object.__setattr__(self, "_unit_cell_jl", self.to_julia())

        return self._unit_cell_jl

    @property
    def basis(self):
        """
        Return a list of basis vectors for the unit cell.
        """
        if self._basis is None:
            # Compute basis
            basis_jl = _JL.basis(self.unit_cell_jl)

            # Convert basis vectors to NumPy arrays
            basis = []
            for vector_jl in basis_jl:
                basis.append(numpy.array([element for element in vector_jl]))

            object.__setattr__(self, "_basis", basis)

        return self._basis

    @property
    def volume(self):
        """
        Return the volume of the unit cell.
        """
        if self._volume is None:
            object.__setattr__(self, "_volume", _JL.volume(self.unit_cell_jl))

        return self._volume

    @property
    def surface_area(self):
        """
        Return the surface area of the unit cell.
        """
        if self._surface_area is None:
            object.__setattr__(
                self, "_surface_area", _JL.surface_area(self.unit_cell_jl)
            )

        return self._surface_area

    # --- Methods

    @abstractmethod
    def to_julia(self):
        """
        Convert Python UnitCell object to a Julia UnitCell object.
        """

    @classmethod
    def from_julia(cls, unit_cell_jl: _JL.UnitCell):
        """
        Convert a Julia UnitCell object to a Python UnitCell object.

        Parameters
        ----------
        unit_cell_jl: Julia UnitCell object

        Return value
        ------------
        Python UnitCell object
        """
        # Dynamically defined in unit_cell/__init__.py

    @abstractmethod
    def __repr__(self):
        """
        Return string representation of UnitCell.
        """

    def __eq__(self, other):
        """
        Return True if `self` and `other`are identical unit cells; otherwise, return
        False.

        Parameters
        ----------
        other: UnitCell object to compare against
        """
        if not isinstance(other, type(self)):
            return False

        for var in vars(self):
            if getattr(self, var) != getattr(other, var):
                return False

        return True

    def isclose(self, other, rtol: Optional[float] = None, atol: float = 0):
        """
        Return True if `self` and `other` are approximately equal unit cells; otherwise,
        return False.

        Parameters
        ----------
        other: UnitCell object to compare against

        rtol: relative tolerance

        atol: absolute tolerance

        Note
        ----
        * The default values for `rtol` and `atol` are set using the same logic as the
          Julia `isapprox()` method.
        """
        # --- Check arguments

        if atol < 0:
            raise ValueError(f"`atol` must be nonnegative. (atol={atol})")

        if rtol is None:
            if atol > 0:
                rtol = 0
            else:
                rtol = math.sqrt(sys.float_info.epsilon)

        if rtol < 0:
            raise ValueError(f"`rtol` must be nonnegative. (rtol={rtol})")

        # --- Compare UnitCell objects

        # Compare types
        if not isinstance(other, type(self)):
            return False

        # Compare lattice_system
        if self.lattice_system != other.lattice_system:
            return False

        # Compare symmetry
        if self.symmetry != other.symmetry:
            return False

        # Compare lattice constants
        for field in fields(self.lattice_constants):
            if not math.isclose(
                getattr(self.lattice_constants, field.name),
                getattr(other.lattice_constants, field.name),
                rel_tol=rtol,
                abs_tol=atol,
            ):
                return False

        return True

    def compute_reduced_cell(self):
        """
        Compute the primitive reduced cell for the unit cell. The Selling-Delaunay
        reduction algorithm is used to compute the reduced basis.

        Return value
        ------------
        UnitCell object representing the reduced cell
        """
        return UnitCell.from_julia(_JL.reduced_cell(self.unit_cell_jl))

    def compute_conventional_cell(self):
        """
        Return the IUCr conventional cell that is equivalent to the unit cell.

        Return values
        -------------
        UnitCell object representing the IUCr conventional cell for the unit cell
        """
        return UnitCell.from_julia(_JL.conventional_cell(self.unit_cell_jl))
