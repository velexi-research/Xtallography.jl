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
Core types and functions that support lattice and unit cell computations
"""

# --- Imports

# Standard library
from abc import abstractmethod, ABC
from collections import namedtuple
from enum import auto, StrEnum
import math
import sys
from typing import Optional

# Local packages/modules
from .. import _JL


# --- Types


# ----- Enumerations


class LatticeSystem(StrEnum):
    TRICLINIC = auto()
    MONOCLINIC = auto()
    ORTHORHOMBIC = auto()
    TETRAGONAL = auto()
    RHOMBOHEDRAL = auto()
    HEXAGONAL = auto()
    CUBIC = auto()

    @classmethod
    def values(cls):
        """
        Return full list of LatticeSystem values.
        """
        return list(map(lambda c: c.value, cls))


class Centering(StrEnum):
    PRIMITIVE = auto()
    BASE_CENTERED = auto()
    BODY_CENTERED = auto()
    FACE_CENTERED = auto()

    def to_julia(self):
        """
        Convert a Python Centering object to a Julia Centering object.
        """
        if self.value == "primitive":
            return _JL.primitive

        elif self.value == "base_centered":
            return _JL.base_centered

        elif self.value == "body_centered":
            return _JL.body_centered

        elif self.value == "face_centered":
            return _JL.face_centered

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
        if centering_jl == _JL.primitive:
            centering = Centering.PRIMITIVE
        elif centering_jl == _JL.base_centered:
            centering = Centering.BASE_CENTERED
        elif centering_jl == _JL.body_centered:
            centering = Centering.BODY_CENTERED
        elif centering_jl == _JL.face_centered:
            centering = Centering.FACE_CENTERED
        else:
            raise ValueError(
                f"Unsupported Centering type. (centering_jl={centering_jl})"
            )

        return centering

    @classmethod
    def values(cls):
        """
        Return full list of Centering values.
        """
        return list(map(lambda c: c.value, cls))


# ----- Lattice


Lattice = namedtuple("Lattice", "lattice_system centering")


# Constants
BRAVAIS_LATTICES = (
    Lattice(LatticeSystem.TRICLINIC, Centering.PRIMITIVE),
    Lattice(LatticeSystem.MONOCLINIC, Centering.PRIMITIVE),
    Lattice(LatticeSystem.MONOCLINIC, Centering.BASE_CENTERED),
    Lattice(LatticeSystem.MONOCLINIC, Centering.BODY_CENTERED),
    Lattice(LatticeSystem.ORTHORHOMBIC, Centering.PRIMITIVE),
    Lattice(LatticeSystem.ORTHORHOMBIC, Centering.BASE_CENTERED),
    Lattice(LatticeSystem.ORTHORHOMBIC, Centering.BODY_CENTERED),
    Lattice(LatticeSystem.ORTHORHOMBIC, Centering.FACE_CENTERED),
    Lattice(LatticeSystem.TETRAGONAL, Centering.PRIMITIVE),
    Lattice(LatticeSystem.TETRAGONAL, Centering.BODY_CENTERED),
    Lattice(LatticeSystem.RHOMBOHEDRAL, Centering.PRIMITIVE),
    Lattice(LatticeSystem.HEXAGONAL, Centering.PRIMITIVE),
    Lattice(LatticeSystem.CUBIC, Centering.PRIMITIVE),
    Lattice(LatticeSystem.CUBIC, Centering.BODY_CENTERED),
    Lattice(LatticeSystem.CUBIC, Centering.FACE_CENTERED),
)


# Functions
def create_lattice(
    lattice_system: str | None = None,
    centering: str | None = None,
) -> Lattice:
    """
    Create a Lattice object.

    Parameters
    ----------
    lattice_system: lattice system

    centering: lattice centering

    Return value
    ------------
    lattice: Lattice object with the specified `lattice_system` and `centering`
    """
    # --- Check arguments

    # ------ `lattice_system`

    if not isinstance(lattice_system, str):
        raise ValueError(
            f"`lattice_system` must be a string. (lattice_system={lattice_system})"
        )

    # Enforce that `lattice_system` is lowercase
    lattice_system = lattice_system.lower()

    # ------ `centering`

    if not isinstance(centering, str):
        raise ValueError(f"`centering` must be a string. (centering={centering})")

    # Enforce that `centering` is lowercase
    centering = centering.lower()

    # --- Construct Lattice

    try:
        lattice = Lattice(LatticeSystem(lattice_system), Centering(centering))

    except ValueError as error:
        if "is not a valid LatticeSystem" in str(error):
            raise ValueError(
                "`lattice_system` must be one of "
                f"{LatticeSystem.values()}. (lattice_system='{lattice_system}')."
            )

        if "is not a valid Centering" in str(error):
            raise ValueError(
                "`centering` must be one of "
                f"{Centering.values()}. (centering='{centering}')."
            )

    return lattice


def standardize_lattice(lattice: Lattice) -> Lattice:
    """
    Standardize and validate field values for `lattice`.

    Parameters
    ----------
    lattice: Lattice to standardize

    Return value
    ------------
    lattice: standardized Lattice
    """
    # --- Validate types

    for field, value in zip(lattice._fields, lattice):
        if not isinstance(value, str):
            raise ValueError(
                f"`lattice.{field}` must be a string. (lattice.{field}={value})"
            )

    # --- Standardize field values

    standardized_field_values = {}
    for field, value in zip(lattice._fields, lattice):
        standardized_field_values[field] = value.lower()

    # --- Construct Lattice with standardized field values

    lattice = create_lattice(
        standardized_field_values["lattice_system"],
        standardized_field_values["centering"],
    )

    return lattice


def is_bravais_lattice(lattice: Lattice) -> bool:
    """
    Return `True` if `lattice` is a valid Bravais lattice; return `False` otherwise.

    Parameters
    ----------
    lattice: Lattice to check
    """
    return lattice in BRAVAIS_LATTICES


# ------ UnitCell


class UnitCell(ABC):
    """
    Abstract base class for unit cell classes.
    """

    # --- Initializer

    def __init__(
        self, lattice_system: LatticeSystem, centering: Centering = Centering.PRIMITIVE
    ):
        """
        Initialize core of UnitCell object.

        Parameters
        ----------
        `lattice_system`: lattice system

        `centering`: centering
        """
        # --- Check arguments

        # Validate that (lattice_system, centering) is a valid Bravais lattice
        if not is_bravais_lattice(Lattice(lattice_system, centering)):
            raise ValueError(
                f"('{lattice_system}', '{centering}') is not a valid Bravais lattice."
            )

        # Enforce that lattice_system is a LatticeSystem object
        if not isinstance(lattice_system, LatticeSystem):
            lattice_system = LatticeSystem(lattice_system)

        # Enforce that centering is a Centering object
        if not isinstance(centering, Centering):
            centering = Centering(centering)

        # --- Initialize attributes

        self._lattice_system = lattice_system
        self._centering = centering

    # --- Properties

    @property
    def lattice_system(self):
        """
        Return the lattice system of the unit cell.
        """
        return self._lattice_system

    @property
    def centering(self):
        """
        Return the unit cell centering.
        """
        return self._centering

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
        unit_cell: Python UnitCell object
        """

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

        # Compare lattice constants
        for var in vars(self):
            if var == "_lattice_system" or var == "_centering":
                # Compare lattice_system or centering
                if getattr(self, var) != getattr(other, var):
                    return False

            else:
                # Compare lattice constant
                if not math.isclose(
                    getattr(self, var), getattr(other, var), rel_tol=rtol, abs_tol=atol
                ):
                    return False

        return True
