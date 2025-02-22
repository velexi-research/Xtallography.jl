# Copyright (c) 2024 Velexi Corporation
#
# Use of this software is governed by the Business Source License included
# in the LICENSE file and at www.mariadb.com/bsl11.
#
# Change Date: Four years from the date the Licensed Work is published
#
# Change License: Mozilla Public License 2.0
"""
Core types and functions that support lattice and unit cell computations
"""

# --- Imports

# Standard library
from abc import abstractmethod, ABC
from collections import namedtuple
from enum import auto, StrEnum
import math

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
        Convert Centering object to Julia struct.
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
        Convert UnitCell object to a Julia struct.
        """


# Functions
# DEPRECATED
def is_valid_unit_cell(
    lattice: Lattice, unit_cell: dict, allow_edge_cases: bool = False
) -> bool:
    """
    Return `True` if `unit_cell` defines lattice constants for `lattice`; otherwise,
    raise an exception with relevant message.

    Parameters
    ----------
    unit_cell: lattice parameters for unit cell
    """
    # --- Check arguments

    # Check `lattice` is a `Lattice` type
    if not isinstance(lattice, Lattice):
        raise ValueError(
            f"`lattice` must be a `Lattice` object. (lattice={lattice})"
        )

    # Check `lattice` is a valid Bravais lattice
    if not is_bravais_lattice(lattice):
        raise ValueError(f"`lattice` must be a Bravais lattice. (lattice={lattice})")

    # Check `unit_cell`
    if not isinstance(unit_cell, dict):
        raise ValueError(f"`unit_cell` must be a dict. (unit_cell={unit_cell})")

    # --- Check that `unit_cell` contains parameters for `lattice`

    # if set(unit_cell.keys()) != set(LATTICE_CONSTANTS[lattice.lattice_system]):
    #    raise ValueError(
    #        "`unit_cell` is not compatible with `lattice`. "
    #        f"(lattice={lattice},unit_cell.keys={tuple(unit_cell.keys())})"
    #    )

    # --- Validate values contained in `unit_cell`

    for key, value in unit_cell.items():
        if not isinstance(value, (int, float)):
            raise ValueError(
                f"`unit_cell['{key}']` must be a number. (unit_cell['{key}']={value})"
            )

        if not allow_edge_cases:
            if value <= 0:
                raise ValueError(
                    f"`unit_cell['{key}']` must be positive. (unit_cell['{key}']={value})"
                )
            if math.isinf(value):
                raise ValueError(
                    f"`unit_cell['{key}']` must be finite. (unit_cell['{key}']={value})"
                )
        else:
            if value < 0:
                raise ValueError(
                    f"`unit_cell['{key}']` must be nonnegative. "
                    f"(unit_cell['{key}']={value})"
                )

    return True
