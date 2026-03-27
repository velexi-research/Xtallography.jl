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
Lattice type and functions
"""
# --- Imports

# Standard library
from dataclasses import dataclass
from typing import Union

# Local packages/modules
from .centering import Centering
from .lattice_system import LatticeSystem


# --- Classes


@dataclass(frozen=True)
class Lattice:
    """
    Class representing a lattice

    Fields
    ------
    * `lattice_system` (LatticeSystem): lattice system

    * `centering` (Centering): centering of lattice
    """

    lattice_system: LatticeSystem
    centering: Centering = Centering.PRIMITIVE

    # --- Initializer

    def __init__(
        self,
        lattice_system: Union[LatticeSystem, str],
        centering: Union[Centering, str] = Centering.PRIMITIVE,
    ):

        # --- Check arguments

        # ------ `lattice_system`

        if not isinstance(lattice_system, Union[LatticeSystem, str]):
            raise ValueError(
                "`lattice_system` must be a LatticeSystem object or a string. "
                f"(lattice_system={lattice_system})"
            )

        # Convert `lattice_system` to LatticeSystem object
        if isinstance(lattice_system, str):
            lattice_system = lattice_system.lower()
            try:
                lattice_system = LatticeSystem(lattice_system)
            except ValueError as error:
                if "is not a valid LatticeSystem" in str(error):
                    raise ValueError(
                        "`lattice_system` must be one of "
                        f"{LatticeSystem.values()}. (lattice_system='{lattice_system}')."
                    )

        # ------ `centering`

        if not isinstance(centering, Union[Centering, str]):
            raise ValueError(
                "`centering` must be a Centering object or a string. "
                f"(centering={centering})"
            )

        # Convert `centering` to Centering object
        if isinstance(lattice_system, str):
            centering = centering.lower()
            try:
                centering = Centering(centering)
            except ValueError as error:
                if "is not a valid Centering" in str(error):
                    raise ValueError(
                        "`centering` must be one of "
                        f"{Centering.values()}. (centering='{centering}')."
                    )

        # --- Initialize field values

        # for frozen DataClasses, field values cannot be set directly
        object.__setattr__(self, "lattice_system", lattice_system)
        object.__setattr__(self, "centering", centering)

    def is_bravais_lattice(self) -> bool:
        """
        Return `True` if Lattice object is a valid Bravais lattice; return `False`
        otherwise.
        """
        return self in BRAVAIS_LATTICES


# --- Constants


BRAVAIS_LATTICES = (
    Lattice(LatticeSystem.TRICLINIC, Centering.PRIMITIVE),
    Lattice(LatticeSystem.MONOCLINIC, Centering.PRIMITIVE),
    Lattice(LatticeSystem.MONOCLINIC, Centering.BASE),
    Lattice(LatticeSystem.MONOCLINIC, Centering.BODY),
    Lattice(LatticeSystem.ORTHORHOMBIC, Centering.PRIMITIVE),
    Lattice(LatticeSystem.ORTHORHOMBIC, Centering.BASE),
    Lattice(LatticeSystem.ORTHORHOMBIC, Centering.BODY),
    Lattice(LatticeSystem.ORTHORHOMBIC, Centering.FACE),
    Lattice(LatticeSystem.TETRAGONAL, Centering.PRIMITIVE),
    Lattice(LatticeSystem.TETRAGONAL, Centering.BODY),
    Lattice(LatticeSystem.RHOMBOHEDRAL, Centering.PRIMITIVE),
    Lattice(LatticeSystem.HEXAGONAL, Centering.PRIMITIVE),
    Lattice(LatticeSystem.CUBIC, Centering.PRIMITIVE),
    Lattice(LatticeSystem.CUBIC, Centering.BODY),
    Lattice(LatticeSystem.CUBIC, Centering.FACE),
)
