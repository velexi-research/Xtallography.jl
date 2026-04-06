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
Lattice constant classes
"""

# --- Imports

# Standard library
from abc import ABC
from dataclasses import dataclass
from typing import Union


# --- Classes


class LatticeConstants(ABC):
    """
    Abstract base class for lattice constant classes
    """


@dataclass(frozen=True)
class TriclinicLatticeConstants(LatticeConstants):
    """
    Class representing the lattice constants for a triclinic unit cell

    TODO: add fields
    """

    a: Union[float, int]
    b: Union[float, int]
    c: Union[float, int]
    alpha: Union[float, int]
    beta: Union[float, int]
    gamma: Union[float, int]


@dataclass(frozen=True)
class MonoclinicLatticeConstants(LatticeConstants):
    """
    Class representing the lattice constants for a monoclinic unit cell

    TODO: add fields
    """

    a: Union[float, int]
    b: Union[float, int]
    c: Union[float, int]
    beta: Union[float, int]


@dataclass(frozen=True)
class OrthorhombicLatticeConstants(LatticeConstants):
    """
    Class representing the lattice constants for a orthorhombic unit cell

    TODO: add fields
    """

    a: Union[float, int]
    b: Union[float, int]
    c: Union[float, int]


@dataclass(frozen=True)
class HexagonalLatticeConstants(LatticeConstants):
    """
    Class representing the lattice constants for a hexagonal unit cell

    TODO: add fields
    """

    a: Union[float, int]
    c: Union[float, int]


@dataclass(frozen=True)
class RhombohedralLatticeConstants(LatticeConstants):
    """
    Class representing the lattice constants for a rhombohedral unit cell

    TODO: add fields
    """

    a: Union[float, int]
    alpha: Union[float, int]


@dataclass(frozen=True)
class TetragonalLatticeConstants(LatticeConstants):
    """
    Class representing the lattice constants for a tetragonal unit cell

    TODO: add fields
    """

    a: Union[float, int]
    c: Union[float, int]


@dataclass(frozen=True)
class CubicLatticeConstants(LatticeConstants):
    """
    Class representing the lattice constants for a cubic unit cell

    TODO: add fields
    """

    a: Union[float, int]
