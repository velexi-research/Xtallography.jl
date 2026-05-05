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
The triclinic module defines classes and methods specific to triclinic lattices.
"""

# --- Imports

# Standard library
import copy
import math
from typing import Optional, Union

# Local packages/modules
from xtallography.symmetry import Centering, LatticeSystem

from .. import _JL
from .lattice_constants import TriclinicLatticeConstants
from .unit_cell import UnitCell
from .unit_cell_symmetry import UnitCellSymmetry


# --- Classes


class TriclinicUnitCell(UnitCell):
    """
    Class representing a triclinic unit cell
    """

    # --- Initializer

    def __init__(
        self,
        a: float,
        b: float,
        c: float,
        alpha: float,
        beta: float,
        gamma: float,
        centering: Optional[Centering] = Centering.PRIMITIVE,
        symmetry_elements: Optional[set] = copy.deepcopy(set()),
        check_angle_constraints: Optional[bool] = True,
    ):
        """
        Initialize TriclinicUnitCell object.

        Parameters
        ----------
        `a`, `b`, `c`, `alpha`, `beta`, `gamma`: lattice constants

        `centering`: centering of unit cell

        `symmetry_elements`: symmetry elements that unit cell is invariant under
        """
        # --- Check arguments

        if a <= 0:
            raise ValueError(f"`a` must be positive. (a={a})")

        if b <= 0:
            raise ValueError(f"`b` must be positive. (b={b})")

        if c <= 0:
            raise ValueError(f"`c` must be positive. (c={c})")

        if alpha <= 0 or alpha >= 2 * math.pi:
            raise ValueError(
                f"`alpha` must lie in the interval (0, 2 pi). (alpha={alpha})"
            )

        if beta <= 0 or beta >= 2 * math.pi:
            raise ValueError(
                f"`beta` must lie in the interval (0, 2 pi). (beta={beta})"
            )

        if gamma <= 0 or gamma >= 2 * math.pi:
            raise ValueError(
                f"`gamma` must lie in the interval (0, 2 pi). (gamma={gamma})"
            )

        if check_angle_constraints and not self.satisfies_angle_constraints(
            alpha, beta, gamma
        ):
            raise ValueError(
                "`alpha`, `beta`, and `gamma` do not satisfy the angle constraints for "
                f"a triclinic unit cell. (alpha={alpha},beta={beta},gamma={gamma})"
            )

        # --- Initialize parent class

        lattice_constants = TriclinicLatticeConstants(a, b, c, alpha, beta, gamma)
        symmetry = UnitCellSymmetry(centering, symmetry_elements)
        super().__init__(LatticeSystem.TRICLINIC, lattice_constants, symmetry)

    # --- Properties

    @property
    def a(self):
        """
        Return `a`.
        """
        return self.lattice_constants.a

    @property
    def b(self):
        """
        Return `b`.
        """
        return self.lattice_constants.b

    @property
    def c(self):
        """
        Return `c`.
        """
        return self.lattice_constants.c

    @property
    def alpha(self):
        """
        Return `alpha`.
        """
        return self.lattice_constants.alpha

    @property
    def beta(self):
        """
        Return `beta`.
        """
        return self.lattice_constants.beta

    @property
    def gamma(self):
        """
        Return `gamma`.
        """
        return self.lattice_constants.gamma

    # --- Methods

    def to_julia(self):
        """
        Convert TriclinicUnitCell object to a Julia UnitCell object.
        """
        return _JL.TriclinicUnitCell(
            self.a,
            self.b,
            self.c,
            self.alpha,
            self.beta,
            self.gamma,
            centering=self.centering.to_julia(),
            symmetry_elements=_JL.Vector(
                [element.to_julia() for element in self.symmetry_elements]
            ),
        )

    @classmethod
    def from_julia(cls, unit_cell_jl: _JL.UnitCell):
        """
        Convert a Julia UnitCell object to a TriclinicUnitCell object.
        """
        # --- Check arguments

        if not _JL.isa(unit_cell_jl, _JL.TriclinicUnitCell):
            raise ValueError(
                "`unit_cell_jl` must be a Julia `TriclinicUnitCell` object. "
                f"(unit_cell_jl={unit_cell_jl})."
            )

        # --- Convert unit_cell_jl to a TriclinicUnitCell object

        unit_cell_symmetry = UnitCellSymmetry.from_julia(unit_cell_jl.symmetry)
        unit_cell = TriclinicUnitCell(
            unit_cell_jl.lattice_constants.a,
            unit_cell_jl.lattice_constants.b,
            unit_cell_jl.lattice_constants.c,
            unit_cell_jl.lattice_constants.α,
            unit_cell_jl.lattice_constants.β,
            unit_cell_jl.lattice_constants.γ,
            centering=unit_cell_symmetry.centering,
            symmetry_elements=unit_cell_symmetry.symmetry_elements,
        )

        return unit_cell

    def __repr__(self):
        """
        Return string representation of TriclinicUnitCell.
        """
        return (
            f"TriclinicUnitCell(a={self.a},b={self.b},c={self.c},"
            f"alpha={self.alpha},beta={self.beta},gamma={self.gamma},"
            f"centering={self.centering},"
            "symmetry_elements=["
            f"{','.join(sorted([str(item) for item in self.symmetry_elements]))}"
            "])"
        )

    @classmethod
    def satisfies_angle_constraints(
        cls, alpha: Union[float, int], beta: Union[float, int], gamma: Union[float, int]
    ) -> bool:
        """
        Determine whether `alpha`, `beta`, and `beta` satisfy the angle constraints for
        a triclinic unit cell:

        * 0 <  alpha + beta + gamma < 2 * pi
        * 0 <  alpha + beta - gamma < 2 * pi
        * 0 <  alpha - beta + gamma < 2 * pi
        * 0 < -alpha + beta + gamma < 2 * pi

        Return values
        =============
        - `True` if (`alpha`, `beta`, `gamma`) form a valid triple of angles for a
          triclinic unit cell; `False` otherwise

        >>> TriclinicUnitCell.satisfies_angle_constraints(pi/4, pi/5, pi/6)
        True
        >>> TriclinicUnitCell.satisfies_angle_constraints(3*pi/4, 4*pi/5, 5*pi/6)
        False
        """
        return (
            (0 < alpha + beta + gamma < 2 * math.pi)
            and (0 < alpha + beta - gamma < 2 * math.pi)
            and (0 < alpha - beta + gamma < 2 * math.pi)
            and (0 < -alpha + beta + gamma < 2 * math.pi)
        )
