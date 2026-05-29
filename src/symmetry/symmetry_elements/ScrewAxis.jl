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
ScrewAxis type and functions
"""
# --- Exports

# Types
export ScrewAxis

# Constants
export a_2_1, a_4_1, a_4_2, a_4_3
export b_2_1, b_4_1, b_4_2, b_4_3
export c_2_1, c_4_1, c_4_2, c_4_3
export c_3_1, c_3_2
export c_6_1, c_6_2, c_6_3, c_6_4, c_6_5

# --- Types

"""
    ScrewAxis

Type representing a screw axis

Supertype: [`SymmetryElement`](@ref)
"""
struct ScrewAxis <: SymmetryElement
    # --- Fields

    # rotation order
    n::Int

    # number of translation steps of size 1/n following rotation by 2π/n
    m::Int

    # direction of rotation axis in unit cell fractional coordinates
    direction::Tuple{Rational,Rational,Rational}

    # a point on the rotation axis in unit cell fractional coordinates
    location::Tuple{Rational,Rational,Rational}

    # --- Constructor

    function ScrewAxis(
        n::Int,
        m::Int,
        direction::Tuple{<:Real,<:Real,<:Real},
        location::Tuple{<:Real,<:Real,<:Real},
    )

        # Enforce constraints
        if n ≤ 0
            throw(ArgumentError("`n` must be positive (n=$n)"))
        end

        if m ≤ 0
            throw(ArgumentError("`m` must be positive (m=$m)"))
        end

        if m ≥ n
            throw(ArgumentError("`m` must be less than `n` (n=$n,m=$m)"))
        end

        # Return new ScrewAxis
        return new(n, m, direction, location)
    end
end

# Outer constructor
"""
    ScrewAxis(
        n::Int,
        m::Int,
        direction::Tuple{<:Real,<:Real,<:Real};
        location::Tuple{<:Real,<:Real,<:Real}=(0, 0, 0)
    )

Construct a ScrewAxis object that

* has an axis by `direction` and `location`

and

* couples a rotation by `2π/n` and a translation along the rotation axis by `m/n` (in unit
  cell fractional coordinates).

Keyword Arguments
=================
- `location`: a point on the rotation axis in unit cell fractional coordinates
"""
function ScrewAxis(
    n::Int,
    m::Int,
    direction::Tuple{<:Real,<:Real,<:Real};
    location::Tuple{<:Real,<:Real,<:Real}=(0, 0, 0),
)
    return ScrewAxis(n, m, direction, location)
end

# --- Functions/Methods

import Base.:(==)
import LinearAlgebra: dot

function Base.:(==)(x::ScrewAxis, y::ScrewAxis)
    # Check that rotation orders and translation steps are equal
    if x.n != y.n
        return false
    end

    if x.m != y.m
        return false
    end

    # Check that directions are the same
    if dot(x.direction, y.direction)^2 !=
        dot(x.direction, x.direction) * dot(y.direction, y.direction)
        return false
    end

    # Check that line through both locations is in the same direction as both lines
    delta = (x.location[i] - y.location[i] for i in 1:3)
    return dot(delta, x.direction)^2 == dot(delta, delta) * dot(x.direction, x.direction)
end

# --- Constants

const a_2_1 = ScrewAxis(2, 1, (1, 0, 0))
const a_4_1 = ScrewAxis(4, 1, (1, 0, 0))
const a_4_2 = ScrewAxis(4, 2, (1, 0, 0))
const a_4_3 = ScrewAxis(4, 3, (1, 0, 0))

const b_2_1 = ScrewAxis(2, 1, (0, 1, 0))
const b_4_1 = ScrewAxis(4, 1, (0, 1, 0))
const b_4_2 = ScrewAxis(4, 2, (0, 1, 0))
const b_4_3 = ScrewAxis(4, 3, (0, 1, 0))

const c_2_1 = ScrewAxis(2, 1, (0, 0, 1))
const c_3_1 = ScrewAxis(3, 1, (0, 0, 1))
const c_3_2 = ScrewAxis(3, 2, (0, 0, 1))
const c_4_1 = ScrewAxis(4, 1, (0, 0, 1))
const c_4_2 = ScrewAxis(4, 2, (0, 0, 1))
const c_4_3 = ScrewAxis(4, 3, (0, 0, 1))
const c_6_1 = ScrewAxis(6, 1, (0, 0, 1))
const c_6_2 = ScrewAxis(6, 2, (0, 0, 1))
const c_6_3 = ScrewAxis(6, 3, (0, 0, 1))
const c_6_4 = ScrewAxis(6, 4, (0, 0, 1))
const c_6_5 = ScrewAxis(6, 5, (0, 0, 1))
