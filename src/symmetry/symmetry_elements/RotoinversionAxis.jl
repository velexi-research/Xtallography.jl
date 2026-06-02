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
RotoinversionAxis type and functions
"""
# --- Exports

# Types
export RotoinversionAxis

# --- Types

"""
    RotoinversionAxis

Type representing a rotoinversion axis

Supertype: [`SymmetryElement`](@ref)
"""
struct RotoinversionAxis <: SymmetryElement
    # --- Fields

    # rotoinversion order
    n::Int

    # direction of rotoinversion axis
    direction::Tuple{Rational,Rational,Rational}

    # location of inversion center in unit cell fractional coordinates
    center::Tuple{Rational,Rational,Rational}

    # Constructor
    function RotoinversionAxis(
        n::Int, direction::Tuple{<:Real,<:Real,<:Real}, center::Tuple{<:Real,<:Real,<:Real}
    )
        # --- Check arguments

        # n be positive
        if n ≤ 0
            throw(ArgumentError("`n` must be positive (n=$n)"))
        end

        # Return new RotoinversionAxis
        return new(n, direction, center)
    end
end

# Outer constructor
"""
    RotoinversionAxis(
        n::Int,
        direction::Tuple{<:Real,<:Real,<:Real};
        center::Tuple{<:Real,<:Real,<:Real}=(0, 0, 0)
    )

Construct a RotoinversionAxis object with order `n` about the axis specified by `direction`.

Keyword Arguments
=================
- `center`: location of inversion center in unit cell fractional coordinates
"""
function RotoinversionAxis(
    n::Int,
    direction::Tuple{<:Real,<:Real,<:Real};
    center::Tuple{<:Real,<:Real,<:Real}=(0, 0, 0),
)
    return RotoinversionAxis(n, direction, center)
end

# --- Functions/Methods

import Base.:(==)
import LinearAlgebra: dot

function Base.:(==)(x::RotoinversionAxis, y::RotoinversionAxis)
    # Check that orders are equal
    if x.n != y.n
        return false
    end

    # Check inversion centers are the same
    if x.center != y.center
        return false
    end

    # Check that directions are the same
    return dot(x.direction, y.direction)^2 ==
           dot(x.direction, x.direction) * dot(y.direction, y.direction)
end
