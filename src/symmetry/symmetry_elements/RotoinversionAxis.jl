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

    # rotation order
    n::Int

    # direction of rotoinversion axis
    direction::Tuple{Rational,Rational,Rational}

    # a point on the rotation axis in unit cell fractional coordinates
    location::Tuple{Rational,Rational,Rational}

    # Constructor
    function RotoinversionAxis(
        n::Int,
        direction::Tuple{<:Real,<:Real,<:Real},
        location::Tuple{<:Real,<:Real,<:Real},
    )
        # Enforce constraints
        if n ≤ 0
            throw(ArgumentError("`n` must be positive (n=$n)"))
        end

        # Return new RotoinversionAxis
        return new(n, direction, location)
    end
end

# Outer constructor
"""
    RotoinversionAxis(
        n::Int,
        direction::Tuple{<:Real,<:Real,<:Real};
        location::Tuple{<:Real,<:Real,<:Real}=(0, 0, 0)
    )

Construct a RotoinversionAxis object with order `n` about the axis specified by `direction`.

Keyword Arguments
=================
- `location`: a point on the rotoinversion axis in unit cell fractional coordinates
"""
function RotoinversionAxis(
    n::Int,
    direction::Tuple{<:Real,<:Real,<:Real};
    location::Tuple{<:Real,<:Real,<:Real}=(0, 0, 0),
)
    return RotoinversionAxis(n, direction, location)
end

# --- Functions/Methods

import Base.:(==)
import LinearAlgebra: dot

function Base.:(==)(x::RotoinversionAxis, y::RotoinversionAxis)
    # Check that orders are equal
    if x.n != y.n
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
