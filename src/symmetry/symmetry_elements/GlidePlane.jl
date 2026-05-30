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
GlidePlane type and functions
"""
# --- Exports

# ------ Types

export GlidePlane

# Constants
export b_perp_a, c_perp_a, n_perp_a, d_perp_a
export a_perp_b, c_perp_b, n_perp_b, d_perp_b
export a_perp_c, b_perp_c, n_perp_c, d_perp_c
export c_perp_110, c_perp_120, d_perp_110

# --- Types

import LinearAlgebra: dot

"""
    GlidePlane

Type representing a glide plane

Supertype: [`SymmetryElement`](@ref)
"""
struct GlidePlane <: SymmetryElement
    # --- Fields

    # glide direction
    glide::Tuple{Rational,Rational,Rational}

    # normal to mirror plane
    normal::Tuple{Rational,Rational,Rational}

    # a point on the mirror plane in unit cell fractional coordinates
    location::Tuple{Rational,Rational,Rational}

    # --- Constructor

    function GlidePlane(
        glide::Tuple{<:Real,<:Real,<:Real},
        normal::Tuple{<:Real,<:Real,<:Real},
        location::Tuple{<:Real,<:Real,<:Real},
    )

        # Enforce constraints
        if dot(glide, normal) != 0
            throw(
                ArgumentError(
                    "`glide` must be orthogonal to `normal` " *
                    "(glide=$glide,normal=$normal)",
                ),
            )
        end

        # Return new GlidePlane
        return new(glide, normal, location)
    end
end

# Outer constructor
"""
    GlidePlane(
        glide::Tuple{<:Real,<:Real,<:Real},
        normal::Tuple{<:Real,<:Real,<:Real};
        location::Tuple{<:Real,<:Real,<:Real}=(0, 0, 0)
    )

Construct a GlidePlane object that

* has an axis by `direction` and `location`

and

* couples a rotation by `2π/n` and a translation along the rotation axis by `m/n` (in unit
  cell fractional coordinates).

Keyword Arguments
=================
- `location`: a point on the rotation axis in unit cell fractional coordinates
"""
function GlidePlane(
    glide::Tuple{<:Real,<:Real,<:Real},
    normal::Tuple{<:Real,<:Real,<:Real};
    location::Tuple{<:Real,<:Real,<:Real}=(0, 0, 0),
)
    return GlidePlane(glide, normal, location)
end

# --- Functions/Methods

import Base.:(==)
import LinearAlgebra: dot

function Base.:(==)(x::GlidePlane, y::GlidePlane)
    # Check that directions of the glides are the same
    if dot(x.glide, y.glide)^2 != dot(x.glide, x.glide) * dot(y.glide, y.glide)
        return false
    end

    # Check that directions of the plane normal vectors are the same
    if dot(x.normal, y.normal)^2 != dot(x.normal, x.normal) * dot(y.normal, y.normal)
        return false
    end

    # Check that line through both locations lies a plane orthogonal to the plane normal
    # vectors
    delta = (x.location[i] - y.location[i] for i in 1:3)
    return dot(delta, x.normal) == 0
end

# --- Constants

# TODO: check d and n glide planes
const b_perp_a = GlidePlane((0, 1, 0), (1, 0, 0))
const c_perp_a = GlidePlane((0, 0, 1), (1, 0, 0))
const n_perp_a = GlidePlane((0, 1/2, 1/2), (1, 0, 0))
const d_perp_a = GlidePlane((0, 1/4, 1/4), (1, 0, 0))

const a_perp_b = GlidePlane((1, 0, 0), (0, 1, 0))
const c_perp_b = GlidePlane((0, 0, 1), (0, 1, 0))
const n_perp_b = GlidePlane((1/2, 0, 1/2), (0, 1, 0))
const d_perp_b = GlidePlane((1/4, 0, 1/4), (0, 1, 0))

const a_perp_c = GlidePlane((1, 0, 0), (0, 0, 1))
const b_perp_c = GlidePlane((0, 1, 0), (0, 0, 1))
const n_perp_c = GlidePlane((1/2, 1/2, 0), (0, 0, 1))
const d_perp_c = GlidePlane((1/4, 1/4, 0), (0, 0, 1))

const c_perp_110 = GlidePlane((0, 0, 1), (1, 1, 0))
const c_perp_120 = GlidePlane((0, 0, 1), (1, 2, 0))
const d_perp_110 = GlidePlane((1/4, -1/4, 0), (1, 1, 0))
