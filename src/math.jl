#   Copyright 2023 Velexi Corporation
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
Math utility functions and constants
"""
# --- Exports

# Functions
export asin_, acos_
export is_basis, volume, surface_area

# --- Constants

# Numerical comparison constants
const ALMOST_ZERO = 1e-9
const COS_APPROX_ZERO = 1e-9
const LINEAR_INDEPENDENCE_ZERO = 1e-9

# Common trigonometric computations
const SIN_PI_OVER_THREE = sqrt(3) / 2
const SIN_PI_OVER_FOUR = 1 / sqrt(2)
const ACOS_MINUS_ONE_THIRD = acos(-1 / 3)

# --- Functions/Methods

# ------ basic mathematical functions

"""
    asin_(x::Real; rtol::Real) -> Float64

Compute the arcsin of `x` with a tolerance for values of `x` that are slightly outside of
the mathematical domain [-1, 1].

When the value of `x` is approximately equal to 1, `asin_(x)` returns π / 2; when the value
of `x` is approximately equal to -1, `asin_(x)` returns -π / 2.

Return values
=============
- arcsin of `x`

Examples
========
```jldoctest
julia> asin_(0.5) ≈ π / 6
true
julia> asin_(1 + eps(1.0)) == π / 2
true
julia> asin_(-1 - eps(1.0)) == -π / 2
true
```
"""
function asin_(x::Real; rtol::Real=√eps(1.0))
    if x > 1 && isapprox(x, 1; rtol=rtol)
        return π / 2
    elseif x < -1 && isapprox(x, -1; rtol=rtol)
        return -π / 2
    else
        return asin(x)
    end
end

"""
    acos_(x::Real; rtol::Real) -> Float64

Compute the arccos of `x` with a tolerance for values of `x` that are slightly outside of
the mathematical domain [-1, 1].

When the value of `x` is approximately equal to 1, `acos_(x)` returns 0; when the value
of `x` is approximately equal to -1, `acos_(x)` returns π.

Return values
=============
- arccos of `x`

Examples
========
```jldoctest
julia> acos_(0.5) ≈ π / 3
true
julia> acos_(1 + eps(1.0)) == 0
true
julia> acos_(-1 - eps(1.0)) == π
true
```
"""
function acos_(x::Real; rtol::Real=√eps(1.0))
    if x > 1 && isapprox(x, 1; rtol=rtol)
        return 0
    elseif x < -1 && isapprox(x, -1; rtol=rtol)
        return π
    else
        return acos(x)
    end
end

# ------ linear algebra

"""
    is_basis(
        v1::Vector{Real}, v2::Vector{Real}, v3::Vector{Real}
    ) -> Bool

Determine if the vectors `v1`, `v2`, and `v3` are a basis for a three-dimensional lattice
(i.e., `v1`, `v2`, and `v3` are linearly independent).

Return values
=============
- `true` if `v1`, `v2`, and `v3` are a basis; `false` otherwise

Examples
========
TODO
"""
function is_basis(v1::Vector{<:Real}, v2::Vector{<:Real}, v3::Vector{<:Real})
    return volume(v1, v2, v3) > norm(v1) * norm(v2) * norm(v3) * LINEAR_INDEPENDENCE_ZERO
end

# ------ geometry

function volume(basis_a::T, basis_b::T, basis_c::T) where {T<:Vector{<:Real}}
    return abs(det(reduce(hcat, [basis_a, basis_b, basis_c])))
end

function surface_area(basis_a::T, basis_b::T, basis_c::T) where {T<:Vector{<:Real}}

    # Compute the cross products of the basis vectors
    cross_ab = cross(basis_a, basis_b)
    cross_ac = cross(basis_a, basis_c)
    cross_bc = cross(basis_b, basis_c)

    # Compute the surface area
    return 2 * (
        sqrt(dot(cross_ab, cross_ab)) +
        sqrt(dot(cross_ac, cross_ac)) +
        sqrt(dot(cross_bc, cross_bc))
    )
end
