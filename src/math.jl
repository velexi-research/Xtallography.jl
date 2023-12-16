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
Math utility functions
"""
# --- Exports

# Functions
export is_basis, volume, surface_area

const LINEAR_INDEPENDENCE_ZERO = 1e-9

"""
    is_basis(
        v1::Vector{Real}, v2::Vector{Real}, v3::Vector{Real}
    ) -> Bool

Determine if the vectors `v1`, `v2`, and `v3` are a basis for a three-dimensional lattice.

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
