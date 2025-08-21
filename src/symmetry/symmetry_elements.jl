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
Symmetry element types and functions
"""
# --- Exports

# ------ Types

export SymmetryElement
export GlidePlane
export ScrewAxis

# ------ Constants

# Glide planes
export b_perp_a, c_perp_a, n_perp_a, d_perp_a
export a_perp_b, c_perp_b, n_perp_b, d_perp_b
export a_perp_c, b_perp_c, n_perp_c, d_perp_c
export c_perp_110, c_perp_120, d_perp_110

# Screw axes
export a_2_1, a_4_1, a_4_2, a_4_3
export b_2_1, b_4_1, b_4_2, b_4_3
export c_2_1, c_4_1, c_4_2, c_4_3
export c_3_1, c_3_2
export c_6_1, c_6_2, c_6_3, c_6_4, c_6_5

# --- Types

# ------ SymmetryElement

"""
    SymmetryElement

Supertype for symmetry elements in 3D

Subtypes
========
[`GlidePlane`](@ref), [`ScrewAxis`](@ref)
"""
abstract type SymmetryElement end

# ------ GlidePlane

"""
    GlidePlane

Type representing a glide plane

Supertype: [`SymmetryElement`](@ref)
"""
struct GlidePlane <: SymmetryElement
    # Fields
    translation::String  # glide translation
    reflection_plane::String  # reflection plane
end

# Constants
const b_perp_a = GlidePlane("0,1,0", "1,0,0")
const c_perp_a = GlidePlane("0,0,1", "1,0,0")
const n_perp_a = GlidePlane("0,1/2,1/2", "1,0,0")
const d_perp_a = GlidePlane("0,1/4,1/4", "1,0,0")

const a_perp_b = GlidePlane("1,0,0", "0,1,0")
const c_perp_b = GlidePlane("0,0,1", "0,1,0")
const n_perp_b = GlidePlane("1/2,0,1/2", "0,1,0")
const d_perp_b = GlidePlane("1/4,0,1/4", "0,1,0")

const a_perp_c = GlidePlane("1,0,0", "0,0,1")
const b_perp_c = GlidePlane("0,1,0", "0,0,1")
const n_perp_c = GlidePlane("1/2,1/2,0", "0,0,1")
const d_perp_c = GlidePlane("1/4,1/4,0", "0,0,1")

const c_perp_110 = GlidePlane("0,0,1", "1,1,0")
const c_perp_120 = GlidePlane("0,0,1", "1,2,0")
const d_perp_110 = GlidePlane("1/4,1/4,1/4", "1,1,0")

# ------ ScrewAxis

"""
    ScrewAxis

Type representing a screw axis

Supertype: [`SymmetryElement`](@ref)
"""
struct ScrewAxis <: SymmetryElement
    # Fields
    axis::String  # direction of rotation axis
    n::Int  # rotation order
    m::Int  # number of translation steps of size 1/n following rotation by 2π/n

    # Constructor
    function ScrewAxis(axis::String, n::Int, m::Int)

        # Enforce constraints
        if m > n
            throw(ArgumentError("`m` be no greater than `n` (n=$n,m=$m)"))
        end

        # Return new ScrewAxis
        return new(axis, n, m)
    end
end

# Constants
const a_2_1 = ScrewAxis("1,0,0", 2, 1)
const a_4_1 = ScrewAxis("1,0,0", 4, 1)
const a_4_2 = ScrewAxis("1,0,0", 4, 2)
const a_4_3 = ScrewAxis("1,0,0", 4, 3)

const b_2_1 = ScrewAxis("0,1,0", 2, 1)
const b_4_1 = ScrewAxis("0,1,0", 4, 1)
const b_4_2 = ScrewAxis("0,1,0", 4, 2)
const b_4_3 = ScrewAxis("0,1,0", 4, 3)

const c_2_1 = ScrewAxis("0,0,1", 2, 1)
const c_3_1 = ScrewAxis("0,0,1", 3, 1)
const c_3_2 = ScrewAxis("0,0,1", 3, 2)
const c_4_1 = ScrewAxis("0,0,1", 4, 1)
const c_4_2 = ScrewAxis("0,0,1", 4, 2)
const c_4_3 = ScrewAxis("0,0,1", 4, 3)
const c_6_1 = ScrewAxis("0,0,1", 6, 1)
const c_6_2 = ScrewAxis("0,0,1", 6, 2)
const c_6_3 = ScrewAxis("0,0,1", 6, 3)
const c_6_4 = ScrewAxis("0,0,1", 6, 4)
const c_6_5 = ScrewAxis("0,0,1", 6, 5)
