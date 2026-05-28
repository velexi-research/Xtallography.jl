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
export RotationAxis, MirrorPlane, InversionCenter, RotoinversionAxis
export GlidePlane, ScrewAxis

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

import Base.:(==)
import LinearAlgebra: dot

# ------ SymmetryElement

"""
    SymmetryElement

Supertype for symmetry elements in 3D

Subtypes
========
[`GlidePlane`](@ref), [`ScrewAxis`](@ref)
"""
abstract type SymmetryElement end

# ------ RotationAxis

"""
    RotationAxis

Type representing a rotation axis

Supertype: [`SymmetryElement`](@ref)
"""
struct RotationAxis <: SymmetryElement
    # --- Fields

    # direction of rotation axis in unit cell fractional coordinates
    direction::Tuple{Rational,Rational,Rational}

    # a point on the rotation axis in unit cell fractional coordinates
    location::Tuple{Rational,Rational,Rational}

    # rotation order
    n::Int
end

function Base.:(==)(x::RotationAxis, y::RotationAxis)
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

# ------ MirrorPlane

"""
    MirrorPlane

Type representing a reflection through a plane

Supertype: [`SymmetryElement`](@ref)
"""
struct MirrorPlane <: SymmetryElement
    # --- Fields

    # normal to mirror plane
    normal::Tuple{Rational,Rational,Rational}

    # a point on the mirror plane in unit cell fractional coordinates
    location::Tuple{Rational,Rational,Rational}
end

function Base.:(==)(x::MirrorPlane, y::MirrorPlane)
    # Check that directions of the plane normal are the same
    if dot(x.normal, y.normal)^2 != dot(x.normal, x.normal) * dot(y.normal, y.normal)
        return false
    end

    # Check that line through both locations lies a plane orthogonal to the normal vectors
    delta = (x.location[i] - y.location[i] for i in 1:3)
    return dot(delta, x.normal) == 0
end

# ------ InversionCenter

"""
    InversionCenter

Type representing an inversion through a point

Supertype: [`SymmetryElement`](@ref)
"""
struct InversionCenter <: SymmetryElement
    # --- Fields

    # location of inversion center
    center::Tuple{Rational,Rational,Rational}
end

# ------ RotoinversionAxis

"""
    RotoinversionAxis

Type representing a rotoinversion axis

Supertype: [`SymmetryElement`](@ref)
"""
struct RotoinversionAxis <: SymmetryElement
    # --- Fields

    # direction of rotoinversion axis in unit cell fractional coordinates
    direction::Tuple{Rational,Rational,Rational}

    # a point on the rotation axis in unit cell fractional coordinates
    location::Tuple{Rational,Rational,Rational}

    # rotation order
    n::Int
end

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

# ------ GlidePlane

"""
    GlidePlane

Type representing a glide plane

Supertype: [`SymmetryElement`](@ref)
"""
struct GlidePlane <: SymmetryElement
    # --- Fields

    # glide translation
    translation::String

    # reflection plane
    reflection_plane::String
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
    # --- Fields

    # direction of rotation axis
    axis::String

    # rotation order
    n::Int

    # number of translation steps of size 1/n following rotation by 2π/n
    m::Int

    # Constructor
    function ScrewAxis(axis::String, n::Int, m::Int)

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
