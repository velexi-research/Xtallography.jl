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
MirrorPlane type and functions
"""
# --- Exports

# Types
export MirrorPlane

# --- Types

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

# --- Functions/Methods

import Base.:(==)
import LinearAlgebra: dot

function Base.:(==)(x::MirrorPlane, y::MirrorPlane)
    # Check that directions of the plane normal vectors are the same
    if dot(x.normal, y.normal)^2 != dot(x.normal, x.normal) * dot(y.normal, y.normal)
        return false
    end

    # Check that line through both locations lies a plane orthogonal to the plane normal
    # vectors
    delta = (x.location[i] - y.location[i] for i in 1:3)
    return dot(delta, x.normal) == 0
end
