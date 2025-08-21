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
Centering types and functions
"""
# --- Exports

# ------ Types

# Centering
export Centering, PrimitiveCentering, BaseCentering, BodyCentering, FaceCentering
export primitive_centering, P_centering
export base_centering
export body_centering, I_centering
export face_centering, F_centering

# --- Types

"""
    Centering

Supertype for the four lattice centerings in 3D

Subtypes
========
[`PrimitiveCentering`](@ref), [`BaseCentering`](@ref), [`BodyCentering`](@ref),
[`FaceCentering`](@ref)
"""
abstract type Centering end

"""
    PrimitiveCentering

Type representing no centering that is the type of [`primitive_centering`](@ref) and
[`P_centering`](@ref)

Supertype: [`Centering`](@ref)
"""
struct PrimitiveCentering <: Centering end

"""
    primitive_centering, P_centering

The singleton instance of type [`PrimitiveCentering`](@ref)
"""
const primitive_centering = PrimitiveCentering()
const P_centering = primitive_centering

"""
    BaseCentering

Type representing base centering that is the type of [`base_centering`](@ref)

!!! note

    By convention, base-centering is

    * on the C-face of the unit cell for orthorhombic lattice systems

    and

    * on the B-face of the unit cell for monoclinic lattice systems.

Supertype: [`Centering`](@ref)
"""
struct BaseCentering <: Centering end

"""
    base_centering

The singleton instance of type [`BaseCentering`](@ref)
"""
const base_centering = BaseCentering()

"""
    BodyCentering

Type representing body centering that is the type of [`body_centering`](@ref) and
[`I_centering`](@ref)

Supertype: [`Centering`](@ref)
"""
struct BodyCentering <: Centering end

"""
    body_centering, I_centering

The singleton instance of type [`BodyCentering`](@ref)
"""
const body_centering = BodyCentering()
const I_centering = body_centering

"""
    FaceCentering

Type representing face centering that is the type of [`face_centering`](@ref) and
[`F_centering`](@ref)

Supertype: [`Centering`](@ref)
"""
struct FaceCentering <: Centering end

"""
    face_centering, F_centering

The singleton instance of type [`FaceCentering`](@ref)
"""
const face_centering = FaceCentering()
const F_centering = face_centering
