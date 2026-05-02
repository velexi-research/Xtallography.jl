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
Constants and functions related to Bravais lattices
"""
# --- Exports

# ------ Functions/Methods

export is_bravais_lattice

# ------ Constants

export BRAVAIS_LATTICES

# --- Constants

"""
    BRAVAIS_LATTICES

List of valid Bravais lattices.
"""
const BRAVAIS_LATTICES = (
    (lattice_system=triclinic, centering=primitive_centering),
    (lattice_system=monoclinic, centering=primitive_centering),
    (lattice_system=monoclinic, centering=body_centering),
    (lattice_system=monoclinic, centering=base_centering),
    (lattice_system=orthorhombic, centering=primitive_centering),
    (lattice_system=orthorhombic, centering=body_centering),
    (lattice_system=orthorhombic, centering=face_centering),
    (lattice_system=orthorhombic, centering=base_centering),
    (lattice_system=tetragonal, centering=primitive_centering),
    (lattice_system=tetragonal, centering=body_centering),
    (lattice_system=rhombohedral, centering=primitive_centering),
    (lattice_system=hexagonal, centering=primitive_centering),
    (lattice_system=cubic, centering=primitive_centering),
    (lattice_system=cubic, centering=body_centering),
    (lattice_system=cubic, centering=face_centering),
)

# --- Functions/Methods

"""
    is_bravais_lattice(lattice_system::LatticeSystem, centering::Centering) -> Bool

Determine if `lattice_system` and `centering` define a valid Bravais lattice type.

Return values
=============
- `true` if `lattice_system` and `centering` define a valid Bravais lattice type; `false`
  otherwise

Examples
========
```jldoctest
julia> is_bravais_lattice(cubic, body_centering)
true

julia> is_bravais_lattice(cubic, base_centering)
false
```
"""
@inline function is_bravais_lattice(lattice_system::LatticeSystem, centering::Centering)
    return (lattice_system=lattice_system, centering=centering) in BRAVAIS_LATTICES
end
