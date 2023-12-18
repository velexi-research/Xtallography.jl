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
The XtallographyUtils package provides utility functions for analysis of 3D crystals.
"""
module XtallographyUtils

# --- Package Metadata

using TOML: TOML
const VERSION = TOML.parsefile(joinpath(pkgdir(@__MODULE__), "Project.toml"))["version"]

# --- Core types and methods

include("constants.jl")
include("math.jl")
include("lattice.jl")
include("unit_cell.jl")

# --- Lattice-specific types and methods

include("lattices/triclinic.jl")
include("lattices/monoclinic.jl")
include("lattices/orthorhombic.jl")
include("lattices/tetragonal.jl")
include("lattices/rhombohedral.jl")
include("lattices/hexagonal.jl")
include("lattices/cubic.jl")

# --- Constants

# Lattice Types
const BRAVAIS_LATTICES = [
    (lattice=Triclinic, centering=Primitive()),
    (lattice=Monoclinic, centering=Primitive()),
    (lattice=Monoclinic, centering=BodyCentered()),
    (lattice=Monoclinic, centering=BaseCentered()),
    (lattice=Orthorhombic, centering=Primitive()),
    (lattice=Orthorhombic, centering=BodyCentered()),
    (lattice=Orthorhombic, centering=FaceCentered()),
    (lattice=Orthorhombic, centering=BaseCentered()),
    (lattice=Tetragonal, centering=Primitive()),
    (lattice=Tetragonal, centering=BodyCentered()),
    (lattice=Rhombohedral, centering=Primitive()),
    (lattice=Hexagonal, centering=Primitive()),
    (lattice=Cubic, centering=Primitive()),
    (lattice=Cubic, centering=BodyCentered()),
    (lattice=Cubic, centering=FaceCentered()),
]

end  # End of XtallographyUtils module
