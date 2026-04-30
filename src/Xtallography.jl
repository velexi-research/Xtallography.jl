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
The Xtallography package provides utility functions for analysis of 3D crystals.
"""
module Xtallography

# --- Package Metadata

using TOML: TOML
const VERSION = TOML.parsefile(joinpath(pkgdir(@__MODULE__), "Project.toml"))["version"]

# --- Core types and methods

include("math.jl")

# --- Symmetry types and methods

include("symmetry/lattice_systems.jl")
include("symmetry/centerings.jl")
include("symmetry/symmetry_elements.jl")
include("symmetry/bravais_lattices.jl")

# --- Unit cell types and methods

include("unit_cell/unit_cell_symmetry.jl")
include("unit_cell/unit_cell.jl")
include("unit_cell/unit_cell_delta.jl")
include("unit_cell/triclinic.jl")
include("unit_cell/monoclinic.jl")
include("unit_cell/orthorhombic.jl")
include("unit_cell/tetragonal.jl")
include("unit_cell/rhombohedral.jl")
include("unit_cell/hexagonal.jl")
include("unit_cell/cubic.jl")

end  # End of Xtallography module
