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
Tests for symmetry/lattice_system.jl
"""
# --- Imports

# Standard library
using Test

# Xtallography package
using Xtallography

# --- Tests

# ------ Types

@testset "LatticeSystem subtypes" begin
    expected_types = [
        Triclinic, Monoclinic, Orthorhombic, Tetragonal, Rhombohedral, Hexagonal, Cubic
    ]
    for type in expected_types
        @test type <: LatticeSystem
    end
end

# ------ Constants

@testset "LatticeSystem constants" begin
    @test triclinic === Triclinic()
    @test monoclinic === Monoclinic()
    @test orthorhombic === Orthorhombic()
    @test tetragonal === Tetragonal()
    @test rhombohedral === Rhombohedral()
    @test hexagonal === Hexagonal()
    @test cubic === Cubic()
end
