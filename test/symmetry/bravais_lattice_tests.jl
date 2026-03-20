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
Tests for symmetry/bravais_lattices.jl
"""
# --- Imports

# Standard library
using Test

# Xtallography package
using Xtallography

# --- Tests

# ------ Constants

@testset "BRAVAIS_LATTICES" begin
    @test BRAVAIS_LATTICES isa Tuple
    @test length(BRAVAIS_LATTICES) == 15
end

# ------ Methods

@testset "is_bravais_lattice(::LatticeSystem, ::Centering)" begin
    # --- Tests

    # cubic
    @test is_bravais_lattice(cubic, primitive_centering)
    @test is_bravais_lattice(cubic, body_centering)
    @test is_bravais_lattice(cubic, face_centering)
    @test !is_bravais_lattice(cubic, base_centering)

    # tetragonal
    @test is_bravais_lattice(tetragonal, primitive_centering)
    @test is_bravais_lattice(tetragonal, body_centering)
    @test !is_bravais_lattice(tetragonal, face_centering)
    @test !is_bravais_lattice(tetragonal, base_centering)

    # orthorhombic
    @test is_bravais_lattice(orthorhombic, primitive_centering)
    @test is_bravais_lattice(orthorhombic, body_centering)
    @test is_bravais_lattice(orthorhombic, face_centering)
    @test is_bravais_lattice(orthorhombic, base_centering)

    # hexagonal
    @test is_bravais_lattice(hexagonal, primitive_centering)
    @test !is_bravais_lattice(hexagonal, body_centering)
    @test !is_bravais_lattice(hexagonal, face_centering)
    @test !is_bravais_lattice(hexagonal, base_centering)

    # rhombohedral
    @test is_bravais_lattice(rhombohedral, primitive_centering)
    @test !is_bravais_lattice(rhombohedral, body_centering)
    @test !is_bravais_lattice(rhombohedral, face_centering)
    @test !is_bravais_lattice(rhombohedral, base_centering)

    # monoclinic
    @test is_bravais_lattice(monoclinic, primitive_centering)
    @test is_bravais_lattice(monoclinic, body_centering)
    @test !is_bravais_lattice(monoclinic, face_centering)
    @test is_bravais_lattice(monoclinic, base_centering)

    # triclinic
    @test is_bravais_lattice(triclinic, primitive_centering)
    @test !is_bravais_lattice(triclinic, body_centering)
    @test !is_bravais_lattice(triclinic, face_centering)
    @test !is_bravais_lattice(triclinic, base_centering)
end
