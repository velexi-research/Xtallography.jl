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
core_tests.jl defines tests for lattices/core.jl
"""
# --- Imports

# Standard library
using LinearAlgebra: dot, qr, I
using Test

# External packages
using Combinatorics: permutations

# XtallographyUtils package
using XtallographyUtils

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

@testset "Centering subtypes" begin
    expected_types = [Primitive, BaseCentered, BodyCentered, FaceCentered]
    for type in expected_types
        @test type <: Centering
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

@testset "Centering constants" begin
    @test primitive === Primitive()
    @test base_centered === BaseCentered()
    @test body_centered === BodyCentered()
    @test face_centered === FaceCentered()
end

@testset "BRAVAIS_LATTICES" begin
    @test BRAVAIS_LATTICES isa Tuple
    @test length(BRAVAIS_LATTICES) == 15
end

# ------ Methods

@testset "is_bravais_lattice(::LatticeSystem, ::Centering)" begin
    # --- Tests

    # cubic
    @test is_bravais_lattice(cubic, primitive)
    @test is_bravais_lattice(cubic, body_centered)
    @test is_bravais_lattice(cubic, face_centered)
    @test !is_bravais_lattice(cubic, base_centered)

    # tetragonal
    @test is_bravais_lattice(tetragonal, primitive)
    @test is_bravais_lattice(tetragonal, body_centered)
    @test !is_bravais_lattice(tetragonal, face_centered)
    @test !is_bravais_lattice(tetragonal, base_centered)

    # orthorhombic
    @test is_bravais_lattice(orthorhombic, primitive)
    @test is_bravais_lattice(orthorhombic, body_centered)
    @test is_bravais_lattice(orthorhombic, face_centered)
    @test is_bravais_lattice(orthorhombic, base_centered)

    # hexagonal
    @test is_bravais_lattice(hexagonal, primitive)
    @test !is_bravais_lattice(hexagonal, body_centered)
    @test !is_bravais_lattice(hexagonal, face_centered)
    @test !is_bravais_lattice(hexagonal, base_centered)

    # rhombohedral
    @test is_bravais_lattice(rhombohedral, primitive)
    @test !is_bravais_lattice(rhombohedral, body_centered)
    @test !is_bravais_lattice(rhombohedral, face_centered)
    @test !is_bravais_lattice(rhombohedral, base_centered)

    # monoclinic
    @test is_bravais_lattice(monoclinic, primitive)
    @test is_bravais_lattice(monoclinic, body_centered)
    @test !is_bravais_lattice(monoclinic, face_centered)
    @test is_bravais_lattice(monoclinic, base_centered)

    # triclinic
    @test is_bravais_lattice(triclinic, primitive)
    @test !is_bravais_lattice(triclinic, body_centered)
    @test !is_bravais_lattice(triclinic, face_centered)
    @test !is_bravais_lattice(triclinic, base_centered)
end
