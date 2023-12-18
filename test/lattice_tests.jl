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
lattice_tests.jl defines tests for lattice.jl
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

@testset "LatticeSystem Subtypes" begin
    expected_types = [
        Triclinic, Monoclinic, Orthorhombic, Tetragonal, Rhombohedral, Hexagonal, Cubic
    ]
    for type in expected_types
        @test type <: LatticeSystem
    end
end

@testset "Centering Subtypes" begin
    expected_types = [Primitive, BaseCentered, BodyCentered, FaceCentered]
    for type in expected_types
        @test type <: Centering
    end
end

# ------ Constnats

@testset "BRAVAIS_LATTICES" begin
    @test BRAVAIS_LATTICES isa Tuple
    @test length(BRAVAIS_LATTICES) == 15
end

# ------ Methods

@testset "is_bravais_lattice(::LatticeSystem, ::Centering)" begin
    # --- Tests

    # Cubic
    @test is_bravais_lattice(Cubic(), Primitive())
    @test is_bravais_lattice(Cubic(), BodyCentered())
    @test is_bravais_lattice(Cubic(), FaceCentered())
    @test !is_bravais_lattice(Cubic(), BaseCentered())

    # Tetragonal
    @test is_bravais_lattice(Tetragonal(), Primitive())
    @test is_bravais_lattice(Tetragonal(), BodyCentered())
    @test !is_bravais_lattice(Tetragonal(), FaceCentered())
    @test !is_bravais_lattice(Tetragonal(), BaseCentered())

    # Orthorhombic
    @test is_bravais_lattice(Orthorhombic(), Primitive())
    @test is_bravais_lattice(Orthorhombic(), BodyCentered())
    @test is_bravais_lattice(Orthorhombic(), FaceCentered())
    @test is_bravais_lattice(Orthorhombic(), BaseCentered())

    # Hexagonal
    @test is_bravais_lattice(Hexagonal(), Primitive())
    @test !is_bravais_lattice(Hexagonal(), BodyCentered())
    @test !is_bravais_lattice(Hexagonal(), FaceCentered())
    @test !is_bravais_lattice(Hexagonal(), BaseCentered())

    # Rhombohedral
    @test is_bravais_lattice(Rhombohedral(), Primitive())
    @test !is_bravais_lattice(Rhombohedral(), BodyCentered())
    @test !is_bravais_lattice(Rhombohedral(), FaceCentered())
    @test !is_bravais_lattice(Rhombohedral(), BaseCentered())

    # Monoclinic
    @test is_bravais_lattice(Monoclinic(), Primitive())
    @test is_bravais_lattice(Monoclinic(), BodyCentered())
    @test !is_bravais_lattice(Monoclinic(), FaceCentered())
    @test is_bravais_lattice(Monoclinic(), BaseCentered())

    # Triclinic
    @test is_bravais_lattice(Triclinic(), Primitive())
    @test !is_bravais_lattice(Triclinic(), BodyCentered())
    @test !is_bravais_lattice(Triclinic(), FaceCentered())
    @test !is_bravais_lattice(Triclinic(), BaseCentered())
end

@testset "is_bravais_lattice(::Type{<:LatticeSystem}, ::Centering)" begin
    # --- Tests

    # Cubic
    @test is_bravais_lattice(Cubic, Primitive())
    @test is_bravais_lattice(Cubic, BodyCentered())
    @test is_bravais_lattice(Cubic, FaceCentered())
    @test !is_bravais_lattice(Cubic, BaseCentered())

    # Tetragonal
    @test is_bravais_lattice(Tetragonal, Primitive())
    @test is_bravais_lattice(Tetragonal, BodyCentered())
    @test !is_bravais_lattice(Tetragonal, FaceCentered())
    @test !is_bravais_lattice(Tetragonal, BaseCentered())

    # Orthorhombic
    @test is_bravais_lattice(Orthorhombic, Primitive())
    @test is_bravais_lattice(Orthorhombic, BodyCentered())
    @test is_bravais_lattice(Orthorhombic, FaceCentered())
    @test is_bravais_lattice(Orthorhombic, BaseCentered())

    # Hexagonal
    @test is_bravais_lattice(Hexagonal, Primitive())
    @test !is_bravais_lattice(Hexagonal, BodyCentered())
    @test !is_bravais_lattice(Hexagonal, FaceCentered())
    @test !is_bravais_lattice(Hexagonal, BaseCentered())

    # Rhombohedral
    @test is_bravais_lattice(Rhombohedral, Primitive())
    @test !is_bravais_lattice(Rhombohedral, BodyCentered())
    @test !is_bravais_lattice(Rhombohedral, FaceCentered())
    @test !is_bravais_lattice(Rhombohedral, BaseCentered())

    # Monoclinic
    @test is_bravais_lattice(Monoclinic, Primitive())
    @test is_bravais_lattice(Monoclinic, BodyCentered())
    @test !is_bravais_lattice(Monoclinic, FaceCentered())
    @test is_bravais_lattice(Monoclinic, BaseCentered())

    # Triclinic
    @test is_bravais_lattice(Triclinic, Primitive())
    @test !is_bravais_lattice(Triclinic, BodyCentered())
    @test !is_bravais_lattice(Triclinic, FaceCentered())
    @test !is_bravais_lattice(Triclinic, BaseCentered())
end
