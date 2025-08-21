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
Tests for symmetry/centering.jl
"""
# --- Imports

# Standard library
using Test

# Xtallography package
using Xtallography

# --- Tests

# ------ Types

@testset "Centering subtypes" begin
    expected_types = [PrimitiveCentering, BaseCentering, BodyCentering, FaceCentering]
    for type in expected_types
        @test type <: Centering
    end
end

# ------ Constants

@testset "Centering constants" begin
    @test primitive_centering === PrimitiveCentering()
    @test P_centering === PrimitiveCentering()

    @test base_centering === BaseCentering()

    @test body_centering === BodyCentering()
    @test I_centering === BodyCentering()

    @test face_centering === FaceCentering()
    @test F_centering === FaceCentering()
end
