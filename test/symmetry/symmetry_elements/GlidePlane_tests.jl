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
Tests for GlidePlane type (defined in symmetry/symmetry_elements.jl)
"""
# --- Imports

# Standard library
using Test

# External packages
using Xtallography

# --- Tests

@testset "GlidePlane constructor" begin
    # --- Tests

    # ------ Valid arguments

    translation = "1,1,1"
    reflection_plane = "0,1,0"
    symmetry_element = GlidePlane(translation, reflection_plane)

    @test symmetry_element.translation == translation
    @test symmetry_element.reflection_plane == reflection_plane
end

@testset "GlidePlane constants" begin
    expected_constants = [
        b_perp_a,
        c_perp_a,
        n_perp_a,
        d_perp_a,
        a_perp_b,
        c_perp_b,
        n_perp_b,
        d_perp_b,
        a_perp_c,
        b_perp_c,
        n_perp_c,
        d_perp_c,
    ]
    for constant in expected_constants
        @test constant isa GlidePlane
    end
end
