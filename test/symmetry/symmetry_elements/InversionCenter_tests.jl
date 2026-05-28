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
Tests for InversionCenter type (defined in symmetry/symmetry_elements.jl)
"""
# --- Imports

# Standard library
using Test

# External packages
using Xtallography

# --- Tests

@testset "InversionCenter constructor" begin
    # --- Valid arguments

    center = (1, 0, 0)

    symmetry_element = InversionCenter(center)

    @test symmetry_element.center == center

    # --- Invalid arguments

    # center does not contain exactly 3 elements
    center = (1, 0, 0, 0)

    @test_throws MethodError InversionCenter(center)
end

@testset ":(==)(::InversionCenter,::InversionCenter)" begin
    # --- Identical inversion centers

    symmetry_element_1 = InversionCenter((1, 0, 0))
    symmetry_element_2 = InversionCenter((1, 0, 0))

    @test symmetry_element_1 == symmetry_element_2

    # --- Inequivalent inversion centers: centers differ

    symmetry_element_1 = InversionCenter((1, 1, 0))
    symmetry_element_2 = InversionCenter((1, 0, 0))

    @test symmetry_element_1 != symmetry_element_2
end
