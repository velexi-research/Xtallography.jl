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
Tests for SymmetryElement type (defined in symmetry/symmetry_elements.jl)
"""
# --- Imports

# Standard library
using Test

# External packages
using Xtallography

# --- Tests

@testset "SymmetryElement Subtypes" begin
    expected_types = [
        RotationAxis, MirrorPlane, InversionCenter, RotoinversionAxis, GlidePlane, ScrewAxis
    ]

    for type in expected_types
        @test type <: SymmetryElement
    end
end

@testset ":(==)(::SymmetryElement,::SymmetryElement): mismatched types" begin
    # --- RotationAxis != MirrorPlane

    symmetry_element_1 = RotationAxis((1, 0, 0), (0, 0, 0), 2)
    symmetry_element_2 = MirrorPlane((1, 0, 0), (0, 0, 0))

    @test symmetry_element_1 != symmetry_element_2

    # --- RotationAxis != InversionCenter

    symmetry_element_1 = RotationAxis((1, 0, 0), (0, 0, 0), 2)
    symmetry_element_2 = InversionCenter((0, 0, 0))

    @test symmetry_element_1 != symmetry_element_2

    # --- RotationAxis != RotoinversionAxis

    symmetry_element_1 = RotationAxis((1, 0, 0), (0, 0, 0), 2)
    symmetry_element_2 = RotoinversionAxis((1, 0, 0), (0, 0, 0), 3)

    @test symmetry_element_1 != symmetry_element_2

    # --- RotationAxis != GlidePlane

    symmetry_element_1 = RotationAxis((1, 0, 0), (0, 0, 0), 2)
    symmetry_element_2 = GlidePlane("0,1,0", "1,0,0")

    @test symmetry_element_1 != symmetry_element_2

    # --- RotationAxis != ScrewAxis

    symmetry_element_1 = RotationAxis((1, 0, 0), (0, 0, 0), 2)
    symmetry_element_2 = ScrewAxis("1,0,0", 2, 1)

    @test symmetry_element_1 != symmetry_element_2

    # --- MirrorPlane != InversionCenter

    symmetry_element_1 = MirrorPlane((1, 0, 0), (0, 0, 0))
    symmetry_element_2 = InversionCenter((0, 0, 0))

    @test symmetry_element_1 != symmetry_element_2

    # --- MirrorPlane != RotoinversionAxis

    symmetry_element_1 = MirrorPlane((1, 0, 0), (0, 0, 0))
    symmetry_element_2 = RotoinversionAxis((1, 0, 0), (0, 0, 0), 3)

    @test symmetry_element_1 != symmetry_element_2

    # --- MirrorPlane != GlidePlane

    symmetry_element_1 = MirrorPlane((1, 0, 0), (0, 0, 0))
    symmetry_element_2 = GlidePlane("0,1,0", "1,0,0")

    @test symmetry_element_1 != symmetry_element_2

    # --- MirrorPlane != ScrewAxis

    symmetry_element_1 = MirrorPlane((1, 0, 0), (0, 0, 0))
    symmetry_element_2 = ScrewAxis("1,0,0", 2, 1)

    @test symmetry_element_1 != symmetry_element_2

    # --- InversionCenter != RotoinversionAxis

    symmetry_element_1 = InversionCenter((0, 0, 0))
    symmetry_element_2 = RotoinversionAxis((1, 0, 0), (0, 0, 0), 3)

    @test symmetry_element_1 != symmetry_element_2

    # --- InversionCenter != GlidePlane

    symmetry_element_1 = InversionCenter((0, 0, 0))
    symmetry_element_2 = GlidePlane("0,1,0", "1,0,0")

    @test symmetry_element_1 != symmetry_element_2

    # --- InversionCenter != ScrewAxis

    symmetry_element_1 = InversionCenter((0, 0, 0))
    symmetry_element_2 = ScrewAxis("1,0,0", 2, 1)

    @test symmetry_element_1 != symmetry_element_2

    # --- RotoinversionAxis != GlidePlane

    symmetry_element_1 = RotoinversionAxis((1, 0, 0), (0, 0, 0), 3)
    symmetry_element_2 = GlidePlane("0,1,0", "1,0,0")

    @test symmetry_element_1 != symmetry_element_2

    # --- RotoinversionAxis != ScrewAxis

    symmetry_element_1 = RotoinversionAxis((1, 0, 0), (0, 0, 0), 3)
    symmetry_element_2 = ScrewAxis("1,0,0", 2, 1)

    @test symmetry_element_1 != symmetry_element_2

    # --- GlidePlane != ScrewAxis

    symmetry_element_1 = GlidePlane("0,1,0", "1,0,0")
    symmetry_element_2 = ScrewAxis("1,0,0", 2, 1)

    @test symmetry_element_1 != symmetry_element_2
end
