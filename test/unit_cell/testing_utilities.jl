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
Testing utilities
"""
# --- Imports

# Standard library
using LinearAlgebra: qr, I
using Test

# External packages
using Combinatorics: permutations

# Xtallography package
using Xtallography

# --- Utility functions

# Helper function for testing UnitCell(::Vector, ::Vector, ::Vector) constructor
function test_basis_rotations_and_permutations(
    rotations::Vector,
    expected_unit_cell::UnitCell,
    basis_a::Vector{<:Real},
    basis_b::Vector{<:Real},
    basis_c::Vector{<:Real};
    centering=nothing,
)
    # --- Preparations

    # Get type of expected unit cell
    expected_type = typeof(expected_unit_cell)

    # --- Test UnitCell(::Vector, ::Vector, ::Vector) for rotations of all basis
    #     permutations

    for rotation in rotations
        # Generate bases to check
        if centering == base_centering
            # Do not permute basis vectors for base-centering
            bases_to_test = ([basis_a, basis_b, basis_c],)
        else
            bases_to_test = permutations([basis_a, basis_b, basis_c])
        end

        # Exercise functionality and check results
        for basis in bases_to_test
            local unit_cell

            # Construct UnitCell object
            if isnothing(centering)
                unit_cell = UnitCell(basis[1], basis[2], basis[3])
            else
                unit_cell = UnitCell(basis[1], basis[2], basis[3]; centering=centering)
            end

            # Check results
            @test unit_cell isa expected_type
            @test unit_cell ≈ expected_unit_cell
        end
    end
end
