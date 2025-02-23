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
math_tests.jl defines tests for math.jl
"""
# --- Imports

# Standard library
using LinearAlgebra: dot, qr, I
using Test

# External packages
using Combinatorics: permutations

# Xtallography package
using Xtallography

# --- Tests

# ------ basic mathematical functions

@testset "asin_()" begin
    # --- Preparations

    delta = eps(1.0)

    # --- Tests

    # ------ x in [-1, 1], x far from endpoints

    @test asin_(0.5) == asin(0.5)

    # ------ x near 1

    # x < 1
    @test asin_(1 - delta) == asin(1 - delta)

    # x = 1
    @test asin_(1) == asin(1)

    # x > 1
    @test asin_(1 + delta) == π / 2

    # x > 1, within atol
    @test asin_(1 + 0.5; atol=1) == π / 2

    # ------ x near -1

    # x > -1
    @test asin_(-1 + delta) == asin(-1 + delta)

    # x = -1
    @test asin_(-1) == asin(-1)

    # x < -1
    @test asin_(-1 - delta) == -π / 2

    # x < -1, within atol
    @test asin_(-1 - 0.5; atol=1) == -π / 2

    # ------ x not in [-1, 1], x far from endpoints

    # |x| >> 1, x > 0
    local error = nothing
    local error_message = ""
    try
        asin_(2.0)
    catch error
        bt = catch_backtrace()
        error_message = sprint(showerror, error, bt)
    end

    @test error isa DomainError

    expected_error = "DomainError with 2.0:\nasin(x) is not defined for |x| > 1"
    @test startswith(error_message, expected_error)

    # |x| >> 1, x < 0
    local error = nothing
    local error_message = ""
    try
        asin_(-3.0)
    catch error
        bt = catch_backtrace()
        error_message = sprint(showerror, error, bt)
    end

    @test error isa DomainError

    expected_error = "DomainError with -3.0:\nasin(x) is not defined for |x| > 1"
    @test startswith(error_message, expected_error)
end

@testset "acos_()" begin
    # --- Preparations

    delta = eps(1.0)

    # --- Tests

    # ------ x in [-1, 1], x far from endpoints

    @test acos_(0.5) == acos(0.5)

    # ------ x near 1

    # x < 1
    @test acos_(1 - delta) == acos(1 - delta)

    # x = 1
    @test acos_(1) == acos(1)

    # x > 1
    @test acos_(1 + delta) == 0

    # x > 1, within atol
    @test acos_(1 + 0.5; atol=1) == 0

    # ------ x near -1

    # x < -1
    @test acos_(-1 - delta) == π

    # x = -1
    @test acos_(-1) == acos(-1)

    # x < -1
    @test acos_(-1 - delta) == π

    # x < -1, within atol
    @test acos_(-1 - 0.5; atol=1) == π

    # ------ x not in [-1, 1], x far from endpoints

    # |x| >> 1, x > 0
    local error = nothing
    local error_message = ""
    try
        acos_(2.0)
    catch error
        bt = catch_backtrace()
        error_message = sprint(showerror, error, bt)
    end

    @test error isa DomainError

    expected_error = "DomainError with 2.0:\nacos(x) not defined for |x| > 1"
    @test startswith(error_message, expected_error)

    # |x| >> 1, x < 0
    local error = nothing
    local error_message = ""
    try
        acos_(-3.0)
    catch error
        bt = catch_backtrace()
        error_message = sprint(showerror, error, bt)
    end

    @test error isa DomainError

    expected_error = "DomainError with -3.0:\nacos(x) not defined for |x| > 1"
    @test startswith(error_message, expected_error)
end

# ------ linear algebra

@testset "is_basis()" begin
    # --- Tests

    # ------ v1, v2, v3 are linearly independent

    # Case 1
    v1 = [1, 0, 0]
    v2 = [0, 1, 0]
    v3 = [0, 0, 1]
    @test is_basis(v1, v2, v3)

    # Case 2
    v1 = [1, 2, 3]
    v2 = [4, 5, 6]
    v3 = [7, 8, 10]
    @test is_basis(v1, v2, v3)

    # Case 3: norms of vectors is small
    v1 = 1e-50 * [1, 2, 3]
    v2 = 1e-50 * [4, 5, 6]
    v3 = 1e-50 * [7, 8, 10]
    @test is_basis(v1, v2, v3)

    # ------ v1, v2, v3 are linearly dependent

    # Case 1
    v1 = [1, 0, 0]
    v2 = [0, 1, 0]
    v3 = [1, 1, 0]
    @test !is_basis(v1, v2, v3)

    # Case 2
    v1 = [1, 2, 3]
    v2 = [4, 5, 6]
    v3 = [7, 8, 9]
    @test !is_basis(v1, v2, v3)

    # Case 3: norms of vectors is small
    v1 = 1e-50 * [1, 2, 3]
    v2 = 1e-50 * [4, 5, 6]
    v3 = 1e-50 * [7, 8, 9]
    @test !is_basis(v1, v2, v3)
end

# --- geometry

@testset "volume(::Vector,::Vector,::Vector)" begin
    # --- Tests

    # Case #1
    basis_a = [1, 0, 0]
    basis_b = [0, 2, 0]
    basis_c = [0, 0, 3]

    V = volume(basis_a, basis_b, basis_c)
    @test V ≈ 6

    # Case #2
    basis_a = [1, 0, 0]
    basis_b = [1, 1, 0]
    basis_c = [1, 0, 1]

    V = volume(basis_a, basis_b, basis_c)
    @test V ≈ 1
end

@testset "surface_area(::Vector,::Vector,::Vector)" begin
    # --- Tests

    # Case #1
    basis_a = [1, 0, 0]
    basis_b = [0, 1, 0]
    basis_c = [0, 0, 1]

    S = surface_area(basis_a, basis_b, basis_c)
    @test S ≈ 6

    # Case #2
    basis_a = [1, 0, 0]
    basis_b = [1, 1, 0]
    basis_c = [1, 0, 1]

    S = surface_area(basis_a, basis_b, basis_c)
    @test S ≈ 4 * (1 + sqrt(3) / 2)
end
