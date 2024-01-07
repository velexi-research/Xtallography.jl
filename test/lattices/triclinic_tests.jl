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
Tests for methods in lattice/triclinic.jl (except for cell standardization methods)
"""
# --- Imports

# Standard library
using Test
using LinearAlgebra: det, norm, cross

# XtallographyUtils package
using XtallographyUtils

# --- Tests

# ------ Types

@testset "TriclinicLatticeConstants constructor: valid arguments" begin
    # --- Preparations

    a = 1
    b = 2
    c = 3
    α = π / 4
    β = π / 2
    γ = 3 * π / 4

    # --- Exercise functionality and check results

    # ------ basic test

    lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)

    @test lattice_constants.a == a
    @test lattice_constants.b == b
    @test lattice_constants.c == c
    @test lattice_constants.α == α
    @test lattice_constants.β == β
    @test lattice_constants.γ == γ

    # ------ edge cases for angles

    # α = 0
    lattice_constants = TriclinicLatticeConstants(a, b, c, 0, β, γ)

    @test lattice_constants.a == a
    @test lattice_constants.b == b
    @test lattice_constants.c == c
    @test lattice_constants.α == 0
    @test lattice_constants.β == β
    @test lattice_constants.γ == γ

    # α = π
    lattice_constants = TriclinicLatticeConstants(a, b, c, π, β, γ)

    @test lattice_constants.a == a
    @test lattice_constants.b == b
    @test lattice_constants.c == c
    @test lattice_constants.α ≈ π
    @test lattice_constants.β == β
    @test lattice_constants.γ == γ

    # β = 0
    lattice_constants = TriclinicLatticeConstants(a, b, c, α, 0, γ)

    @test lattice_constants.a == a
    @test lattice_constants.b == b
    @test lattice_constants.c == c
    @test lattice_constants.α == α
    @test lattice_constants.β == 0
    @test lattice_constants.γ == γ

    # β = π
    lattice_constants = TriclinicLatticeConstants(a, b, c, α, π, γ)

    @test lattice_constants.a == a
    @test lattice_constants.b == b
    @test lattice_constants.c == c
    @test lattice_constants.α == α
    @test lattice_constants.β ≈ π
    @test lattice_constants.γ == γ

    # γ = 0
    lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, 0)

    @test lattice_constants.a == a
    @test lattice_constants.b == b
    @test lattice_constants.c == c
    @test lattice_constants.α == α
    @test lattice_constants.β == β
    @test lattice_constants.γ == 0

    # γ = π
    lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, π)

    @test lattice_constants.a == a
    @test lattice_constants.b == b
    @test lattice_constants.c == c
    @test lattice_constants.α == α
    @test lattice_constants.β == β
    @test lattice_constants.γ ≈ π
end

@testset "TriclinicLatticeConstants constructor: invalid arguments" begin
    # --- Preparations

    # Valid arguments
    a = 1
    b = 2
    c = 3
    α = π / 4
    β = π / 2
    γ = 3 * π / 4

    # --- Tests

    # ------ a

    # a = 0
    local error = nothing
    local error_message = ""
    try
        lattice_constants = TriclinicLatticeConstants(0, b, c, α, β, γ)
    catch error
        bt = catch_backtrace()
        error_message = sprint(showerror, error, bt)
    end

    @test error isa ArgumentError

    expected_error = "ArgumentError: `a` must be positive"
    @test startswith(error_message, expected_error)

    # a < 0
    local error = nothing
    local error_message = ""
    try
        lattice_constants = TriclinicLatticeConstants(-1.0, b, c, α, β, γ)
    catch error
        bt = catch_backtrace()
        error_message = sprint(showerror, error, bt)
    end

    @test error isa ArgumentError

    expected_error = "ArgumentError: `a` must be positive"
    @test startswith(error_message, expected_error)

    # ------ b

    # b = 0
    local error = nothing
    local error_message = ""
    try
        lattice_constants = TriclinicLatticeConstants(a, 0, c, α, β, γ)
    catch error
        bt = catch_backtrace()
        error_message = sprint(showerror, error, bt)
    end

    @test error isa ArgumentError

    expected_error = "ArgumentError: `b` must be positive"
    @test startswith(error_message, expected_error)

    # b < 0
    local error = nothing
    local error_message = ""
    try
        lattice_constants = TriclinicLatticeConstants(a, -1.0, c, α, β, γ)
    catch error
        bt = catch_backtrace()
        error_message = sprint(showerror, error, bt)
    end

    @test error isa ArgumentError

    expected_error = "ArgumentError: `b` must be positive"
    @test startswith(error_message, expected_error)

    # ------ c

    # c = 0
    local error = nothing
    local error_message = ""
    try
        lattice_constants = TriclinicLatticeConstants(a, b, 0, α, β, γ)
    catch error
        bt = catch_backtrace()
        error_message = sprint(showerror, error, bt)
    end

    @test error isa ArgumentError

    expected_error = "ArgumentError: `c` must be positive"
    @test startswith(error_message, expected_error)

    # c < 0
    local error = nothing
    local error_message = ""
    try
        lattice_constants = TriclinicLatticeConstants(a, b, -1.0, α, β, γ)
    catch error
        bt = catch_backtrace()
        error_message = sprint(showerror, error, bt)
    end

    @test error isa ArgumentError

    expected_error = "ArgumentError: `c` must be positive"
    @test startswith(error_message, expected_error)

    # ------ α

    # α < 0
    local error = nothing
    local error_message = ""
    try
        lattice_constants = TriclinicLatticeConstants(a, b, c, -1.0, β, γ)
    catch error
        bt = catch_backtrace()
        error_message = sprint(showerror, error, bt)
    end

    @test error isa ArgumentError

    expected_error = "ArgumentError: `α` must lie in the interval [0, π]"
    @test startswith(error_message, expected_error)

    # α > π
    local error = nothing
    local error_message = ""
    try
        lattice_constants = TriclinicLatticeConstants(a, b, c, 200, β, γ)
    catch error
        bt = catch_backtrace()
        error_message = sprint(showerror, error, bt)
    end

    @test error isa ArgumentError

    expected_error = "ArgumentError: `α` must lie in the interval [0, π]"
    @test startswith(error_message, expected_error)

    # ------ β

    # β < 0
    local error = nothing
    local error_message = ""
    try
        lattice_constants = TriclinicLatticeConstants(a, b, c, α, -1.0, γ)
    catch error
        bt = catch_backtrace()
        error_message = sprint(showerror, error, bt)
    end

    @test error isa ArgumentError

    expected_error = "ArgumentError: `β` must lie in the interval [0, π]"
    @test startswith(error_message, expected_error)

    # β > π
    local error = nothing
    local error_message = ""
    try
        lattice_constants = TriclinicLatticeConstants(a, b, c, α, 200, γ)
    catch error
        bt = catch_backtrace()
        error_message = sprint(showerror, error, bt)
    end

    @test error isa ArgumentError

    expected_error = "ArgumentError: `β` must lie in the interval [0, π]"
    @test startswith(error_message, expected_error)

    # ------ γ

    # γ < 0
    local error = nothing
    local error_message = ""
    try
        lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, -1.0)
    catch error
        bt = catch_backtrace()
        error_message = sprint(showerror, error, bt)
    end

    @test error isa ArgumentError

    expected_error = "ArgumentError: `γ` must lie in the interval [0, π]"
    @test startswith(error_message, expected_error)

    # γ > π
    local error = nothing
    local error_message = ""
    try
        lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, 200)
    catch error
        bt = catch_backtrace()
        error_message = sprint(showerror, error, bt)
    end

    @test error isa ArgumentError

    expected_error = "ArgumentError: `γ` must lie in the interval [0, π]"
    @test startswith(error_message, expected_error)
end

@testset "TriclinicLatticeConstantDeltas constructor" begin
    # --- Tests

    Δa = 1
    Δb = 3
    Δc = 5
    Δα = π / 7
    Δβ = 2π / 7
    Δγ = 3π / 7
    lattice_constants_Δ = TriclinicLatticeConstantDeltas(Δa, Δb, Δc, Δα, Δβ, Δγ)

    @test lattice_constants_Δ.Δa == Δa
    @test lattice_constants_Δ.Δb == Δb
    @test lattice_constants_Δ.Δc == Δc
    @test lattice_constants_Δ.Δα == Δα
    @test lattice_constants_Δ.Δβ == Δβ
    @test lattice_constants_Δ.Δγ == Δγ
end

# ------ LatticeConstants functions

@testset "isapprox(::LatticeConstants)" begin
    # --- Preparations

    x = TriclinicLatticeConstants(1.0, 2.0, 3.0, π / 5, π / 4, 2π / 5)
    y = TriclinicLatticeConstants(1.5, 2.5, 3.5, π / 5 + 0.5, π / 4 + 0.5, 2π / 5 + 0.5)

    # --- Exercise functionality and check results

    # x ≈ (x + delta)
    @test x ≈ TriclinicLatticeConstants(1.0 + 1e-9, 2.0, 3.0, π / 5, π / 4, 2π / 5)
    @test x ≈ TriclinicLatticeConstants(1.0, 2.0 + 1e-9, 3.0, π / 5, π / 4, 2π / 5)
    @test x ≈ TriclinicLatticeConstants(1.0, 2.0, 3.0 - 1e-9, π / 5, π / 4, 2π / 5)
    @test x ≈ TriclinicLatticeConstants(1.0, 2.0, 3.0, π / 5 - 1e-9, π / 4, 2π / 5)
    @test x ≈ TriclinicLatticeConstants(1.0, 2.0, 3.0, π / 5, π / 4 + 1e-9, 2π / 5)
    @test x ≈ TriclinicLatticeConstants(1.0, 2.0, 3.0, π / 5, π / 4, 2π / 5 - 1e-9)

    # x !≈ y
    @test !(x ≈ y)

    # x ≈ y: atol = 1
    @test isapprox(x, y; atol=1)

    # x ≈ y: rtol = 1
    @test isapprox(x, y; rtol=1)

    # x ≈ y: atol = 0.01, rtol = 1
    @test isapprox(x, y; atol=0.01, rtol=1)

    # x ≈ y: atol = 1, rtol = 0.01
    @test isapprox(x, y; atol=1, rtol=0.01)

    # x !≈ y: atol = 0.01, rtol = 0.01
    @test !isapprox(x, y; atol=0.01, rtol=0.01)
end

@testset "lattice_system()" begin
    # --- Tests

    lattice_constants = TriclinicLatticeConstants(1, 2, 3, π / 5, 2π / 5, 3π / 5)
    @test lattice_system(lattice_constants) === Triclinic()
end

@testset "standardize(): Type I cell" begin
    # --- Tests

    # ------ lattice constants already in standard form

    a = 1.0
    b = 5.0
    c = 10.0
    α = 1π / 10
    β = 2π / 10
    γ = 3π / 10
    lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)
    @test is_triclinic_type_I_cell(lattice_constants)

    # Exercise functionality and check results
    standardized_lattice_constants, standardized_centering = standardize(
        lattice_constants, primitive
    )

    expected_lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)
    @test standardized_lattice_constants ≈ expected_lattice_constants

    @test standardized_centering === primitive

    # ------ edge lengths not sorted

    a = 1.0
    b = 10.0
    c = 5.0
    α = 1π / 10
    β = 2π / 10
    γ = 3π / 10
    lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)
    @test is_triclinic_type_I_cell(lattice_constants)

    # Exercise functionality and check results
    standardized_lattice_constants, standardized_centering = standardize(
        lattice_constants, primitive
    )

    expected_lattice_constants = TriclinicLatticeConstants(a, c, b, α, γ, β)
    @test standardized_lattice_constants ≈ expected_lattice_constants

    @test standardized_centering === primitive

    # ------ origin not at "homogeneous corner"

    a = 1.0
    b = 5.0
    c = 10.0
    α = 1π / 10
    β = 8π / 10
    γ = 7π / 10
    lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)
    @test is_triclinic_type_I_cell(lattice_constants)

    # Exercise functionality and check results
    standardized_lattice_constants, standardized_centering = standardize(
        lattice_constants, primitive
    )

    expected_lattice_constants = TriclinicLatticeConstants(a, b, c, α, π - β, π - γ)
    @test standardized_lattice_constants ≈ expected_lattice_constants

    @test standardized_centering === primitive

    # ------ edge lengths not sorted and origin not at "homogeneous corner"

    a = 10.0
    b = 1.0
    c = 5.0
    α = 1π / 10
    β = 8π / 10
    γ = 7π / 10
    lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)
    @test is_triclinic_type_I_cell(lattice_constants)

    # Exercise functionality and check results
    standardized_lattice_constants, standardized_centering = standardize(
        lattice_constants, primitive
    )

    expected_lattice_constants = TriclinicLatticeConstants(b, c, a, π - β, π - γ, α)
    @test standardized_lattice_constants ≈ expected_lattice_constants

    @test standardized_centering === primitive

    # ------ a ≈ b ≈ c, only need to swap α and β

    a = b = c = 1.0
    α = 2π / 10
    β = 1π / 10
    γ = 3π / 10
    lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)
    @test is_triclinic_type_I_cell(lattice_constants)

    # Exercise functionality and check results
    standardized_lattice_constants, standardized_centering = standardize(
        lattice_constants, primitive
    )

    expected_lattice_constants = TriclinicLatticeConstants(a, b, c, β, α, γ)
    @test standardized_lattice_constants ≈ expected_lattice_constants

    @test standardized_centering === primitive

    # ------ a ≈ b ≈ c, only need to swap β and γ

    a = b = c = 1.0
    α = 1π / 10
    β = 3π / 10
    γ = 2π / 10
    lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)
    @test is_triclinic_type_I_cell(lattice_constants)

    # Exercise functionality and check results
    standardized_lattice_constants, standardized_centering = standardize(
        lattice_constants, primitive
    )

    expected_lattice_constants = TriclinicLatticeConstants(a, b, c, α, γ, β)
    @test standardized_lattice_constants ≈ expected_lattice_constants

    @test standardized_centering === primitive

    # ------ a ≈ b ≈ c, only need to reverse order of angles

    a = b = c = 1.0
    α = 3π / 10
    β = 2π / 10
    γ = 1π / 10
    lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)
    @test is_triclinic_type_I_cell(lattice_constants)

    # Exercise functionality and check results
    standardized_lattice_constants, standardized_centering = standardize(
        lattice_constants, primitive
    )

    expected_lattice_constants = TriclinicLatticeConstants(a, b, c, γ, β, α)
    @test standardized_lattice_constants ≈ expected_lattice_constants

    @test standardized_centering === primitive

    # ------ a ≈ b (after sorting a, b, c)

    a = b = 1.0
    c = 5
    α = 4π / 10
    β = 8π / 10
    γ = 7π / 10
    lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)
    @test is_triclinic_type_I_cell(lattice_constants)

    # Exercise functionality and check results
    standardized_lattice_constants, standardized_centering = standardize(
        lattice_constants, primitive
    )

    expected_lattice_constants = TriclinicLatticeConstants(a, b, c, π - β, α, π - γ)
    @test standardized_lattice_constants ≈ expected_lattice_constants

    @test standardized_centering === primitive

    # ------ b ≈ c (after sorting a, b, c)

    a = b = 5.0
    c = 1
    α = 8π / 10
    β = 7π / 10
    γ = 4π / 10
    lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)
    @test is_triclinic_type_I_cell(lattice_constants)

    # Exercise functionality and check results
    standardized_lattice_constants, standardized_centering = standardize(
        lattice_constants, primitive
    )

    expected_lattice_constants = TriclinicLatticeConstants(c, a, b, γ, π - α, π - β)
    @test standardized_lattice_constants ≈ expected_lattice_constants

    @test standardized_centering === primitive
end

@testset "standardize(): Type II cell" begin
    # --- Tests

    # ------ lattice constants already in standard form

    a = 1.0
    b = 5.0
    c = 10.0
    α = 6π / 10
    β = 7π / 10
    γ = 8π / 10
    lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)
    @test !is_triclinic_type_I_cell(lattice_constants)

    # Exercise functionality and check results
    standardized_lattice_constants, standardized_centering = standardize(
        lattice_constants, primitive
    )

    expected_lattice_constants = lattice_constants
    @test standardized_lattice_constants ≈ expected_lattice_constants

    @test standardized_centering === primitive

    # ------ edge lengths not sorted

    a = 1.0
    b = 10.0
    c = 5.0
    α = 6π / 10
    β = 7π / 10
    γ = 8π / 10
    lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)
    @test !is_triclinic_type_I_cell(lattice_constants)

    # Exercise functionality and check results
    standardized_lattice_constants, standardized_centering = standardize(
        lattice_constants, primitive
    )

    expected_lattice_constants = TriclinicLatticeConstants(a, c, b, α, γ, β)
    @test standardized_lattice_constants ≈ expected_lattice_constants

    @test standardized_centering === primitive

    # ------ origin not at "homogeneous corner"

    a = 1.0
    b = 5.0
    c = 10.0
    α = 4π / 10
    β = 3π / 10
    γ = 8π / 10
    lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)
    @test !is_triclinic_type_I_cell(lattice_constants)

    # Exercise functionality and check results
    standardized_lattice_constants, standardized_centering = standardize(
        lattice_constants, primitive
    )

    expected_lattice_constants = TriclinicLatticeConstants(a, b, c, π - α, π - β, γ)
    @test standardized_lattice_constants ≈ expected_lattice_constants

    @test standardized_centering === primitive

    # ------ edge lengths not sorted and origin not at "homogeneous corner"

    a = 10.0
    b = 1.0
    c = 5.0
    α = 4π / 10
    β = 3π / 10
    γ = 8π / 10
    lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)
    @test !is_triclinic_type_I_cell(lattice_constants)

    # Exercise functionality and check results
    standardized_lattice_constants, standardized_centering = standardize(
        lattice_constants, primitive
    )

    expected_lattice_constants = TriclinicLatticeConstants(b, c, a, π - β, γ, π - α)
    @test standardized_lattice_constants ≈ expected_lattice_constants

    @test standardized_centering === primitive

    # ------ a ≈ b ≈ c, only need to swap α and β

    a = b = c = 1.0
    α = 8π / 10
    β = 7π / 10
    γ = 9π / 10
    lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)
    @test !is_triclinic_type_I_cell(lattice_constants)

    # Exercise functionality and check results
    standardized_lattice_constants, standardized_centering = standardize(
        lattice_constants, primitive
    )

    expected_lattice_constants = TriclinicLatticeConstants(a, b, c, β, α, γ)
    @test standardized_lattice_constants ≈ expected_lattice_constants

    @test standardized_centering === primitive

    # ------ a ≈ b ≈ c, only need to swap β and γ

    a = b = c = 1.0
    α = 7π / 10
    β = 9π / 10
    γ = 8π / 10
    lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)
    @test !is_triclinic_type_I_cell(lattice_constants)

    # Exercise functionality and check results
    standardized_lattice_constants, standardized_centering = standardize(
        lattice_constants, primitive
    )

    expected_lattice_constants = TriclinicLatticeConstants(a, b, c, α, γ, β)
    @test standardized_lattice_constants ≈ expected_lattice_constants

    @test standardized_centering === primitive

    # ------ a ≈ b ≈ c, only need to reverse order of angles

    a = b = c = 1.0
    α = 9π / 10
    β = 8π / 10
    γ = 7π / 10
    lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)
    @test !is_triclinic_type_I_cell(lattice_constants)

    # Exercise functionality and check results
    standardized_lattice_constants, standardized_centering = standardize(
        lattice_constants, primitive
    )

    expected_lattice_constants = TriclinicLatticeConstants(a, b, c, γ, β, α)
    @test standardized_lattice_constants ≈ expected_lattice_constants

    @test standardized_centering === primitive

    # ------ a ≈ b (after sorting a, b, c)

    a = b = 1.0
    c = 5
    α = 7π / 10
    β = 4π / 10
    γ = 4π / 10
    lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)
    @test !is_triclinic_type_I_cell(lattice_constants)

    # Exercise functionality and check results
    standardized_lattice_constants, standardized_centering = standardize(
        lattice_constants, primitive
    )

    expected_lattice_constants = TriclinicLatticeConstants(a, b, c, π - β, α, π - γ)
    @test standardized_lattice_constants ≈ expected_lattice_constants

    @test standardized_centering === primitive

    # ------ b ≈ c (after sorting a, b, c)

    a = b = 5.0
    c = 1
    α = 8π / 10
    β = 4π / 10
    γ = 3π / 10
    lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)
    @test !is_triclinic_type_I_cell(lattice_constants)

    # Exercise functionality and check results
    standardized_lattice_constants, standardized_centering = standardize(
        lattice_constants, primitive
    )

    expected_lattice_constants = TriclinicLatticeConstants(c, a, b, π - γ, π - β, α)
    @test standardized_lattice_constants ≈ expected_lattice_constants

    @test standardized_centering === primitive
end

@testset "standardize(): invalid arguments" begin
    # --- Preparations

    lattice_constants = TriclinicLatticeConstants(1, 2, 3, 2π / 5, 3π / 5, 4π / 5)

    # --- Tests

    for centering in (body_centered, face_centered, base_centered)
        local error = nothing
        local error_message = ""
        try
            standardize(lattice_constants, centering)
        catch error
            bt = catch_backtrace()
            error_message = sprint(showerror, error, bt)
        end

        @test error isa ArgumentError

        expected_error =
            "ArgumentError: " *
            "Invalid Bravais lattice: " *
            "(lattice_system=Triclinic, centering=$(nameof(typeof(centering))))"

        @test startswith(error_message, expected_error)
    end
end

@testset "satisfies_triclinic_angle_constraints()" begin
    # --- Tests

    # valid angles for a triclinic unit cell
    @test satisfies_triclinic_angle_constraints(π / 4, π / 5, π / 6)

    # invalid angles for a triclinic unit cell
    @test !satisfies_triclinic_angle_constraints(3π / 4, 4π / 5, 5π / 6)
end

# ------ LatticeConstantDeltas functions

@testset "isapprox(::LatticeConstantDeltas)" begin
    # --- Preparations

    x = TriclinicLatticeConstantDeltas(1.0, 2.0, 3.0, π / 5, π / 4, 2π / 5)
    y = TriclinicLatticeConstantDeltas(
        1.5, 2.5, 3.5, π / 5 + 0.5, π / 4 + 0.5, 2π / 5 + 0.5
    )

    # --- Exercise functionality and check results

    # x ≈ (x + delta)
    @test x ≈ TriclinicLatticeConstantDeltas(1.0 + 1e-9, 2.0, 3.0, π / 5, π / 4, 2π / 5)
    @test x ≈ TriclinicLatticeConstantDeltas(1.0, 2.0 + 1e-9, 3.0, π / 5, π / 4, 2π / 5)
    @test x ≈ TriclinicLatticeConstantDeltas(1.0, 2.0, 3.0 - 1e-9, π / 5, π / 4, 2π / 5)
    @test x ≈ TriclinicLatticeConstantDeltas(1.0, 2.0, 3.0, π / 5 - 1e-9, π / 4, 2π / 5)
    @test x ≈ TriclinicLatticeConstantDeltas(1.0, 2.0, 3.0, π / 5, π / 4 + 1e-9, 2π / 5)
    @test x ≈ TriclinicLatticeConstantDeltas(1.0, 2.0, 3.0, π / 5, π / 4, 2π / 5 - 1e-9)

    # x !≈ y
    @test !(x ≈ y)

    # x ≈ y: atol = 1
    @test isapprox(x, y; atol=1)

    # x ≈ y: rtol = 1
    @test isapprox(x, y; rtol=1)

    # x ≈ y: atol = 0.01, rtol = 1
    @test isapprox(x, y; atol=0.01, rtol=1)

    # x ≈ y: atol = 1, rtol = 0.01
    @test isapprox(x, y; atol=1, rtol=0.01)

    # x !≈ y: atol = 0.01, rtol = 0.01
    @test !isapprox(x, y; atol=0.01, rtol=0.01)
end

# ------ Unit cell computations

@testset "basis()" begin
    # --- Preparations

    a = 2
    b = 3
    c = 5
    α = π / 3
    β = π / 4
    γ = π / 5
    lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)

    # --- Exercise functionality

    basis_a, basis_b, basis_c = basis(lattice_constants)

    # --- Check results

    @test basis_a ≈ [a, 0, 0]
    @test basis_b ≈ [b * cos(γ), b * sin(γ), 0]

    V = volume(lattice_constants)
    @test basis_c ≈
        [c * cos(β), c / sin(γ) * (cos(α) - cos(β) * cos(γ)), V / sin(γ) / a / b]
end

@testset "volume()" begin
    # --- Preparations

    # Construct basis vectors for unit cell
    a = 2
    b = 3
    c = 5
    α = π / 3
    β = π / 4
    γ = π / 5
    lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)

    basis_a, basis_b, basis_c = basis(lattice_constants)

    # --- Exercise functionality and check results

    @test volume(lattice_constants) ≈ abs(det(hcat(basis_a, basis_b, basis_c)))
end

@testset "surface_area(::LatticeConstants)" begin
    # --- Preparations

    # Construct basis vectors for unit cell
    a = 2
    b = 3
    c = 5
    α = π / 3
    β = π / 4
    γ = π / 5
    lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)

    basis_a, basis_b, basis_c = basis(lattice_constants)

    # --- Exercise functionality and check results

    @test surface_area(lattice_constants) ≈
        2 * norm(cross(basis_a, basis_b)) +
          2 * norm(cross(basis_b, basis_c)) +
          2 * norm(cross(basis_c, basis_a))
end

@testset "reduced_cell()" begin
    # --- Preparations

    a = 5
    b = 8
    c = 10
    α = 2π / 5
    β = 3π / 5
    γ = 4π / 5
    lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)
    basis_a, basis_b, basis_c = basis(lattice_constants)

    # --- Exercise functionality and check results

    # primitive unit cell defined by [basis_a, basis_b, basis_c]
    unit_cell = UnitCell(lattice_constants, primitive)

    expected_reduced_cell = reduced_cell(UnitCell(lattice_constants, primitive))

    reduced_cell_ = reduced_cell(unit_cell)
    @test reduced_cell_.lattice_constants isa TriclinicLatticeConstants
    @test volume(reduced_cell_) ≈ volume(unit_cell)
    @test reduced_cell_ ≈ expected_reduced_cell

    # primitive unit cell defined by linear combination of [basis_a, basis_b, basis_c]
    unit_cell = UnitCell(
        LatticeConstants(basis_a + basis_b + basis_c, basis_b, basis_c), primitive
    )

    expected_reduced_cell = reduced_cell(UnitCell(lattice_constants, primitive))

    reduced_cell_ = reduced_cell(unit_cell)
    @test reduced_cell_.lattice_constants isa TriclinicLatticeConstants
    @test volume(reduced_cell_) ≈ volume(unit_cell)
    @test reduced_cell_ ≈ expected_reduced_cell

    #=
    # TODO: Delaunay reduction does not terminate
    a = 5
    b = 8
    c = 10
    α = π / 5
    β = 2π / 5
    γ = 3π / 5
    lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)
    basis_a, basis_b, basis_c = basis(lattice_constants)
    unit_cell = UnitCell(
        LatticeConstants(basis_a + basis_b + basis_c, basis_b, basis_c), primitive
    )

    expected_reduced_cell = reduced_cell(UnitCell(lattice_constants, primitive))

    reduced_cell_ = reduced_cell(unit_cell)
    =#
end

@testset "is_equivalent_unit_cell(::UnitCell, ::UnitCell)" begin
    # --- Preparations

    a = 2
    b = 3
    c = 5
    α = 3π / 7
    β = 4π / 7
    γ = 5π / 7
    lattice_constants = TriclinicLatticeConstants(a, b, c, α, β, γ)
    basis_a, basis_b, basis_c = basis(lattice_constants)

    # --- Tests

    # identical triclinic unit cells
    unit_cell_ref = UnitCell(lattice_constants, primitive)
    unit_cell_test = UnitCell(
        LatticeConstants(basis_a, basis_b, basis_c; identify_lattice_system=false),
        primitive,
    )
    @test is_equivalent_unit_cell(unit_cell_test, unit_cell_ref)

    # primitive unit cell defined by linear combination of [basis_a, basis_b, basis_c]
    unit_cell_ref = reduced_cell(UnitCell(lattice_constants, primitive))
    unit_cell_test = UnitCell(
        LatticeConstants(basis_a + basis_b + basis_c, basis_b, basis_c), primitive
    )
    @test is_equivalent_unit_cell(unit_cell_test, unit_cell_ref)

    # body-centered unit cell
    body_centered_unit_cell = UnitCell(lattice_constants, body_centered)
    primitive_unit_cell = UnitCell(
        LatticeConstants(
            basis_a,
            basis_b,
            0.5 * (basis_a + basis_b + basis_c);
            identify_lattice_system=false,
        ),
        primitive,
    )
    @test is_equivalent_unit_cell(body_centered_unit_cell, primitive_unit_cell)
end

@testset "is_equivalent_unit_cell(::LatticeConstants, ::LatticeConstants)" begin
    # --- Preparations

    a_ref = 6
    b_ref = 10
    c_ref = 8
    α_ref = π / 3
    β_ref = π / 4
    γ_ref = π / 5
    lattice_constants_ref = TriclinicLatticeConstants(
        a_ref, b_ref, c_ref, α_ref, β_ref, γ_ref
    )
    basis_a, basis_b, basis_c = basis(lattice_constants_ref)

    # --- Exercise functionality and check results

    # equivalent primitive unit cell defined by linear combination of
    # [basis_a, basis_b, basis_c]
    lattice_constants_test = LatticeConstants(basis_a + basis_b + basis_c, basis_b, basis_c)
    @test is_equivalent_unit_cell(lattice_constants_test, lattice_constants_ref)

    # equivalent primitive unit cell defined by a different linear combination of
    # [basis_a, basis_b, basis_c]
    lattice_constants_test = LatticeConstants(
        basis_a + basis_b + basis_c, basis_b + basis_c, basis_c
    )
    @test is_equivalent_unit_cell(lattice_constants_test, lattice_constants_ref)

    # test unit cell and reference unit cell are for different lattice systems
    lattice_constants_test = CubicLatticeConstants(1)
    @test !is_equivalent_unit_cell(lattice_constants_test, lattice_constants_ref)
end

@testset "is_triclinic_type_I_cell()" begin
    # --- Tests

    # ------ Type I

    # no angles equal to π/2
    a = 1
    b = 1
    c = 1
    α = 2π / 10
    β = 3π / 10
    γ = 4π / 10
    @test is_triclinic_type_I_cell(TriclinicLatticeConstants(a, b, c, α, β, γ))

    # α equal to π/2
    a = 1
    b = 1
    c = 1
    α = π / 2
    β = 3π / 10
    γ = 4π / 10
    @test is_triclinic_type_I_cell(TriclinicLatticeConstants(a, b, c, α, β, γ))

    # α equal to π/2 + ϵ (so that cos(α) < 0)
    a = 1
    b = 1
    c = 1
    α = π / 2 + 1e-9
    β = 3π / 10
    γ = 4π / 10
    @test is_triclinic_type_I_cell(TriclinicLatticeConstants(a, b, c, α, β, γ))

    # ------ Type II

    # Type II: no angles equal to π/2
    a = 1
    b = 1
    c = 1
    α = 6π / 10
    β = 7π / 10
    γ = 8π / 10
    @test !is_triclinic_type_I_cell(TriclinicLatticeConstants(a, b, c, α, β, γ))

    # γ equal to π/2
    a = 1
    b = 1
    c = 1
    α = 8π / 10
    β = 7π / 10
    γ = π / 2
    @test is_triclinic_type_I_cell(TriclinicLatticeConstants(a, b, c, α, β, γ))

    # β equal to π/2 + ϵ (so that cos(β) < 0)
    a = 1
    b = 1
    c = 1
    α = 7π / 10
    β = π / 2 + 1e-9
    γ = 8π / 10
    @test is_triclinic_type_I_cell(TriclinicLatticeConstants(a, b, c, α, β, γ))
end
