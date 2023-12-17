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
Types and functions that support lattice computations
"""
# --- Exports

# Types
export LatticeSystem, Centering
export LatticeConstants
export UnitCell

# Functions
export isapprox, lattice_system, standardize
export isempty, resize
export is_bravais_lattice
export basis, volume, surface_area
export iucr_conventional_cell, reduced_cell
export is_supercell, is_equivalent_unit_cell

# --- Types

"""
    LatticeSystem

Supertype for the seven lattice systems in 3D

Subtypes
========
[`Triclinic`](@ref), [`Monoclinic`](@ref), [`Orthorhombic`](@ref), [`Tetragonal`](@ref),
[`Rhombohedral`](@ref), [`Hexagonal`](@ref), [`Cubic`](@ref)
"""
abstract type LatticeSystem end

""" 
    Centering
    
Enumerated type representing the four lattice centerings in 3D
"""
@enum Centering begin
    PRIMITIVE
    BODY
    FACE
    BASE
end

"""
    LatticeConstants

Supertype for lattice constants for the seven lattice systems in 3D

Subtypes
========
[`TriclinicLatticeConstants`](@ref), [`MonoclinicLatticeConstants`](@ref),
[`OrthorhombicLatticeConstants`](@ref), [`TetragonalLatticeConstants`](@ref),
[`RhombohedralLatticeConstants`](@ref), [`HexagonalLatticeConstants`](@ref),
[`CubicLatticeConstants`](@ref)
"""
abstract type LatticeConstants end

using AngleBetweenVectors: angle
using LinearAlgebra: norm, cross, dot
"""
    LatticeConstants(
        basis_a::Vector{<:Real},
        basis_b::Vector{<:Real},
        basis_c::Vector{<:Real};
        identify_lattice_system=true,
        centering=PRIMITIVE
    ) -> LatticeConstants

Construct a LatticeConstants object from a set of basis vectors.

Keyword Arguments
=================
- `identify_lattice_system`: if set to `true`, return a LatticeConstants object for the
  highest symmetry lattice system that is consistent with the basis. Otherwise, return a
  TriclinicLatticeConstants object.

- `centering`: centering of unit cell

Examples
========
TODO
"""
function LatticeConstants(
    basis_a::Vector{<:Real},
    basis_b::Vector{<:Real},
    basis_c::Vector{<:Real};
    identify_lattice_system=true,
    centering=PRIMITIVE,
)
    # Convert basis vectors to Vector{AbstractFloat}
    basis_a = convert(Vector{Float64}, basis_a)
    basis_b = convert(Vector{Float64}, basis_b)
    basis_c = convert(Vector{Float64}, basis_c)

    # Compute lattice constants
    a = norm(basis_a)
    b = norm(basis_b)
    c = norm(basis_c)
    α = angle(basis_b, basis_c)
    β = angle(basis_c, basis_a)
    γ = angle(basis_a, basis_b)

    if identify_lattice_system
        # Construct LatticeConstants object for the highest symmetry lattice system that
        # is consistent with the basis

        if α ≈ β ≈ γ
            # --- Case: orthorhombic, tetragonal, cubic, or rhombohedral

            if α ≈ β ≈ γ ≈ π / 2
                # --- Case: orthorhombic, tetragonal, or cubic

                if a ≈ b ≈ c
                    # Case: cubic
                    return CubicLatticeConstants(a)

                elseif a ≈ b
                    # Case: tetragonal
                    return TetragonalLatticeConstants(a, c)

                elseif b ≈ c
                    # Case: tetragonal
                    return TetragonalLatticeConstants(b, a)

                elseif c ≈ a
                    # Case: tetragonal
                    return TetragonalLatticeConstants(c, b)

                else
                    # Case: orthorhombic
                    lattice_constants, _ = standardize(
                        OrthorhombicLatticeConstants(a, b, c), centering
                    )
                    return lattice_constants
                end
            else
                if a ≈ b ≈ c
                    # Case: rhombohedral
                    return RhombohedralLatticeConstants(a, α)
                end
            end

        else
            # --- Case: hexagonal, monoclinic, or triclinic

            if a ≈ b && α ≈ β ≈ π / 2 && (γ ≈ 2π / 3 || γ ≈ π / 3)
                # Case: hexagonal
                return HexagonalLatticeConstants(a, c)

            elseif b ≈ c && β ≈ γ ≈ π / 2 && (α ≈ 2π / 3 || α ≈ π / 3)
                # Case: hexagonal
                return HexagonalLatticeConstants(b, a)

            elseif c ≈ a && γ ≈ α ≈ π / 2 && (β ≈ 2π / 3 || β ≈ π / 3)
                # Case: hexagonal
                return HexagonalLatticeConstants(c, b)

            elseif α ≈ γ ≈ π / 2
                # Case: monoclinic
                lattice_constants, _ = standardize(
                    MonoclinicLatticeConstants(a, b, c, β), centering
                )
                return lattice_constants

            elseif β ≈ α ≈ π / 2
                # Case: monoclinic
                lattice_constants, _ = standardize(
                    MonoclinicLatticeConstants(b, c, a, γ), centering
                )
                return lattice_constants

            elseif γ ≈ β ≈ π / 2
                # Case: monoclinic
                lattice_constants, _ = standardize(
                    MonoclinicLatticeConstants(c, a, b, α), centering
                )
                return lattice_constants
            end
        end
    end

    # Case: triclinic
    return standardize(TriclinicLatticeConstants(a, b, c, α, β, γ))
end

"""
    UnitCell

Unit cell for a lattice

Fields
======
* `lattice_constants`: lattice constants of unit cell
* `centering`: centering of unit cell
"""
struct UnitCell
    # Fields
    lattice_constants::LatticeConstants
    centering::Centering

    # TODO: add standardization of lattice constants?
end

# --- Functions/Methods

# ------ LatticeConstants functions

import Base.isapprox

# Default isapprox() implementation to allow comparison between lattice constants of
# different lattice systems.
function isapprox(
    x::LatticeConstants, y::LatticeConstants; atol::Real=0, rtol::Real=atol > 0 ? 0 : √eps()
)
    return false
end

"""
    lattice_system(lattice_constants::LatticesConstants) -> LatticeSystem

Return the lattice system for a set of `lattice_constants`.

Return values
=============
- lattice system
"""
function lattice_system end

"""
    standardize(unit_cell::UnitCell) -> LatticeConstants

Standardize the lattice constants and centering for `unit_cell`.

!!! note

    This function _only_ enforces lattice constant constraints. It _does not modify_ the
    Bravais lattice type. To find an equivalent Bravais lattice with higher symmetry (if
    one exists), use `iucr_conventional_cell()`.

!!! note

    Lattice constant standardizations are based on the conventions provided in the Table
    3.1.4.1. of the International Tables for Crystallography (2016). For triclinic
    lattices, the lattice constants are standardized using the following conventions:

    - `a` ≤ `b` ≤ `c`

    - all three angles are acute (Type I cell) or all three angles are non-acute (Type II
      cell)

    - angles sorted in increasing order when edge lengths are equal
        - `α` ≤ `β` when `a` = `b`
        - `β` ≤ `γ` when `b` = `c`
        - `α` ≤ `β` ≤ `γ` when `a` = `b` = `c`

!!! note

    Except for monoclinic lattices, `standardize()` does not modify the unit cell.

    For monoclinic lattices, the unit cell may be potentially modified in two ways:

    - the 2D unit cell in the plane normal to the __b__-axis may be reduced in order to
      satisfy the IUCr conventions for `a`, `c`, and `β`;

    - base-centered unit cells are transformed to equivalent body-centered unit cells.

Return values
=============
- unit cell with standardized lattice constants and centering

Examples
========
TODO
"""
function standardize(unit_cell::UnitCell)
    return UnitCell(standardize(unit_cell.lattice_constants, unit_cell.centering)...)
end

"""
    standardize(
        lattice_constants::LatticeConstants, centering::Centering
    ) -> (LatticeConstants, Centering)

Standardize the lattice constants and centering for the unit cell defined by
`lattice_constants` and `centering`.

Return values
=============
- standardized lattice constants and centering

Examples
========
TODO
"""
function standardize(lattice_constants::LatticeConstants, centering::Centering)
    # --- Check arguments

    standardize_arg_checks(lattice_constants, centering)

    # --- By default, no standardization is performed

    return lattice_constants, centering
end

"""
    standardize(lattice_constants::LatticeConstants) -> LatticeConstants

Standardize the lattice constants the primitive unit cell defined by `lattice_constants`.

Return values
=============
- standardized lattice constants for primitive unit cell

Examples
========
TODO
"""
function standardize(lattice_constants::LatticeConstants)
    # --- Check arguments

    standardize_arg_checks(lattice_constants, PRIMITIVE)

    # Return standardized lattice constants for primitive unit cell
    standardized_lattice_constants, _ = standardize(lattice_constants, PRIMITIVE)
    return standardized_lattice_constants
end

function standardize_arg_checks(lattice_constants::LatticeConstants, centering::Centering)
    # --- Check arguments

    if !is_bravais_lattice(lattice_system(lattice_constants), centering)
        throw(
            ArgumentError(
                "Invalid Bravais lattice: " *
                "(lattice_system=$(nameof(lattice_system(lattice_constants))), " *
                "centering=$centering)",
            ),
        )
    end
end

# ------ UnitCell functions

function isapprox(x::UnitCell, y::UnitCell; atol::Real=0, rtol::Real=atol > 0 ? 0 : √eps())
    return (
        x.centering == y.centering &&
        isapprox(x.lattice_constants, y.lattice_constants; atol=atol, rtol=rtol)
    )
end

# ------ Unit cell computations

using LinearAlgebra: dot
using Combinatorics: combinations

"""
    is_bravais_lattice(lattice_system::LatticeSystem, centering::Centering) -> Bool

    is_bravais_lattice(lattice_system::Type{<:LatticeSystem}, centering::Centering) -> Bool

    is_bravais_lattice(unit_cell::UnitCell) -> Bool

Determine if the unit cell defined by `unit_cell` or `lattice_system` and `centering` is
a valid Bravais lattice type.

Return values
=============
- `true` if `lattice_system` and `centering` define a valid Bravais lattice type; `false`
  otherwise

Examples
========
TODO
"""
function is_bravais_lattice(lattice_system_::LatticeSystem, centering::Centering)
    if lattice_system_ == Cubic() && centering in (PRIMITIVE, BODY, FACE)
        return true
    elseif lattice_system_ == Tetragonal() && centering in (PRIMITIVE, BODY)
        return true
    elseif lattice_system_ == Orthorhombic() && centering in (PRIMITIVE, BODY, FACE, BASE)
        return true
    elseif lattice_system_ == Hexagonal() && centering == PRIMITIVE
        return true
    elseif lattice_system_ == Rhombohedral() && centering == PRIMITIVE
        return true
    elseif lattice_system_ == Monoclinic() && centering in (PRIMITIVE, BODY, BASE)
        return true
    elseif lattice_system_ == Triclinic() && centering == PRIMITIVE
        return true
    end

    return false
end

function is_bravais_lattice(lattice_system_::Type{<:LatticeSystem}, centering::Centering)
    return is_bravais_lattice(lattice_system_(), centering)
end

function is_bravais_lattice(unit_cell::UnitCell)
    return is_bravais_lattice(
        lattice_system(unit_cell.lattice_constants), unit_cell.centering
    )
end

"""
    basis(unit_cell::UnitCell) -> (Vector{Float64}, Vector{Float64}, Vector{Float64})

    basis(
        lattice_constants::LatticeConstants
    ) -> (Vector{Float64}, Vector{Float64}, Vector{Float64})

Return a set of basis vectors ``\\vec{a}, \\vec{b}, \\vec{c}`` for the unit cell defined by
`unit_cell` or `lattice_constants`.

Return values
=============
- basis vectors ``\\vec{a}``, ``\\vec{b}``, ``\\vec{c}``

Examples
========
TODO
"""
function basis end

function basis(unit_cell::UnitCell)
    return basis(unit_cell.lattice_constants)
end

"""
    volume(unit_cell::UnitCell) -> Float64

    volume(lattice_constants::LatticeConstants) -> Float64

    volume(
        basis_a::Vector{<:Real}, basis_b::Vector{<:Real}, basis_c::Vector{<:Real}
    ) -> Float64

Compute the volume of the unit cell defined by `unit_cell`, `lattice_constants`, or the
basis vectors [`basis_a`, `basis_b`, `basis_c`].

Return values
=============
- volume of the unit cell

Examples
========
TODO
"""
function volume(unit_cell::UnitCell)
    return volume(unit_cell.lattice_constants)
end

"""
    surface_area(unit_cell::UnitCell) -> Float64

    surface_area(lattice_constants::LatticeConstants) -> Float64

    surface_area(
        basis_a::Vector{<:Real}, basis_b::Vector{<:Real}, basis_c::Vector{<:Real}
    ) -> Float64

Compute the surface area of the unit cell defined by `unit_cell`, `lattice_constants`, or
the basis vectors [`basis_a`, `basis_b`, `basis_c`].

Return values
=============
- surface area of the unit cell

Examples
========
```jldoctest
julia> S = surface_area([1, 0, 0], [1, 1, 0], [1, 0, 1]);

julia> S ≈ 4 * (1 + sqrt(3)/2)
true
```
"""
function surface_area(unit_cell::UnitCell)
    return surface_area(unit_cell.lattice_constants)
end

"""
    iucr_conventional_cell(unit_cell::UnitCell) -> UnitCell

Return the IUCr conventional cell that is equivalent to `unit_cell`.

Return values
=============
- IUCr conventional cell for `unit_cell`

Examples
========
TODO
"""
function iucr_conventional_cell(unit_cell::UnitCell)
    # --- Check arguments

    iucr_conventional_cell_arg_checks(unit_cell)

    # --- Compute IUCr conventional cells for lattice systems with limiting cases that
    #     change the Bravais lattice type

    # Get lattice system
    lattice_system_ = lattice_system(unit_cell.lattice_constants)

    # TODO
    if lattice_system_ == Triclinic
        return iucr_conventional_cell(Triclinic(), unit_cell)
    elseif lattice_system_ == Monoclinic
        return iucr_conventional_cell(Monoclinic(), unit_cell)
    elseif lattice_system_ == Orthorhombic
        return iucr_conventional_cell(Orthorhombic(), unit_cell)
    elseif lattice_system_ == Tetragonal
        return iucr_conventional_cell(Tetragonal(), unit_cell)
    elseif lattice_system_ == Rhombohedral
        return iucr_conventional_cell(Rhombohedral(), unit_cell)
    end

    # --- By default, the unit cell is unchanged

    return unit_cell
end

function iucr_conventional_cell_arg_checks(unit_cell::UnitCell)
    # --- Check arguments

    if !is_bravais_lattice(lattice_system(unit_cell.lattice_constants), unit_cell.centering)
        throw(
            ArgumentError(
                "Invalid Bravais lattice: " *
                "(lattice_system=$(nameof(lattice_system(unit_cell.lattice_constants))), " *
                "centering=$(unit_cell.centering))",
            ),
        )
    end
end

"""
    reduced_cell(unit_cell::UnitCell) -> UnitCell

Compute the primitive reduced cell for the lattice defined by `unit_cell`. The
Selling-Delaunay reduction algorithm is used to compute the reduced basis.

Return values
=============
- primitive reduced cell

Examples
========
TODO
"""
function reduced_cell(unit_cell::UnitCell)
    # --- Preparations

    # Get basis of unit cell
    basis_a, basis_b, basis_c = basis(unit_cell.lattice_constants)

    # Construct primitive cell basis for centering
    if unit_cell.centering == BODY
        basis_c = 0.5 * (basis_a + basis_b + basis_c)

    elseif unit_cell.centering == FACE
        basis_a_primitive = 0.5 * (basis_a + basis_b)
        basis_b_primitive = 0.5 * (basis_a - basis_b)
        basis_c_primitive = 0.5 * (basis_a + basis_c)
        basis_a = basis_a_primitive
        basis_b = basis_b_primitive
        basis_c = basis_c_primitive

    elseif unit_cell.centering == BASE
        basis_a = 0.5 * (basis_a + basis_b)
    end

    # Initialize working basis set
    working_basis = [basis_a, basis_b, basis_c, -(basis_a + basis_b + basis_c)]

    # --- Perform Selling-Delaunay reduction

    # ------ Reduce the sum of the squares of the lengths of the working basis set

    # TODO: refactor and test
    while true
        # Compute scalar products
        scalar_products = [
            dot(working_basis[1], working_basis[2]),
            dot(working_basis[1], working_basis[3]),
            dot(working_basis[1], working_basis[4]),
            dot(working_basis[2], working_basis[3]),
            dot(working_basis[2], working_basis[4]),
            dot(working_basis[3], working_basis[4]),
        ]

        # Check if the reduction is complete (when none of the scalar products is positive)
        if all(scalar_products .<= ALMOST_ZERO)
            break
        end

        # Find the largest positive scalar product
        _, idx = findmax(scalar_products)
        if idx == 1
            b_1 = working_basis[1]
            b_2 = working_basis[2]
            b_3 = working_basis[3]
            b_4 = working_basis[4]
        elseif idx == 2
            b_1 = working_basis[1]
            b_2 = working_basis[3]
            b_3 = working_basis[2]
            b_4 = working_basis[4]
        elseif idx == 3
            b_1 = working_basis[1]
            b_2 = working_basis[4]
            b_3 = working_basis[2]
            b_4 = working_basis[3]
        elseif idx == 4
            b_1 = working_basis[2]
            b_2 = working_basis[3]
            b_3 = working_basis[1]
            b_4 = working_basis[4]
        elseif idx == 5
            b_1 = working_basis[2]
            b_2 = working_basis[4]
            b_3 = working_basis[1]
            b_4 = working_basis[3]
        else
            b_1 = working_basis[3]
            b_2 = working_basis[4]
            b_3 = working_basis[1]
            b_4 = working_basis[2]
        end

        # Update working basis set to reduce the sum of the squares of the lengths
        working_basis = [-b_1, b_2, b_1 + b_3, b_1 + b_4]
    end

    # ------ Compute the reduced basis

    # Compute Delaunay set and squared vector lengths
    delaunay_set = [
        (vector=working_basis[1], length_sq=dot(working_basis[1], working_basis[1])),
        (vector=working_basis[2], length_sq=dot(working_basis[2], working_basis[2])),
        (vector=working_basis[3], length_sq=dot(working_basis[3], working_basis[3])),
        (vector=working_basis[4], length_sq=dot(working_basis[4], working_basis[4])),
        (
            vector=working_basis[1] + working_basis[2],
            length_sq=dot(
                working_basis[1] + working_basis[2], working_basis[1] + working_basis[2]
            ),
        ),
        (
            vector=working_basis[1] + working_basis[3],
            length_sq=dot(
                working_basis[1] + working_basis[3], working_basis[1] + working_basis[3]
            ),
        ),
        (
            vector=working_basis[2] + working_basis[3],
            length_sq=dot(
                working_basis[2] + working_basis[3], working_basis[2] + working_basis[3]
            ),
        ),
    ]

    # Sort vectors by squared vector length
    sort!(delaunay_set; by=item -> item.length_sq)

    # TODO: refactor and test
    # Remove vectors that are scalar multiples of shorter vectors
    to_remove = []
    for i in 1:length(delaunay_set)
        basis_1 = delaunay_set[i].vector
        length_basis_1 = sqrt(delaunay_set[i].length_sq)

        for j in (i + 1):length(delaunay_set)
            basis_2 = delaunay_set[j].vector
            length_basis_2 = sqrt(delaunay_set[j].length_sq)

            if dot(basis_1, basis_2) / length_basis_1 / length_basis_2 ≈ 1.0
                push!(to_remove, j)
            end
        end
    end
    delaunay_set = [delaunay_set[i] for i in 1:length(delaunay_set) if !(i in to_remove)]

    # TODO: refactor and test
    # Generate all possible bases that can be formed from the candidate basis vectors
    # and compute the (1) sum of squared lengths of the basis vectors and (2) the
    # surface area of the unit cell
    #candidate_basis_vectors = [ candidate for candidate in delaunay_set]
    basis_choices = [
        (
            vectors=[
                candidate_basis[1].vector,
                candidate_basis[2].vector,
                candidate_basis[3].vector,
            ],
            sum_length_sq=candidate_basis[1].length_sq +
                          candidate_basis[2].length_sq +
                          candidate_basis[3].length_sq,
            surface_area=surface_area(
                candidate_basis[1].vector,
                candidate_basis[2].vector,
                candidate_basis[3].vector,
            ),
        ) for candidate_basis in combinations(delaunay_set, 3)
    ]

    # Eliminate basis choices that are not linearly independent
    basis_choices = [
        candidate_basis for
        candidate_basis in basis_choices if is_basis(candidate_basis.vectors...)
    ]

    # Retain basis choices that have the minimum sum of length squared
    min_sum_length_sq = minimum(item.sum_length_sq for item in basis_choices)
    basis_choices = [
        candidate_basis for candidate_basis in basis_choices if (
            candidate_basis.sum_length_sq < min_sum_length_sq ||
            candidate_basis.sum_length_sq ≈ min_sum_length_sq
        )
    ]

    # Select basis with the maximum surface area
    max_surface_area = maximum(item.surface_area for item in basis_choices)
    reduced_basis = nothing
    for candidate in basis_choices
        if candidate.surface_area == max_surface_area
            reduced_basis = candidate
            break
        end
    end

    # Construct reduced basis (sorted in increasing order of length)
    reduced_basis_a, reduced_basis_b, reduced_basis_c = sort(
        reduced_basis.vectors; by=v -> dot(v, v)
    )

    # TODO: refactor and test
    # Adjust the signs of the basis vectors so that they represent the "homogeneous
    # corner" of the unit cell (i.e., where the three unit cell angles are all acute
    # or all non-acute)
    dot_ab = dot(reduced_basis_a, reduced_basis_b)
    dot_bc = dot(reduced_basis_b, reduced_basis_c)
    dot_ca = dot(reduced_basis_c, reduced_basis_a)

    if dot_ab * dot_bc * dot_ca > 0
        # Unit cell is of Type I, so there is a choice of sign for the basis such that
        # all three unit cell angles are acute

        if dot_ab > 0
            if dot_bc <= 0
                reduced_basis_c *= -1
            end
        elseif dot_bc > 0
            if dot_ca <= 0
                reduced_basis_a *= -1
            end
        else  # dot_ca > 0
            if dot_ab <= 0
                reduced_basis_b *= -1
            end
        end

    else
        # Unit cell is of Type II, so there is a choice of sign for the basis such that
        # all three unit cell angles are non-acute

        if dot_ab <= 0
            if dot_bc > 0
                reduced_basis_c *= -1
            end
        elseif dot_bc <= 0
            if dot_ca > 0
                reduced_basis_a *= -1
            end
        else  # dot_ca <= 0
            if dot_ab > 0
                reduced_basis_b *= -1
            end
        end
    end

    # TODO: refactor and test
    # When there are multiple orderings of the reduced basis vectors by length, reorder
    # the reduced basis vectors so that the angles affected by basis ordering are in
    # increasing
    length_sq_a = dot(reduced_basis_a, reduced_basis_a)
    length_sq_b = dot(reduced_basis_b, reduced_basis_b)
    length_sq_c = dot(reduced_basis_c, reduced_basis_c)

    if length_sq_a ≈ length_sq_b ≈ length_sq_c
        # Use bubble sort to reorder the basis vectors

        if dot(reduced_basis_b, reduced_basis_c) < dot(reduced_basis_c, reduced_basis_a)
            tmp = reduced_basis_a
            reduced_basis_a = reduced_basis_b
            reduced_basis_b = tmp
        end

        if dot(reduced_basis_c, reduced_basis_a) < dot(reduced_basis_a, reduced_basis_b)
            tmp = reduced_basis_b
            reduced_basis_b = reduced_basis_c
            reduced_basis_c = tmp
        end

        if dot(reduced_basis_b, reduced_basis_c) < dot(reduced_basis_c, reduced_basis_a)
            tmp = reduced_basis_a
            reduced_basis_a = reduced_basis_b
            reduced_basis_b = tmp
        end

    elseif length_sq_a ≈ length_sq_b
        if dot(reduced_basis_b, reduced_basis_c) < dot(reduced_basis_c, reduced_basis_a)
            tmp = reduced_basis_a
            reduced_basis_a = reduced_basis_b
            reduced_basis_b = tmp
        end

    elseif length_sq_b ≈ length_sq_c
        if dot(reduced_basis_c, reduced_basis_a) < dot(reduced_basis_a, reduced_basis_b)
            tmp = reduced_basis_b
            reduced_basis_b = reduced_basis_c
            reduced_basis_c = tmp
        end
    end

    # --- Return standardized primitive unit cell defined by reduced basis

    return UnitCell(
        standardize(LatticeConstants(reduced_basis_a, reduced_basis_b, reduced_basis_c)),
        PRIMITIVE,
    )
end

"""
    is_equivalent_unit_cell(
        unit_cell_test::UnitCell,
        unit_cell_ref::UnitCell;
        tol::Real=1e-3
    ) -> Bool

Check if the unit cell defined by `unit_cell_test` is equivalent to the unit cell defined
by `unit_cell_ref`.

Keyword Arguments
=================
- `tol`: absolute tolerance of the deviation between the reduced unit cells defined by
  `unit_cell_test` and `unit_cell_ref`.

Return values
=============
- `true` if the test unit cell is equivalent to the reference unit cell; `false` otherwise

Examples
========
TODO
"""
function is_equivalent_unit_cell(
    unit_cell_test::UnitCell, unit_cell_ref::UnitCell; tol::Real=1e-3
)
    # --- Check arguments

    if tol <= 0
        throw(ArgumentError("`tol` must be positive"))
    end

    # --- Compare reduced cells

    return isapprox(reduced_cell(unit_cell_test), reduced_cell(unit_cell_ref); rtol=tol)
end

"""
    is_equivalent_unit_cell(
        lattice_constants_test::LatticeConstants,
        lattice_constants_ref::LatticeConstants;
        tol::Real=1e-3
    ) -> Bool

Check if the primitive unit cell defined by `lattice_constants_test` is equivalent to the
unit cell defined by `lattice_constants_ref`.

Keyword Arguments
=================
- `tol`: absolute tolerance of the deviation between the reduced unit cells defined by
  `lattice_constants_test` and `lattice_constants_ref`.

Return values
=============
- `true` if the test unit cell is equivalent to the reference unit cell; `false` otherwise

Examples
========
TODO
"""
function is_equivalent_unit_cell(
    lattice_constants_test::LatticeConstants,
    lattice_constants_ref::LatticeConstants;
    tol::Real=1e-3,
)
    # Compare primitive unit cells
    return is_equivalent_unit_cell(
        UnitCell(lattice_constants_test, PRIMITIVE),
        UnitCell(lattice_constants_ref, PRIMITIVE);
        tol=tol,
    )
end

"""
    is_supercell(
        lattice_constants_test::LatticeConstants,
        lattice_constants_ref::LatticeConstants;
        tol::Real=1e-3,
        max_index::Integer=3
    ) -> Bool

Check if the unit cell defined by `lattice_constants_test` is a supercell of the unit cell
defined by `lattice_constants_ref`.

Keyword Arguments
=================
- `tol`: absolute tolerance of the deviation between test unit cell and supercells of the
  reference unit cell

- `max_index`: maximum multiple of the basis vectors of the reference unit cell to include
  when checking whether the test unit cell is a supercell of the reference unit cell

  !!! note

      `max_index` is ignored for Cubic lattice systems.

Return values
=============
- `true` if the test unit cell is a supercell of the reference unit cell; `false` otherwise

Examples
========
TODO
"""
function is_supercell end
