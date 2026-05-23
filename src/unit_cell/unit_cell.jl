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
Unit cell type and functions
"""
# --- Exports

# Types
export UnitCell

# Functions
export lattice_system
export lattice_constants, symmetry, centering, symmetry_elements
export is_bravais_lattice
export basis, volume, surface_area
export standardize, conventional_cell
export reduced_cell, compute_delaunay_set, prune_delaunay_set, find_reduced_basis
export is_equivalent, is_supercell
export isapprox

# --- Types

"""
    UnitCell{T<:LatticeSystem}

Unit cell for the lattice system `T`

Type Aliases
============
[`TriclinicUnitCell`](@ref), [`MonoclinicUnitCell`](@ref), [`OrthorhombicUnitCell`](@ref),
[`HexagonalUnitCell`](@ref), [`RhombohedralUnitCell`](@ref), [`TetragonalUnitCell`](@ref),
[`CubicUnitCell`](@ref)
"""
struct UnitCell{T<:LatticeSystem}
    # --- Fields

    lattice_constants::NamedTuple
    symmetry::UnitCellSymmetry

    # --- Constructors
    #
    # Notes
    # -----
    # - For efficiency of argument checking, a separate constructor is provided for each
    #   lattice system.

    function UnitCell{T}(
        lattice_constants::NamedTuple, symmetry::UnitCellSymmetry
    ) where {T<:LatticeSystem}
        # Default constructor to allow extension to custom lattice systems

        return new(lattice_constants, symmetry)
    end

    function UnitCell{Triclinic}(lattice_constants::NamedTuple, symmetry::UnitCellSymmetry)

        # --- Check arguments

        if Set(keys(lattice_constants)) != Set((:a, :b, :c, :α, :β, :γ))
            throw(
                ArgumentError(
                    "Invalid lattice_constants argument passed to " *
                    "UnitCell{Triclinic} constructor. " *
                    "Expected keys: (:a, :b, :c, :α, :β, :γ). " *
                    "Provided keys: $(keys(lattice_constants)).",
                ),
            )
        end

        # --- Return new UnitCell

        return new(lattice_constants, symmetry)
    end

    function UnitCell{Monoclinic}(lattice_constants::NamedTuple, symmetry::UnitCellSymmetry)

        # --- Check arguments

        if Set(keys(lattice_constants)) != Set((:a, :b, :c, :β))
            throw(
                ArgumentError(
                    "Invalid lattice_constants argument passed to " *
                    "UnitCell{Monoclinic} constructor. " *
                    "Expected keys: (:a, :b, :c, :β). " *
                    "Provided keys: $(keys(lattice_constants)).",
                ),
            )
        end

        # --- Return new UnitCell

        return new(lattice_constants, symmetry)
    end

    function UnitCell{Orthorhombic}(
        lattice_constants::NamedTuple, symmetry::UnitCellSymmetry
    )

        # --- Check arguments

        if Set(keys(lattice_constants)) != Set((:a, :b, :c))
            throw(
                ArgumentError(
                    "Invalid lattice_constants argument passed to " *
                    "UnitCell{Orthorhombic} constructor. " *
                    "Expected keys: (:a, :b, :c). " *
                    "Provided keys: $(keys(lattice_constants)).",
                ),
            )
        end

        # --- Return new UnitCell

        return new(lattice_constants, symmetry)
    end

    function UnitCell{Tetragonal}(lattice_constants::NamedTuple, symmetry::UnitCellSymmetry)

        # --- Check arguments

        if Set(keys(lattice_constants)) != Set((:a, :c))
            throw(
                ArgumentError(
                    "Invalid lattice_constants argument passed to " *
                    "UnitCell{Tetragonal} constructor. " *
                    "Expected keys: (:a, :c). " *
                    "Provided keys: $(keys(lattice_constants)).",
                ),
            )
        end

        # --- Return new UnitCell

        return new(lattice_constants, symmetry)
    end

    function UnitCell{Rhombohedral}(
        lattice_constants::NamedTuple, symmetry::UnitCellSymmetry
    )

        # --- Check arguments

        if Set(keys(lattice_constants)) != Set((:a, :α))
            throw(
                ArgumentError(
                    "Invalid lattice_constants argument passed to " *
                    "UnitCell{Rhombohedral} constructor. " *
                    "Expected keys: (:a, :α). " *
                    "Provided keys: $(keys(lattice_constants)).",
                ),
            )
        end

        # --- Return new UnitCell

        return new(lattice_constants, symmetry)
    end

    function UnitCell{Hexagonal}(lattice_constants::NamedTuple, symmetry::UnitCellSymmetry)

        # --- Check arguments

        if Set(keys(lattice_constants)) != Set((:a, :c))
            throw(
                ArgumentError(
                    "Invalid lattice_constants argument passed to " *
                    "UnitCell{Hexagonal} constructor. " *
                    "Expected keys: (:a, :c). " *
                    "Provided keys: $(keys(lattice_constants)).",
                ),
            )
        end

        # --- Return new UnitCell

        return new(lattice_constants, symmetry)
    end

    function UnitCell{Cubic}(lattice_constants::NamedTuple, symmetry::UnitCellSymmetry)

        # --- Check arguments

        if Set(keys(lattice_constants)) != Set((:a,))
            throw(
                ArgumentError(
                    "Invalid lattice_constants argument passed to " *
                    "UnitCell{Cubic} constructor. " *
                    "Expected keys: (:a,). " *
                    "Provided keys: $(keys(lattice_constants)).",
                ),
            )
        end

        # --- Return new UnitCell

        return new(lattice_constants, symmetry)
    end
end

# ------ Outer constructors

using AngleBetweenVectors: angle
using LinearAlgebra: norm, cross, dot

"""
    UnitCell{T}(
        lattice_constants::NamedTuple;
        centering::Centering=primitive_centering,
        symmetry_elements::Union{Set,Vector,Nothing}=nothing
    ) -> UnitCell{T}

Construct a UnitCell from a set of lattice constants. The lattice system `T` must be
consistent with the fields present in `lattice_constants`.

Keyword Arguments
=================
- `centering`: centering of unit cell

- `symmetry_elements`: symmetry elements of crystal
"""
@inline function UnitCell{T}(
    lattice_constants::NamedTuple;
    centering::Centering=primitive_centering,
    symmetry_elements::Union{Set,Vector,Nothing}=nothing,
) where {T<:LatticeSystem}
    return UnitCell{T}(
        lattice_constants, UnitCellSymmetry(centering; symmetry_elements=symmetry_elements)
    )
end

"""
    UnitCell(unit_cell::UnitCell{T}) -> UnitCell{T}

Copy constructor. Construct a copy of `unit_cell.`
"""
@inline function UnitCell(unit_cell::UnitCell{T}) where {T<:LatticeSystem}
    return UnitCell{T}(lattice_constants(unit_cell), symmetry(unit_cell))
end

"""
    UnitCell(
        basis_a::Vector{<:Real},
        basis_b::Vector{<:Real},
        basis_c::Vector{<:Real};
        identify_lattice_system::Bool=true,
        centering::Centering=primitive_centering,
    ) -> UnitCell

Construct a UnitCell from a set of basis vectors.

Keyword Arguments
=================
- `identify_lattice_system`: if `true`, return a UnitCell with (1) the highest symmetry
  lattice system that is consistent with the basis vectors and (2) standardized lattice
  constants. Otherwise, return a TriclinicUnitCell with lattice constants computed directly
  from the basis vectors (without lattice constants standardization).

- `centering`: centering of unit cell

Examples
========
```jldoctest
julia> UnitCell([1, 0, 0], [0, 2, 0], [0, 0, 1])
TetragonalUnitCell((a = 1.0, c = 2.0), UnitCellSymmetry(PrimitiveCentering(), Set{SymmetryElement}()))
```
"""
function UnitCell(
    basis_a::Vector{<:Real},
    basis_b::Vector{<:Real},
    basis_c::Vector{<:Real};
    identify_lattice_system::Bool=true,
    centering::Centering=primitive_centering,
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
        # Construct a UnitCell object for the highest symmetry lattice system that is
        # consistent with the basis

        # --- Case: orthorhombic, tetragonal, cubic, or rhombohedral

        if α ≈ β ≈ γ ≈ π / 2
            # --- Case: orthorhombic, tetragonal, or cubic

            if a ≈ b ≈ c
                # Case: cubic
                return CubicUnitCell(a; centering=centering)

            elseif a ≈ b
                # Case: tetragonal
                return TetragonalUnitCell(a, c; centering=centering)

            elseif b ≈ c
                # Case: tetragonal
                return TetragonalUnitCell(b, a; centering=centering)

            elseif c ≈ a
                # Case: tetragonal
                return TetragonalUnitCell(c, b; centering=centering)

            else
                # Case: orthorhombic
                return standardize(OrthorhombicUnitCell(a, b, c; centering=centering))
            end

        else
            # --- Case: hexagonal, monoclinic, rhombohedral, or triclinic

            if a ≈ b &&
                α ≈ β ≈ π / 2 &&
                (γ ≈ 2π / 3 || γ ≈ π / 3) &&
                centering === P_centering
                # Case: hexagonal
                return HexagonalUnitCell(a, c)

            elseif b ≈ c &&
                β ≈ γ ≈ π / 2 &&
                (α ≈ 2π / 3 || α ≈ π / 3) &&
                centering === P_centering
                # Case: hexagonal
                return HexagonalUnitCell(b, a)

            elseif c ≈ a &&
                γ ≈ α ≈ π / 2 &&
                (β ≈ 2π / 3 || β ≈ π / 3) &&
                centering === P_centering
                # Case: hexagonal
                return HexagonalUnitCell(c, b)

            elseif α ≈ γ ≈ π / 2
                # Case: monoclinic
                return standardize(MonoclinicUnitCell(a, b, c, β; centering=centering))

            elseif β ≈ α ≈ π / 2
                # Case: monoclinic
                return standardize(MonoclinicUnitCell(b, c, a, γ; centering=centering))

            elseif γ ≈ β ≈ π / 2
                # Case: monoclinic
                return standardize(MonoclinicUnitCell(c, a, b, α; centering=centering))

            else
                # Case: rhombohedral or triclinic

                # Standardize triclinic lattice constants (needed to correctly identify
                # rhombohedral unit cells)
                a, b, c, α, β, γ = standardize(a, b, c, α, β, γ)

                if α ≈ β ≈ γ && a ≈ b ≈ c
                    # Case: rhombohedral
                    return RhombohedralUnitCell(a, α)
                else
                    # Case: triclinic (with lattice constant standardization)
                    return TriclinicUnitCell(a, b, c, α, β, γ)
                end
            end
        end
    end

    # Return triclinic unit cell (without lattice constant standardization)
    return TriclinicUnitCell(a, b, c, α, β, γ)
end

# --- Functions/Methods

import Base.:(-)
import Base.:(==)
import Base.isapprox
using LinearAlgebra: dot
using Combinatorics: combinations

@inline function lattice_system(unit_cell::UnitCell{T}) where {T<:LatticeSystem}
    return T()
end

"""
    lattice_constants(unit_cell::UnitCell) -> NamedTuple

Return the lattice constants for `unit_cell`.

Return values
=============
lattice constants
"""
@inline function lattice_constants(unit_cell::UnitCell)
    return unit_cell.lattice_constants
end

"""
    symmetry(unit_cell::UnitCell) -> UnitCellSymmetry

Return the symmetry for `unit_cell`.

Return values
=============
symmetry
"""
@inline function symmetry(unit_cell::UnitCell)
    return unit_cell.symmetry
end

"""
    centering(unit_cell::UnitCell) -> Centering

Return the centering of `unit_cell`.

Return values
=============
centering
"""
@inline function centering(unit_cell::UnitCell)
    return centering(unit_cell.symmetry)
end

"""
    symmetry_elements(unit_cell::UnitCell) -> Vector{SymmetryElement}

Return the symmetry elements of `unit_cell`.

Return values
=============
symmetry elements
"""
@inline function symmetry_elements(unit_cell::UnitCell)
    return symmetry_elements(unit_cell.symmetry)
end

"""
    is_bravais_lattice(unit_cell::UnitCell) -> Bool

Determine if the unit cell defined by `unit_cell` is a valid Bravais lattice type.

Return values
=============
`true` if the lattice system and centering of `unit_cell` define a valid Bravais lattice
type; `false` otherwise

Examples
========
```jldoctest
julia> is_bravais_lattice(TetragonalUnitCell(2, 3; centering=primitive_centering))
true

julia> is_bravais_lattice(TetragonalUnitCell(2, 3; centering=face_centering))
false
```
"""
@inline function is_bravais_lattice(unit_cell::UnitCell)
    return is_bravais_lattice(lattice_system(unit_cell), centering(unit_cell))
end

"""
    basis(unit_cell::UnitCell) -> (Vector{Float64}, Vector{Float64}, Vector{Float64})

Return a set of basis vectors ``\\vec{a}, \\vec{b}, \\vec{c}`` for `unit_ cell`.

Return values
=============
basis vectors ``\\vec{a}``, ``\\vec{b}``, ``\\vec{c}``

Examples
========
```jldoctest
julia> B = basis(UnitCell([1, 0, 0], [1, 1, 0], [1, 0, 2]; centering=P_centering));

julia> B[1] ≈ [1.0, 0.0, 0.0]
true
julia> B[2] ≈ [1.0, 1.0, 0.0]
true
julia> B[3] ≈ [1.0, 0.0, 2.0]
true
```
"""
function basis end

"""
    volume(unit_cell::UnitCell) -> Float64

Compute the volume of the unit cell defined by `unit_cell`.

Return values
=============
volume of the unit cell

Examples
========
```jldoctest
julia> volume(UnitCell([1, 0, 0], [1, 1, 0], [1, 0, 2]; centering=P_centering))
2.0
```
"""
function volume end

"""
    surface_area(unit_cell::UnitCell) -> Float64

Compute the surface area of the unit cell defined by `unit_cell`.

Return values
=============
surface area of the unit cell

Examples
========
```jldoctest; filter = r"(\\d*)\\.(\\d{4})\\d+" => s"\\1.\\2***"
julia> surface_area(UnitCell([1, 0, 0], [1, 1, 0], [1, 0, 1]; centering=P_centering))
7.464101615137754
```
"""
function surface_area end

"""
    standardize(unit_cell::UnitCell) -> UnitCell

Standardize the lattice constants and centering for `unit_cell`.

!!! note

    This function _only_ enforces lattice constant constraints. It _does not modify_ the
    Bravais lattice type. To find an equivalent Bravais lattice with higher symmetry (if
    one exists), use `conventional_cell()`.

!!! note

    Lattice constant standardizations are based on the conventions provided in the Table
    3.1.4.1. of the International Tables for Crystallography (2016).

    * For triclinic lattices, the lattice constants are standardized using the following
      conventions:

      - `a` ≤ `b` ≤ `c`

      - all three angles are acute (Type I cell) or all three angles are non-acute (Type II
        cell)

      - angles sorted in increasing order when edge lengths are equal
          - `α` ≤ `β` when `a` = `b`
          - `β` ≤ `γ` when `b` = `c`
          - `α` ≤ `β` ≤ `γ` when `a` = `b` = `c`

    * For monoclinic lattices, the lattice constants are standardized using the following
      conventions:

      - `a` ≤ `c`

      - π/2 ≤ β ≤ π

    * For orthorhombic lattices, the lattice constants are standardized using the following
      conventions:

      - Primitive, body-centered, and face-centered unit cells: `a` ≤ `b` ≤ `c`

      - Base-C centered unit cells: `a` ≤ `b`, no constraints on `c`

!!! note

    Except for monoclinic lattices, `standardize()` does not modify the unit cell.

    For monoclinic lattices, the unit cell may be potentially modified in two ways:

    - the 2D unit cell in the plane normal to the __b__-axis may be reduced in order to
      satisfy the IUCr conventions for `a`, `c`, and `β`;

    - base-centered unit cells are transformed to equivalent body-centered unit cells.

Return values
=============
unit cell with standardized lattice constants and centering

Examples
========
```jldoctest
julia> unit_cell = OrthorhombicUnitCell(3, 2, 1; centering=P_centering)
OrthorhombicUnitCell((a = 3, b = 2, c = 1), UnitCellSymmetry(PrimitiveCentering(), Set{SymmetryElement}()))
julia> standardize(unit_cell)
OrthorhombicUnitCell((a = 1, b = 2, c = 3), UnitCellSymmetry(PrimitiveCentering(), Set{SymmetryElement}()))

julia> unit_cell = OrthorhombicUnitCell(3, 2, 1; centering=base_centering)
OrthorhombicUnitCell((a = 3, b = 2, c = 1), UnitCellSymmetry(BaseCentering(), Set{SymmetryElement}()))
julia> standardize(unit_cell)
OrthorhombicUnitCell((a = 2, b = 3, c = 1), UnitCellSymmetry(BaseCentering(), Set{SymmetryElement}()))
```
"""
function standardize(unit_cell::UnitCell)
    # --- Check arguments

    standardize_check_args(unit_cell)

    # --- By default, no standardization is performed

    return unit_cell
end

function standardize_check_args(unit_cell::UnitCell)
    # --- Check arguments

    if !is_bravais_lattice(lattice_system(unit_cell), centering(unit_cell))
        throw(
            ArgumentError(
                "Invalid Bravais lattice: " *
                "(lattice_system=$(nameof(typeof(lattice_system(unit_cell)))), " *
                "centering=$(nameof(typeof(centering(unit_cell)))))",
            ),
        )
    end
end

"""
    conventional_cell(unit_cell::UnitCell) -> UnitCell

Return the IUCr conventional cell that is equivalent to `unit_cell`.

Return values
=============
IUCr conventional cell for `unit_cell`

Examples
========
```jldoctest
julia> unit_cell = UnitCell([1, 1, 0], [1, -1, 0], [0, 1, 1]; centering=P_centering);

julia> lattice_system(unit_cell)
Triclinic()

julia> conventional_cell(unit_cell) ≈ CubicUnitCell(2.0; centering=F_centering)
true
```
"""
function conventional_cell(unit_cell::UnitCell)
    # --- Check arguments

    conventional_cell_check_args(unit_cell)

    # --- Compute IUCr conventional cells for lattice systems with limiting cases that
    #     change the Bravais lattice type

    lattice_system_ = lattice_system(unit_cell)
    if lattice_system_ === triclinic
        return conventional_cell(triclinic, unit_cell)
    elseif lattice_system_ === monoclinic
        return conventional_cell(monoclinic, unit_cell)
    elseif lattice_system_ === orthorhombic
        return conventional_cell(orthorhombic, unit_cell)
    elseif lattice_system_ === rhombohedral
        return conventional_cell(rhombohedral, unit_cell)
    elseif lattice_system_ === tetragonal
        return conventional_cell(tetragonal, unit_cell)
    end

    # --- By default, the unit cell is unchanged

    return unit_cell
end

function conventional_cell_check_args(unit_cell::UnitCell)
    # --- Check arguments

    if !is_bravais_lattice(unit_cell)
        throw(
            ArgumentError(
                "Invalid Bravais lattice: " *
                "(lattice_system=" *
                "$(nameof(typeof(lattice_system(unit_cell)))), " *
                "centering=$(nameof(typeof(centering(unit_cell)))))",
            ),
        )
    end
end

"""
    reduced_cell(unit_cell::UnitCell) -> UnitCell

Compute the primitive reduced cell for `unit_cell`. The Selling-Delaunay reduction
algorithm is used to compute the reduced basis.

Return values
=============
primitive reduced cell

Examples
========
```jldoctest
julia> reduced_cell(UnitCell([1, 0, 0], [1, 1, 0], [0, 0, 2]; centering=P_centering))
TetragonalUnitCell((a = 1.0, c = 2.0), UnitCellSymmetry(PrimitiveCentering(), Set{SymmetryElement}()))
```
"""
function reduced_cell(unit_cell::UnitCell)
    # --- Preparations

    # Get basis of unit cell
    basis_a, basis_b, basis_c = basis(unit_cell)

    # Construct primitive cell basis for centering
    if centering(unit_cell) == base_centering
        basis_a = 0.5 * (basis_a + basis_b)

    elseif centering(unit_cell) == body_centering
        basis_c = 0.5 * (basis_a + basis_b + basis_c)

    elseif centering(unit_cell) == face_centering
        basis_a_primitive = 0.5 * (basis_a + basis_b)
        basis_b_primitive = 0.5 * (basis_a - basis_b)
        basis_c_primitive = 0.5 * (basis_a + basis_c)
        basis_a = basis_a_primitive
        basis_b = basis_b_primitive
        basis_c = basis_c_primitive
    end

    # --- Perform Selling-Delaunay reduction

    # Compute Delaunay set
    delaunay_set = compute_delaunay_set(basis_a, basis_b, basis_c)

    # Remove (1) zero vectors and (2) vectors that are scalar multiples of another vector
    # with shorter length.
    delaunay_set = prune_delaunay_set(delaunay_set)

    # Find reduced basis (smallest sum of squared lengths with maximum surface area)
    reduced_basis_a, reduced_basis_b, reduced_basis_c = find_reduced_basis(delaunay_set)

    # --- Return standardized primitive unit cell defined by reduced basis

    return standardize(
        UnitCell(
            reduced_basis_a, reduced_basis_b, reduced_basis_c; centering=primitive_centering
        ),
    )
end

"""
    compute_delaunay_set(
        basis_a::Vector{<:Real}, basis_b::Vector{<:Real}, basis_c::Vector{<:Real}
    ) -> Vector{Vector{Real}}

Compute Delaunay set associated with `basis_a`, `basis_b`, and `basis_c`.

Return values
=============
Delaunay set
"""
function compute_delaunay_set(
    basis_a::Vector{<:Real}, basis_b::Vector{<:Real}, basis_c::Vector{<:Real}
)

    # Check arguments
    if length(basis_a) != 3
        throw(
            ArgumentError("`basis_a` must contain exactly 3 components (basis_a=$basis_a)")
        )
    end

    if length(basis_b) != 3
        throw(
            ArgumentError("`basis_b` must contain exactly 3 components (basis_b=$basis_b)")
        )
    end

    if length(basis_c) != 3
        throw(
            ArgumentError("`basis_c` must contain exactly 3 components (basis_c=$basis_c)")
        )
    end

    # Initialize reduced vector set
    reduced_set = [basis_a, basis_b, basis_c, -(basis_a + basis_b + basis_c)]

    # Reduce the sum of the squares of the lengths of the candidate basis vectxors
    while true
        # Compute scalar products
        scalar_products = [
            dot(reduced_set[1], reduced_set[2]),
            dot(reduced_set[1], reduced_set[3]),
            dot(reduced_set[1], reduced_set[4]),
            dot(reduced_set[2], reduced_set[3]),
            dot(reduced_set[2], reduced_set[4]),
            dot(reduced_set[3], reduced_set[4]),
        ]

        # Check if the reduction is complete (when none of the scalar products is positive)
        if all(scalar_products .<= ALMOST_ZERO)
            break
        end

        # Find the largest positive scalar product
        _, idx = findmax(scalar_products)
        if idx == 1
            b_1 = reduced_set[1]
            b_2 = reduced_set[2]
            b_3 = reduced_set[3]
            b_4 = reduced_set[4]
        elseif idx == 2
            b_1 = reduced_set[1]
            b_2 = reduced_set[3]
            b_3 = reduced_set[2]
            b_4 = reduced_set[4]
        elseif idx == 3
            b_1 = reduced_set[1]
            b_2 = reduced_set[4]
            b_3 = reduced_set[2]
            b_4 = reduced_set[3]
        elseif idx == 4
            b_1 = reduced_set[2]
            b_2 = reduced_set[3]
            b_3 = reduced_set[1]
            b_4 = reduced_set[4]
        elseif idx == 5
            b_1 = reduced_set[2]
            b_2 = reduced_set[4]
            b_3 = reduced_set[1]
            b_4 = reduced_set[3]
        else
            b_1 = reduced_set[3]
            b_2 = reduced_set[4]
            b_3 = reduced_set[1]
            b_4 = reduced_set[2]
        end

        # Update reduced basis set to reduce the sum of the squares of the lengths
        reduced_set = [-b_1, b_2, b_1 + b_3, b_1 + b_4]
    end

    # Compute Delaunay set
    delaunay_set = [
        reduced_set[1],
        reduced_set[2],
        reduced_set[3],
        reduced_set[4],
        reduced_set[1] + reduced_set[2],
        reduced_set[1] + reduced_set[3],
        reduced_set[2] + reduced_set[3],
    ]

    return delaunay_set
end

"""
    prune_delaunay_set(delaunay_set::Vector{<:Vector{<:Real}})

Remove the following vectors from `delaunay_set`:

* zero vectors

and

* vectors that are scalar multiples of another vector with shorter length.

Return values
=============
pruned Delaunay set with each vector `v` augmented by the squared length of `v`. The
returned list contains NamedTuple objects with the following keys: `vector` and
`length_sq`.

!!! note

    This function converts all vectors in `delaunay_set` to `Vector{Float64}`.
"""
function prune_delaunay_set(delaunay_set::Vector{<:Vector{<:Real}})

    # Compute squared length of each vector
    delaunay_set = [
        (vector=convert(Vector{Float64}, v), length_sq=Float64(dot(v, v))) for
        v in delaunay_set
    ]

    # Sort vectors by squared length
    sort!(delaunay_set; by=item -> item.length_sq)

    # Remove zero vectors
    deleteat!(delaunay_set, findall(x->abs(x.length_sq) < ALMOST_ZERO, delaunay_set))

    # Remove vectors that are scalar multiples of shorter vectors
    to_remove = []
    for i in 1:length(delaunay_set)
        basis_1 = delaunay_set[i].vector
        length_basis_1 = sqrt(delaunay_set[i].length_sq)

        for j in (i + 1):length(delaunay_set)
            basis_2 = delaunay_set[j].vector
            length_basis_2 = sqrt(delaunay_set[j].length_sq)

            if abs(dot(basis_1, basis_2)) / length_basis_1 / length_basis_2 ≈ 1.0
                push!(to_remove, j)
            end
        end
    end

    delaunay_set = [delaunay_set[i] for i in 1:length(delaunay_set) if !(i in to_remove)]
end

"""
    find_reduced_basis(
        delaunay_set::Vector{@NamedTuple{vector::Vector{Float64}, length_sq::Float64}}
    ) -> Vector{Vector{Float64}}

Identify reduced basis for `delaunay_set`.

Return values
=============
reduced basis vectors sorted by length (in ascending order)
"""
function find_reduced_basis(
    delaunay_set::Vector{@NamedTuple{vector::Vector{Float64},length_sq::Float64}}
)
    # Check arguments
    if length(delaunay_set) < 3
        throw(
            ArgumentError(
                "`delaunay_set` must contain at least 3 elements" *
                "(delaunay_set=$delaunay_set)",
            ),
        )
    end

    # Generate all possible bases that can be formed from the candidate basis vectors
    # and compute the (1) sum of squared lengths of the basis vectors and (2) the
    # surface area of the unit cell
    basis_choices = [
        (
            vectors=[
                candidate_basis[1].vector,
                candidate_basis[2].vector,
                candidate_basis[3].vector,
            ],
            sum_length_sq=(
                candidate_basis[1].length_sq +
                candidate_basis[2].length_sq +
                candidate_basis[3].length_sq
            ),
            surface_area=surface_area(
                candidate_basis[1].vector,
                candidate_basis[2].vector,
                candidate_basis[3].vector,
            ),
        ) for candidate_basis in combinations(delaunay_set, 3)
    ]

    # Eliminate basis choices that are degenerate (i.e., not linearly independent)
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

    return reduced_basis_a, reduced_basis_b, reduced_basis_c
end

"""
    is_equivalent(
        unit_cell_test::UnitCell,
        unit_cell_ref::UnitCell;
        atol::Real=1e-3,
        rtol::Real=atol > 0 ? 0 : 1e-3,
        p::Real=2,
    ) -> Bool

Check if the unit cell defined by `unit_cell_test` is equivalent to the unit cell defined
by `unit_cell_ref`.

Keyword Arguments
=================
- `atol`: tolerance of the absolute difference between the reduced unit cells defined by
  `unit_cell_test` and `unit_cell_ref`.

- `rtol`: tolerance of the relative difference between the reduced unit cells defined by
  `unit_cell_test` and `unit_cell_ref`.

Return values
=============
`true` if the test unit cell is equivalent to the reference unit cell; `false` otherwise

Examples
========
```jldoctest
julia> unit_cell_ref = UnitCell([1, 0, 0], [1, 1, 0], [0, 0, 2]);

julia> unit_cell_test = TetragonalUnitCell(1.0, 2.0);

julia> is_equivalent(unit_cell_test, unit_cell_ref)
true
```
"""
function is_equivalent(
    unit_cell_test::UnitCell,
    unit_cell_ref::UnitCell;
    atol::Real=1e-3,
    rtol::Real=atol > 0 ? 0 : 1e-3,
)
    # --- Check arguments

    if atol < 0
        throw(DomainError(atol, "`atol` must be nonnegative"))
    end

    if rtol < 0
        throw(DomainError(rtol, "`rtol` must be nonnegative"))
    end

    # --- Compare reduced cells

    return isapprox(
        reduced_cell(unit_cell_test), reduced_cell(unit_cell_ref); atol=atol, rtol=rtol
    )
end

"""
    is_supercell(
        unit_cell_test::UnitCell,
        unit_cell_ref::UnitCell;
        tol::Real=1e-3,
        max_index::Integer=3
    ) -> Bool

Check if the unit cell defined by `unit_cell_test` is a supercell of the unit cell
defined by `unit_cell_ref`.

Keyword Arguments
=================
- `tol`: absolute tolerance of the deviation between test unit cell and supercells of the
  reference unit cell

- `max_index`: maximum multiple of the basis vectors of the reference unit cell to include
  when checking whether the test unit cell is a supercell of the reference unit cell

  !!! note

      `max_index` is ignored for Cubic unit cells.

Return values
=============
`true` if the test unit cell is a supercell of the reference unit cell; `false` otherwise

Examples
========
```jldoctest
julia> unit_cell_ref = UnitCell([1, 0, 0], [0, 1, 0], [0, 0, 1]);

julia> is_supercell(CubicUnitCell(2), unit_cell_ref)
true

julia> is_supercell(CubicUnitCell(2.5), unit_cell_ref)
false

julia> is_supercell(UnitCell([1, 0, 0], [0, 2, 0], [0, 0, 3]), unit_cell_ref)
false
```
"""
function is_supercell(::UnitCell, ::UnitCell)
    # Default is_supercell() implementation to allow comparison between unit cells for
    # different lattice systems.
    return false
end

function Base.:(==)(x::UnitCell, y::UnitCell)
    # Default :(==) implementation to allow comparison between unit cells for different
    # lattice systems.
    return false
end

function Base.:(==)(x::UnitCell{T}, y::UnitCell{T}) where {T<:LatticeSystem}
    return (
        symmetry(x) == symmetry(y) && all(
            getfield(x.lattice_constants, name) == getfield(y.lattice_constants, name) for
            name in fieldnames(typeof(x.lattice_constants))
        )
    )
end

"""
    isapprox(x::UnitCell, y::UnitCell;
             atol::Real=0, rtol::Real=atol>0 ? 0 : √eps)

Inexact equality comparison between `UnitCell`. Two unit cells are approximately equal if
(1) their lattice constant values are equal to within the tolerance bounds and (2) they
have the same symmetry. For instance, `isapprox` returns `true` for
`TetragonalUnitCell` if
`isapprox(lattice_constants(x).a - lattice_constants(y).a; atol=atol, rtol=rtol)` and
`isapprox(lattice_constants(x).b - lattice_constants(y).b; atol=atol, rtol=rtol)` and
`symmetry(x) == symmetry(y)`.  Returns `false` if `x` and `y` have different types.
"""
function isapprox(x::UnitCell, y::UnitCell; atol::Real=0, rtol::Real=atol > 0 ? 0 : √eps())
    # Default isapprox() implementation to allow comparison between unit cells for
    # different lattice systems.
    return false
end

function isapprox(
    x::UnitCell{T}, y::UnitCell{T}; atol::Real=0, rtol::Real=atol > 0 ? 0 : √eps()
) where {T<:LatticeSystem}
    return (
        symmetry(x) == symmetry(y) && all(
            isapprox(
                getfield(x.lattice_constants, name),
                getfield(y.lattice_constants, name);
                atol=atol,
                rtol=rtol,
            ) for name in fieldnames(typeof(x.lattice_constants))
        )
    )
end

using DataStructures: OrderedDict
function Base.:(-)(x::UnitCell{T}, y::UnitCell{T}) where {T<:LatticeSystem}
    Δlattice_constants = OrderedDict([
        (
            Symbol("Δ$name"),
            getfield(x.lattice_constants, name) - getfield(y.lattice_constants, name),
        ) for name in fieldnames(typeof(x.lattice_constants))
    ])
    return UnitCellDelta{T}(
        NamedTuple{Tuple(keys(Δlattice_constants))}(values(Δlattice_constants))
    )
end
