```@meta
CurrentModule = Xtallography
```

# Unit Cells

Types and functions related to unit cells.

-------------------------------------------------------------------------------------------
## Unit Cells

```@docs
UnitCell
basis
centering(::UnitCell)
conventional_cell
isapprox(::UnitCell, ::UnitCell)
is_equivalent
is_supercell
lattice_constants
reduced_cell
standardize
surface_area
symmetry
volume
```

### Concrete Types

```@docs
CubicUnitCell
CubicUnitCell(
    a::Real;
    centering::Centering=primitive_centering,
    symmetry_elements::Union{Set,Vector,Nothing}=nothing,
)
HexagonalUnitCell
HexagonalUnitCell(
    a::Real,
    c::Real;
    centering::Centering=primitive_centering,
    symmetry_elements::Union{Set,Vector,Nothing}=nothing,
)
OrthorhombicUnitCell
OrthorhombicUnitCell(
    a::Real,
    b::Real,
    c::Real;
    centering::Centering=primitive_centering,
    symmetry_elements::Union{Set,Vector,Nothing}=nothing,
)
MonoclinicUnitCell
MonoclinicUnitCell(
    a::Real,
    b::Real,
    c::Real,
    β::Real;
    centering::Centering=primitive_centering,
    symmetry_elements::Union{Set,Vector,Nothing}=nothing,
)
RhombohedralUnitCell
RhombohedralUnitCell(
    a::Real,
    α::Real;
    centering::Centering=primitive_centering,
    symmetry_elements::Union{Set,Vector,Nothing}=nothing,
)
TetragonalUnitCell
TetragonalUnitCell(
    a::Real,
    c::Real;
    centering::Centering=primitive_centering,
    symmetry_elements::Union{Set,Vector,Nothing}=nothing,
)
TriclinicUnitCell
TriclinicUnitCell(
    a::Real,
    b::Real,
    c::Real,
    α::Real,
    β::Real,
    γ::Real;
    centering::Centering=primitive_centering,
    symmetry_elements::Union{Set,Vector,Nothing}=nothing,
    check_angle_constraints::Bool=true,
)
```

### Lattice-Specific Functions

#### Triclinic Unit Cell Functions

```@docs
convert_to_mC
convert_to_mI
convert_to_mP
is_triclinic_type_I_cell
satisfies_triclinic_angle_constraints
```

#### Monoclinic Unit Cell Functions

```@docs
convert_to_base_centering
convert_to_body_centering
```

-------------------------------------------------------------------------------------------
## Unit Cell Deltas

```@docs
UnitCellDelta
Δlattice_constants
isapprox(::UnitCellDelta, ::UnitCellDelta)
```

### Concrete Types

```@docs
CubicUnitCellDelta
CubicUnitCellDelta(Δa::Real)
HexagonalUnitCellDelta
HexagonalUnitCellDelta(Δa::Real, Δc::Real)
MonoclinicUnitCellDelta
MonoclinicUnitCellDelta(Δa::Real, Δb::Real, Δc::Real, Δβ::Real)
OrthorhombicUnitCellDelta
OrthorhombicUnitCellDelta(Δa::Real, Δb::Real, Δc::Real)
RhombohedralUnitCellDelta
RhombohedralUnitCellDelta(Δa::Real, Δα::Real)
TetragonalUnitCellDelta
TetragonalUnitCellDelta(Δa::Real, Δc::Real)
TriclinicUnitCellDelta
TriclinicUnitCellDelta(Δa::Real, Δb::Real, Δc::Real, Δα::Real, Δβ::Real, Δγ::Real)
```

-------------------------------------------------------------------------------------------
## Unit Cell Symmetry

```@docs
UnitCellSymmetry
centering(::UnitCellSymmetry)
primitive_unit_cell_symmetry
symmetry_elements
```

-------------------------------------------------------------------------------------------
