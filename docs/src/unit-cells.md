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
TriclinicUnitCell
MonoclinicUnitCell
OrthorhombicUnitCell
HexagonalUnitCell
RhombohedralUnitCell
TetragonalUnitCell
CubicUnitCell
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
HexagonalUnitCellDelta
MonoclinicUnitCellDelta
OrthorhombicUnitCellDelta
RhombohedralUnitCellDelta
TetragonalUnitCellDelta
TriclinicUnitCellDelta
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
