```@meta
CurrentModule = Xtallography
```

# Unit Cells

Types and functions related to unit cells.

-------------------------------------------------------------------------------------------
## Unit Cell

```@docs
UnitCell
UnitCellDelta
UnitCellSymmetry
```

### Constants

```@docs
primitive_unit_cell_symmetry
```

### Properties

#### UnitCell

```@docs
centering
lattice_constants
symmetry
```

#### UnitCellDelta

```@docs
Δlattice_constants
```

#### UnitCellSymmetry

```@docs
symmetry_elements
```

### Unit Cell Geometry Functions

```@docs
basis
surface_area
volume
```

### Unit Cell Standardization and Comparison Functions

```@docs
isapprox(::UnitCell, ::UnitCell)
isapprox(::UnitCellDelta, ::UnitCellDelta)
conventional_cell
is_equivalent_unit_cell
is_supercell
reduced_cell
standardize(::UnitCell)
```

-------------------------------------------------------------------------------------------
## Lattice-Specific Unit Cell Types

```@docs
TriclinicUnitCell
TriclinicUnitCellDelta
MonoclinicUnitCell
MonoclinicUnitCellDelta
OrthorhombicUnitCell
OrthorhombicUnitCellDelta
HexagonalUnitCell
HexagonalUnitCellDelta
RhombohedralUnitCell
RhombohedralUnitCellDelta
TetragonalUnitCell
TetragonalUnitCellDelta
CubicUnitCell
CubicUnitCellDelta
```

-------------------------------------------------------------------------------------------
### Lattice-Specific Unit Cell Functions

#### Triclinic Unit Cell Functions

```@docs
convert_to_mP
convert_to_mI
convert_to_mC
is_triclinic_type_I_cell
satisfies_triclinic_angle_constraints
```

#### Monoclinic Unit Cell Functions

```@docs
convert_to_base_centering
convert_to_body_centering
```

-------------------------------------------------------------------------------------------
