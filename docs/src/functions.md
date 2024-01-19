```@meta
CurrentModule = XtallographyUtils
```

# Functions

-------------------------------------------------------------------------------------------
## Lattice System Functions

```@docs
is_bravais_lattice
```

-------------------------------------------------------------------------------------------
## Unit Cell Functions

### Unit Cell Properties

```@docs
basis
surface_area(::UnitCell)
volume(::UnitCell)
```

### Unit Cell Standardization and Comparison Functions

```@docs
isapprox(::LatticeConstants, ::LatticeConstants)
conventional_cell
is_equivalent_unit_cell
is_supercell
reduced_cell
standardize(::UnitCell)
standardize(::LatticeConstants, ::Centering)
standardize(::LatticeConstants)
```

### Other Functions

```@docs
lattice_system
isapprox(::LatticeConstantDeltas, ::LatticeConstantDeltas)
convert
```

### Lattice-Specific Functions

#### Triclinic Systems

```@docs
convert_to_mP
convert_to_mI
convert_to_mC
is_triclinic_type_I_cell
satisfies_triclinic_angle_constraints
```

#### Monoclinic Systems

```@docs
convert_to_base_centering
convert_to_body_centering
```

-------------------------------------------------------------------------------------------
## Math Functions

```@docs
asin_
acos_
is_basis
surface_area(::Vector{<:Real}, ::Vector{<:Real}, ::Vector{<:Real})
volume(::Vector{<:Real}, ::Vector{<:Real}, ::Vector{<:Real})
```

-------------------------------------------------------------------------------------------
