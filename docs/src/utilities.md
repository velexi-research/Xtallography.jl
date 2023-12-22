```@meta
CurrentModule = XtallographyUtils
```

# Utility Functions

-------------------------------------------------------------------------------------------
## Math Functions

```@docs
asin_
acos_
is_basis
volume(::Vector{<:Real}, ::Vector{<:Real}, ::Vector{<:Real})
surface_area(::Vector{<:Real}, ::Vector{<:Real}, ::Vector{<:Real})
```
-------------------------------------------------------------------------------------------
## LatticeConstant Functions

```@docs
lattice_system
standardize(::LatticeConstants, ::Centering)
standardize(::LatticeConstants)
```
-------------------------------------------------------------------------------------------
## Unit Cell Function

```@docs
standardize(::UnitCell)
basis
volume(::UnitCell)
surface_area(::UnitCell)
conventional_cell
reduced_cell
is_equivalent_unit_cell
is_supercell
```

Triclinic Unit Cell Functions
```@docs
is_triclinic_type_I_cell
convert_to_body_centering
convert_to_base_centering
convert_to_mP
convert_to_mI
convert_to_mC
satisfies_triclinic_angle_constraints
```
-------------------------------------------------------------------------------------------
