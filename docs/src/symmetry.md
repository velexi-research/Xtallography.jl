```@meta
CurrentModule = Xtallography
```

# Symmetry

Types and functions related to defining crystal symmetry.

-------------------------------------------------------------------------------------------
## Lattice Systems

!!! note

    All concrete lattice system types are (1) subtypes of the `LatticeSystem` abstract
    type and (2) singleton types. For convenience, a singleton instance is defined for
    each lattice system type.

```@docs
LatticeSystem
lattice_system
Triclinic
triclinic
Monoclinic
monoclinic
Orthorhombic
orthorhombic
Hexagonal
hexagonal
Rhombohedral
rhombohedral
Tetragonal
tetragonal
Cubic
cubic
```

-------------------------------------------------------------------------------------------
## Centerings

!!! note

    All concrete centering types are (1) subtypes of the `Centering` abstract type and
    (2) singleton types. For convenience, a singleton instance is defined for each
    centering type.

```@docs
Centering
CENTERINGS
PrimitiveCentering
primitive_centering
P_centering
BaseCentering
base_centering
BodyCentering
body_centering
I_centering
FaceCentering
face_centering
F_centering
```

-------------------------------------------------------------------------------------------
## Symmetry Elements

!!! note

    All concrete symmetry element types are subtypes of the `SymmetryElement` abstract
    type.

```@docs
SymmetryElement
GlidePlane
ScrewAxis
```

-------------------------------------------------------------------------------------------
## Bravais Lattices

```@docs
BRAVAIS_LATTICES
is_bravais_lattice
```

-------------------------------------------------------------------------------------------
