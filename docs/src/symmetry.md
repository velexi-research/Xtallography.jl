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
LATTICE_SYSTEMS
```

### Concrete Types

```@docs
Cubic
cubic
Hexagonal
hexagonal
Monoclinic
monoclinic
Orthorhombic
orthorhombic
Rhombohedral
rhombohedral
Tetragonal
tetragonal
Triclinic
triclinic
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
```

### Concrete Types

```@docs
BaseCentering
base_centering
BodyCentering
body_centering
I_centering
FaceCentering
face_centering
F_centering
PrimitiveCentering
primitive_centering
P_centering
```

-------------------------------------------------------------------------------------------
## Symmetry Elements

!!! note

    All concrete symmetry element types are subtypes of the `SymmetryElement` abstract
    type.

```@docs
SymmetryElement
```

### Concrete Types

```@docs
GlidePlane
GlidePlane(glide::Tuple{<:Real,<:Real,<:Real},
    normal::Tuple{<:Real,<:Real,<:Real};
    location::Tuple{<:Real,<:Real,<:Real}=(0, 0, 0),
)
InversionCenter
MirrorPlane
MirrorPlane(
    normal::Tuple{<:Real,<:Real,<:Real}; location::Tuple{<:Real,<:Real,<:Real}=(0, 0, 0)
)
RotationAxis
RotationAxis(
    n::Int,
    direction::Tuple{<:Real,<:Real,<:Real};
    location::Tuple{<:Real,<:Real,<:Real}=(0, 0, 0),
)
RotoinversionAxis
RotoinversionAxis(
    n::Int,
    direction::Tuple{<:Real,<:Real,<:Real};
    center::Tuple{<:Real,<:Real,<:Real}=(0, 0, 0),
)
ScrewAxis
ScrewAxis(
    n::Int,
    m::Int,
    direction::Tuple{<:Real,<:Real,<:Real};
    location::Tuple{<:Real,<:Real,<:Real}=(0, 0, 0),
)
```

-------------------------------------------------------------------------------------------
## Bravais Lattices

```@docs
BRAVAIS_LATTICES
is_bravais_lattice
```

-------------------------------------------------------------------------------------------
