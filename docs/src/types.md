```@meta
CurrentModule = Xtallography
```

# Types

__Lattice Types__
* `LatticeSystem`
* `Centering`

__Unit Cell Types__
* `LatticeConstants`
* `LatticeConstantDeltas`
* `UnitCell`

-------------------------------------------------------------------------------------------
## Lattice Types

```@docs
LatticeSystem
Centering
```

### Concrete Types

#### Lattice Systems

All subtypes of `LatticeSystem` are singleton types. For convenience, a singleton instance
is defined for each lattice system type.

```@docs
Triclinic
Monoclinic
Orthorhombic
Hexagonal
Rhombohedral
Tetragonal
Cubic
```

#### Centerings

All subtypes of `Centering` are singleton types. For convenience, a singleton instance is
defined for each centering type.

```@docs
Primitive
BaseCentered
BodyCentered
FaceCentered
```

### Constants

#### Lattice Systems

```@docs
triclinic
monoclinic
orthorhombic
hexagonal
rhombohedral
tetragonal
cubic
```

#### Centerings

```@docs
primitive
base_centered
body_centered
face_centered
```

#### Other Constants

```@docs
BRAVAIS_LATTICES
```

-------------------------------------------------------------------------------------------
## Unit Cell Types

```@docs
LatticeConstants
LatticeConstantDeltas
UnitCell
```

### Lattice Constants

```@docs
TriclinicLatticeConstants
MonoclinicLatticeConstants
OrthorhombicLatticeConstants
HexagonalLatticeConstants
RhombohedralLatticeConstants
TetragonalLatticeConstants
CubicLatticeConstants
```

### Lattice Constant Deltas

```@docs
TriclinicLatticeConstantDeltas
MonoclinicLatticeConstantDeltas
OrthorhombicLatticeConstantDeltas
HexagonalLatticeConstantDeltas
RhombohedralLatticeConstantDeltas
TetragonalLatticeConstantDeltas
CubicLatticeConstantDeltas
```

-------------------------------------------------------------------------------------------
