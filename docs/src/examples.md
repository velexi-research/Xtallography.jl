# Examples

* Create a unit cell object.

  ```julia
  julia> unit_cell = OrthorhombicUnitCell(2, 3, 4; centering=body_centering)
  OrthorhombicUnitCell((a = 2, b = 3, c = 4), UnitCellSymmetry(BodyCentering(), Set{SymmetryElement}()))
  ```

* Compute basic unit cell attributes.

  ```julia
  julia> volume(unit_cell)
  24

  julia> surface_area(unit_cell)
  52

  julia> basis(unit_cell)
  ([2.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 4.0])
  ```

* Standardize the lattice constants for a unit cell to be consistent with IUCr conventions.

  ```julia
  julia> unit_cell = OrthorhombicUnitCell(4, 2, 3; centering=primitive_centering)
  OrthorhombicUnitCell((a = 4, b = 2, c = 3), UnitCellSymmetry(PrimitiveCentering(), Set{SymmetryElement}()))

  julia> standardize(unit_cell)
  OrthorhombicUnitCell((a = 2, b = 3, c = 4), UnitCellSymmetry(PrimitiveCentering(), Set{SymmetryElement}()))
  ```

* Compute the Delaunay reduced cell for a unit cell.

  ```julia
  julia> unit_cell = CubicUnitCell(4; centering=face_centering)
  CubicUnitCell((a = 4,), UnitCellSymmetry(FaceCentering(), Set{SymmetryElement}()))

  julia> r_cell = reduced_cell(unit_cell)
  TriclinicUnitCell((a = 2.8284271247461903, b = 2.8284271247461903, c = 2.8284271247461903, α = 1.0471975511965974, β = 1.0471975511965974, γ = 1.5707963267948966), UnitCellSymmetry(PrimitiveCentering(), Set{SymmetryElement}()))
  ```

* Compute the conventional cell for a unit cell.

  ```julia
  julia> conventional_cell(r_cell)
  CubicUnitCell((a = 3.9999999999999982,), UnitCellSymmetry(FaceCentering(), Set{SymmetryElement}()))
  ```
