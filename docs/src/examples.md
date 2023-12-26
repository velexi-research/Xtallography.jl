# Examples

* Create a unit cell object.

  ```julia
  julia> unit_cell = UnitCell(OrthorhombicLatticeConstants(2, 3, 4), body_centered)
  UnitCell(OrthorhombicLatticeConstants(2.0, 3.0, 4.0), BodyCentered())
  ```

* Compute basic unit cell attributes.

  ```julia
  julia> volume(unit_cell)
  24.0

  julia> surface_area(unit_cell)
  52.0

  julia> basis(unit_cell)
  ([2.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 4.0])
  ```

* Standardize the lattice constants for a unit cell to be consistent with IUCr conventions.

  ```julia
  julia> unit_cell = UnitCell(OrthorhombicLatticeConstants(4, 2, 3), primitive)
  UnitCell(OrthorhombicLatticeConstants(4.0, 2.0, 3.0), Primitive())

  julia> standardize(unit_cell)
  UnitCell(OrthorhombicLatticeConstants(2.0, 3.0, 4.0), Primitive())
  ```

* Compute the Delaunay reduced cell for a unit cell.

  ```julia
  julia> unit_cell = UnitCell(CubicLatticeConstants(4), face_centered)
  UnitCell(CubicLatticeConstants(4.0), FaceCentered())

  julia> r_cell = reduced_cell(unit_cell)
  UnitCell(TriclinicLatticeConstants(2.8284271247461903, 2.8284271247461903, 2.8284271247461903, 1.0471975511965974, 1.0471975511965974, 1.5707963267948966), Primitive())
  ```

* Compute the conventional cell for a unit cell.

  ```julia
  julia> conventional_cell(r_cell)
  UnitCell(CubicLatticeConstants(3.9999999999999982), FaceCentered())
  ```
