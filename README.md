Xtallography Utilities
===============================================================================

[------------------------------------ BADGES: BEGIN ------------------------------------]: #

<table>
  <tr>
    <td>Documentation</td>
    <td>
      <a href="https://velexi-research.github.io/XtallographyUtils.jl/dev/"><img style="vertical-align: bottom;" src="https://img.shields.io/badge/docs-dev-blue.svg"/></a>
      <a href="https://velexi-research.github.io/XtallographyUtils.jl/stable/"><img style="vertical-align: bottom;" src="https://img.shields.io/badge/docs-stable-blue.svg"/></a>
    </td>
  </tr>

  <tr>
    <td>Build Status</td>
    <td>
      <a href="https://github.com/velexi-research/XtallographyUtils.jl/actions/workflows/CI.yml"><img style="vertical-align: bottom;" src="https://github.com/velexi-research/XtallographyUtils.jl/actions/workflows/CI.yml/badge.svg"/></a>
      <a href="https://codecov.io/gh/velexi-research/XtallographyUtils.jl">
        <img style="vertical-align: bottom;" src="https://codecov.io/gh/velexi-research/XtallographyUtils.jl/graph/badge.svg?token=BZKPGII992"/></a>
    </td>
  </tr>

  <!-- Miscellaneous Badges -->
  <tr>
    <td colspan=2 align="center">
      <a href="https://github.com/velexi-research/XtallographyUtils.jl/issues"><img style="vertical-align: bottom;" src="https://img.shields.io/badge/contributions-welcome-brightgreen.svg?style=flat"/></a>
      <a href="https://github.com/invenia/BlueStyle"><img style="vertical-align: bottom;" src="https://img.shields.io/badge/code%20style-blue-4495d1.svg"/></a>
    </td>
  </tr>
</table>

[------------------------------------- BADGES: END -------------------------------------]: #

The Xtallography Utilities package defines a collection of basic types and miscellaneous
utility functions that support crystallography computations. The functionality provided by
this package are _not_ intended to be comprehensive. Rather, functionality is added on an
as-needed basis to support research projects.

Currently, the Xtallography Utilities provides support for:

* types for defining Bravais lattice types,

* basic unit cell computations (e.g., basis, volume, surface area),

* standardization of lattice constants for unit cells, and

* conversions between equivalent unit cells for a lattice.

## Getting Started

* Add the Velexi Julia package registry.

  ```julia
  julia>  # Press ']' to enter the Pkg REPL mode.
  pkg> registry add https://github.com/velexi-research/JuliaRegistry.git
  ```

  ___Notes___

  * `XtallographyUtils` is registered with the Velexi Julia package registry (not the
    General Julia package registry), so the Pkg REPL will be able to find
    `XtallographyUtils` only if the Velexi Julia package registry has been added to your
    Julia installation. For more information about local registries for Julia packages,
    [LocalRegistry.jl](https://github.com/GunnarFarneback/LocalRegistry.jl)

  * This step only needs to be performed once per Julia installation.

* Install the `XtallographyUtils` package via the Pkg REPL. That's it!

  ```julia
  julia>  # Press ']' to enter the Pkg REPL mode.
  pkg> add XtallographyUtils
  ```

## Examples

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

## Related Packages

There are a couple of crystallography packages in the Julia ecosystem that provide support
for various crystall and lattice computations.

* [Crystallography.jl](https://github.com/MineralsCloud/Crystallography.jl) and related
  packages

* [Spglib](https://github.com/singularitti/Spglib.jl)
