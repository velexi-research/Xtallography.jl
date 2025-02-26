# Xtallography

The [Xtallography](https://github.com/velexi-research/Xtallography.jl) package defines
a collection of basic types and functions that support crystallography computations. The
functionality provided by this package are _not_ intended to be comprehensive. Rather,
functionality is added on an as-needed basis to support research projects.

Currently, Xtallography provides support for:

* types for defining Bravais lattice types,

* basic unit cell computations (e.g., basis, volume, surface area),

* standardization of lattice constants for unit cells,

* conversions between equivalent unit cells for a lattice.

In addition, Xtallography provides a basic [Python interface](python-interface/) to support
integration with Python codebases.

## Getting Started

* Add the Velexi Julia package registry.

  ```julia
  julia>  # Press ']' to enter the Pkg REPL mode.
  pkg> registry add https://github.com/velexi-research/JuliaRegistry.git
  ```

  !!! note "Xtallography is registered with a local Julia package registry"

      `Xtallography` is registered with the Velexi Julia package registry (not the General
      Julia package registry), so the Pkg REPL will be able to find `Xtallography` only if
      the Velexi Julia package registry has been added to your Julia installation. For more
      information about local registries for Julia packages,
      [LocalRegistry.jl](https://github.com/GunnarFarneback/LocalRegistry.jl).

  !!! tip "Only needed once"

      This step only needs to be performed once per Julia installation.

* Install the `Xtallography` package via the Pkg REPL. That's it!

  ```julia
  julia>  # Press ']' to enter the Pkg REPL mode.
  pkg> add Xtallography
  ```

## Related Packages

There are a couple of crystallography packages in the Julia ecosystem that provide support
for various crystallography and lattice computations.

* [Crystallography.jl](https://github.com/MineralsCloud/Crystallography.jl) and related
  packages

* [Spglib](https://github.com/singularitti/Spglib.jl)
