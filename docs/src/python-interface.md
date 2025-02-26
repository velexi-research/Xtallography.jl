# Python Interface

[Xtallography](https://github.com/velexi-research/Xtallography.jl) provides a Python
interface to facilitate integration with Python codebases. See the
[xtallography](../python/) documentation for details.

## Installation

* Add the Velexi Julia package registry to the Julia installation.

  ```julia
  julia>  # Press ']' to enter the Pkg REPL mode.
  pkg> registry add https://github.com/velexi-research/JuliaRegistry.git
  ```

  !!! note "`Xtallgraphy` is registered with a local Julia package registry"

      The `xtallography` Python package depends on the `Xtallography` Julia package, which
      is registered with the Velexi Julia package registry (not the General Julia package
      registry). The Julia package manager used by the `xtallography` Python package will
      be able to find `Xtallography` only if the Velexi Julia package registry has been
      added to your Julia installation. For more information about local registries for
      Julia packages,
      [LocalRegistry.jl](https://github.com/GunnarFarneback/LocalRegistry.jl).

  !!! tip "Only needed once"

      This step only needs to be performed once per Julia installation.

* Install the `xtallography` package from the GitHub repository

  !!! tip "No need to manually install the `Xtallography.jl` Julia package"

      The Julia package dependencies for `xtallography` are installed and updated
      automatically whenever the package is used.

  * `poetry`

    ```shell
    $ poetry add git+https://github.com/velexi-research/Xtallography.jl.git
    ```

  * `pip`

    ```shell
    $ pip install git+https://github.com/velexi-research/Xtallography.jl.git
    ```
