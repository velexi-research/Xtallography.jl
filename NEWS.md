XtallographyUtils Release Notes
============================================================================================

--------------------------------------------------------------------------------------------
0.3.0 (2024-02-27)
==================
* Add `LatticeSystem` parameter to `LatticeConstants` and `LatticeConstantDeltas` types.
* Simplify implementation of `lattice_system()` methods.
* Fix minimum angle for rhombohedral unit cells.
* Update package documentation.
* Remove unnecessary imports.

--------------------------------------------------------------------------------------------
0.2.0 (2024-01-19)
==================
* Add `convert(::LatticeConstants)` and `convert(::LatticeConstantDeltas)` methods.
* Change `tol` argument to `is_equivalent_unit_cell()` to `atol` and `rtol` to be more
  consistent with other approximate comparison methods.
* Remove `norm()` methods.
* Update and improve unit tests.
* Update package documentation.
* Add cache for `build-docs` job in CI workflow.

--------------------------------------------------------------------------------------------
0.1.2 (2024-01-18)
==================
* Add `norm(::LatticeConstant)` method.

--------------------------------------------------------------------------------------------
0.1.1 (2024-01-07)
==================
* Add `LatticeConstantDeltas` types and functions.
* Add difference operator for `LatticeConstants` types.

--------------------------------------------------------------------------------------------
0.1.0 (2023-12-26)
==================
* Initial version of package.

--------------------------------------------------------------------------------------------
