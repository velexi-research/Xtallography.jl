Xtallography Release Notes
============================================================================================

--------------------------------------------------------------------------------------------
0.5.3 (2025-04-27)
==================
**Bug Fixes**
- Add Python version constraint to fix `juliacall` incompatibility.

--------------------------------------------------------------------------------------------
0.5.2 (2025-04-27)
==================
**Bug Fixes**
- Fixed version bug in `juliapkg.json`.

--------------------------------------------------------------------------------------------
0.5.1 (2025-03-21)
==================
**Enhancements**
- Add general Julia `UnitCell.from_julia()` method to convert Julia `UnitCell` objects
  of unknown lattice type to the appropriate Python objects.

--------------------------------------------------------------------------------------------
0.5.0 (2025-03-21)
==================
**Enhancements**
- Upgrade Python interface.
  - Add methods for converting Julia `Centering` and `UnitCell` objects to equivalent
    Python objects.
  - Add comparison methods for `UnitCell` objects.
- Improved error messages.

**Developer Updates**
- Improve robustness of unit tests.
- Update developer reference documents.

--------------------------------------------------------------------------------------------
0.4.5 (2025-02-26)
==================
**Developer Updates**
- Add gh-pages workflow to enable manual deployment of GitHub Pages when automatic
  deployment fails.

--------------------------------------------------------------------------------------------
0.4.4 (2025-02-25)
==================
**Developer Updates**
- No change release to fix Github Pages deployment.

--------------------------------------------------------------------------------------------
0.4.3 (2025-02-25)
==================
**Enhancements**
- Added installation instructions for Python interface.
- Simplified installation of Julia dependencies for Python package.

**Bug Fixes**
- Fixed license incantations in Python source code.

**Developer Updates**
- Added juliapkg.json to configure Julia dependencies of Python package.
- Updated Python packaging configuration.
- Updated CI workflow to use the development version of juliapkg.json.
- Added pre-commit hook to check that juliapkg.json is not a development version.
- Moved developer notes to README.

--------------------------------------------------------------------------------------------
0.4.2 (2025-02-23)
==================
**Bug Fixes**
- Fixed links to Python API in documentation.

**Developer Updates**
- Fixed documentation deployment bugs in GitHub Actions workflows.

--------------------------------------------------------------------------------------------
0.4.1 (2025-02-23)
==================
**Enhancements**
- Improved robustness of Makefile targets.
- Bumped Julia and Python version requirements to improve package reliability.

**Bug Fixes**
- Removed function from Python interface that was accidentally included in v0.4.0.
- Fixed bugs in documentation.

**Developer Updates**
- Added developer notes.
- Updated unit tests.
- Improved robustness of doc generation.
- Updated test matrix in CI workflow.
- Fixed bugs in CI workflow.
- Updated pre-commit configuration.

--------------------------------------------------------------------------------------------
0.4.0 (2025-02-22)
==================
**Enhancements**
- Renamed package to `Xtallography`.
- Added Python interface.
- Refactored lattice and unit cell code.
- Updated package documentation.

**Developer Updates**
- Updated `Project.toml` files.
- Improved CI workflows.
- Added `Aqua` checks to unit tests.
- Updated GitHub Actions workflows.
- Updated pre-commit version and configuration.

--------------------------------------------------------------------------------------------
0.3.2 (2024-05-30)
==================
* Fix bugs in angle constraints for `TriclinicLatticeConstants`,
  `MonoclinicLatticeConstants`, `RhombohedralLatticeConstants` types.
* Update package infrastructure files (`Project.toml`, `CI.yaml`, `.ignore`).

--------------------------------------------------------------------------------------------
0.3.1 (2024-03-01)
==================
* Fix bugs in angle constraints for `TriclinicLatticeConstants`,
  `MonoclinicLatticeConstants`, `RhombohedralLatticeConstants` types.

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
