Contributing to `Xtallography.jl`
============================================================================================

## Reporting Issues

When reporting issues please include please provide a detailed description of the issue,
incluing a self-contained code snippet that demonstrates the problem.

--------------------------------------------------------------------------------------------
## Contributing Code

[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor%27s%20Guide-blueviolet)](https://github.com/SciML/ColPrac)

Thanks for your interest in contributing code to `Xtallography.jl`! We roughly follow the
[ColPrac guide for collaborative practices](https://docs.sciml.ai/ColPrac/stable/). We
recommend that new contributors read that guide. Below is some additional information for
contributors.

### Package Contents

```
├── README.md          <- README file for package
├── NEWS.md            <- package release notes
├── LICENSE            <- package license
├── NOTICE             <- package copyright notice
├── Makefile           <- Makefile containing useful shortcuts (`make` rules).
│                         Use `make help` to show the list of available rules.
├── Project.toml       <- Julia package metadata file
├── Manifest.toml      <- Julia environment manifest file
├── pyproject.toml     <- Python dependency and configuration file
├── poetry.lock        <- Poetry lockfile
├── docs/              <- package documentation
├── extras/            <- additional files and references that may be useful
│                         for package development
├── spikes/            <- experimental code snippets, etc.
├── src/               <- package source code
└── tests/             <- package test code
```

### Setting Up a Development Environment

___Note___: this project uses `poetry` to manage Python package dependencies.

1. Prerequisites

   * Install [Git][git].

   * Install [Python][python] 3.8 (or greater).

   * Install [Poetry][poetry] 1.2 (or greater).

   * _Optional_. Install [direnv][direnv].

2. ___Recommended___. Set up a dedicated virtual environment for the package
   development using `direnv` (because manages the environment for both the
   shell and Python).

   * Prerequisite. Install `direnv`.

   * Copy `extras/dot-envrc` to the package root directory, and rename it to
     `.envrc`.

     ```shell
     $ cd $PROJECT_ROOT_DIR
     $ cp extras/dot-envrc .envrc
     ```

   * Grant permission to direnv to execute the .envrc file.

     ```shell
     $ direnv allow
     ```

3. Install the Python package dependencies required for development.

   ```shell
   $ poetry install --no-root
   ```

4. Install the Git pre-commit hooks.

   ```shell
   $ pre-commit install
   ```

### Running Automated Tests

This package is configured to support (1) automated testing of code located in
the `src` directory and (2) analysis of how well the tests cover of the source
code (i.e., coverage analysis).

* Run all of the tests.

  ```shell
  $ make test
  ```

* Run all of the tests in fail-fast mode (i.e., stop after the first failing
  test).

  ```shell
  $ make fast-test
  ```

### Cleaning the Development Directory

* Use `make clean` to automatically remove temporary files and directories
  generated during testing (e.g., temporary directories, coverage files).

  ```shell
  $ make clean
  ```

--------------------------------------------------------------------------------------------

[------------------------------------ EXTERNAL LINKS -----------------------------------]: #

[direnv]: https://direnv.net/

[git]: https://git-scm.com/

[python]: https://www.python.org/

[poetry]: https://python-poetry.org/
