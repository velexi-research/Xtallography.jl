Julia Local Registry Guide
==========================

__Authors__  
Kevin T. Chu `<kevin@velexi.com>`

--------------------------------------------------------------------------------------------

Table of Contents
-----------------

1. [Registry Maintainers][#1]

2. [Package Maintainers][#2]

   2.1. [Setting up][#2.1]

   2.2. [Registering a New Package][#2.2]

   2.3. [Releasing a New Version of a Package][#2.3]

3. [Package Users][#3]

4. [References][#4]

--------------------------------------------------------------------------------------------

## 1. Registry Maintainers

* Set up a Git repository for the local registry.

  * __Recommendation__. Use a remote Git repository (e.g., a GitHub repository).

  * Alternative. It is possible to set up a local registry using a local Git repository.

* Create the local registry.

  ```julia
  julia> using LocalRegistry
  julia> create_registry(REGISTRY_NAME,
                         REGISTRY_URL;
                         description="My local Julia registry")
  ```

  The `create_registry()` function creates a registry directory
  `~/.julia/registries/REGISTRY_NAME`, sets it up as a Git repository, and sets the
  remote for the Git repository to `REGISTRY_URL`. Both `REGISTRY_NAME` and `REGISTRY_URL`
  should be specified as strings.

  __Note__. To use SSH to access the remote Git repository for registry, set
  `REGISTRY_URL` to the SSH URL for the remote Git repository (e.g.,
  `git@github.com:your-name/registry-name.git`). To use HTTPS to access the remote Git
  repository for registry, set `REGISTRY_URL` to the HTTPS URL to the remote Git repository
  (e.g., `https://github.com/your-name/registry-name`).

* Push the local registry to remote Git repository.

  ```shell
  $ cd ~/.julia/registries/REGISTRY_NAME
  $ git push -u origin HEAD
  ```

--------------------------------------------------------------------------------------------

## 2. Package Maintainers

___Note___. Package maintainers must have write access to the GitHub repository for the
registry.

### 2.1. Setting Up

* Add the local registry to the Julia depot.

  ```julia
  julia> ]
  pkg> registry add REGISTRY_URL
  ```

  __Recommendation__. Use Julia >=1.7 and set the environment variable
  `JULIA_PKG_USE_CLI_GIT` to true avoid SSH key compatibility issues that are inherent in
  the `libgit2` and `libssh2` packages that the Julia Pkg package uses by default for git
  operations.

### 2.2. Registering a New Package

* Prepare the package for release.

  * Set the initial package version number in the `Project.toml` file.

  * Add release notes for the initial release in the `NEWS.md` file.

  * Merge all repository updates into the `main` branch.

* Use the `LocalRegistry` package to register the package with the local registry.

  ```julia
  julia> using LocalRegistry
  julia> register(package; registry=registry)
  ```

  Both `package` and `registry` may be specified by name as strings. Alternatively,

  * `package` may be specified as (1) a module that has been imported into the REPL or
    (2) a path (which is distinguished from a package name by having more than one path
    component). A path relative to the current working directory can be specified by
    starting the path with `./`;

  * `registry` should be specfied as an HTTPS URL to the registry's remote Git repository
    so that it is readily accessible by Julia (via `add PKG_NAME`).

* Complete the package release process by following the instructions in the last step of
  the ["Releasing a New Version of a Package"][#2.3] section.

### 2.3. Releasing a New Version of a Package

* Prepare the package for release.

  * Increment the package version number in the `Project.toml` file.

  * Update the release notes in the `NEWS.md` file.

  * Merge all repository updates into the `main` branch.

* From the package development directory, use `LocalRegistry` to register a new release of
  the package.

  ```julia
  julia> using LocalRegistry
  julia> register()
  ```

  __Note__. When `register()` is called without a package, the package in the currently
  active project is registered with the registry. See the documentation for `register()`
  for alternative methods of specifying the package.

* Complete the package release process.

  * Tag the registered version using one of the following two methods.

    * Create a tag and release in GitHub.

    * Create a tag manually from the package Git repository.

      ```shell
      $ git tag -a vX.Y.Z -m "<description of version>" <commit hash>
      $ git push origin vX.Y.Z
      ```

  * (Optional) Create a release in GitHub.

--------------------------------------------------------------------------------------------

## 3. Package Users

___Note___. Package users only require read access to the GitHub repository for the
registry.

* Add the local registry to the Julia depot.

  ```julia
  julia> ]
  pkg> registry add REGISTRY_URL
  ```

* Add packages from the local registry.

  ```julia
  julia> ]
  pkg> add PACKAGE_NAME
  ```

--------------------------------------------------------------------------------------------

## 4. References

* [LocalRegistry.jl: Create and maintain local registries for Julia packages][local-registry]

* [LocalRegistry.jl: Migrating Packages from the General Registry][local-registry-migration]

--------------------------------------------------------------------------------------------

[----------------------------------- INTERNAL LINKS -----------------------------------]: #

[#1]: #1-registry-maintainers

[#2]: #2-package-maintainers
[#2.1]: #21-setting-up
[#2.2]: #22-registering-a-new-package
[#2.3]: #23-releasing-a-new-version-of-a-package

[#3]: #3-package-users

[#4]: #4-references

[------------------------------------- REFERENCES -------------------------------------]: #

[local-registry]: https://github.com/GunnarFarneback/LocalRegistry.jl
[local-registry-migration]: https://github.com/GunnarFarneback/LocalRegistry.jl/blob/master/docs/migration_from_general.md
