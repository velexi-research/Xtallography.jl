Xtallography: Developer Notes
============================================================================================

## Known Issues

* When running in GitHub Actions, the Python unit tests occasionally and inconsistently
  fail because the tests crash (i.e., segmentation fault). These failures appear to be
  related to JuliaCall/PythonCall, but we have not yet tracked down the origin of the
  bug.

  The failed jobs usually pass when they are re-run.
