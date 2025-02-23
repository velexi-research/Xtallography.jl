#   Copyright 2023 Velexi Corporation
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

"""
Unit tests for Xtallography package.
"""

# --- Imports

# Standard library
using Test

# External packages
using Aqua: Aqua
using TestTools: jltest, jltest.EnhancedTestSet

# Xtallography package
using Xtallography

# --- Run tests

jltest.run_tests(@__DIR__)

# --- Run Aqua code quality checks

println()
print("Aqua.jl checks: ")
@testset EnhancedTestSet "Aqua.jl code quality checks" begin
    Aqua.test_all(Xtallography; deps_compat=(ignore=[:LinearAlgebra, :Logging, :Test],))
end
