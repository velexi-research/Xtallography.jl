# Copyright (c) 2025 Velexi Corporation
#
# Use of this software is governed by the Business Source License included
# in the LICENSE file and at www.mariadb.com/bsl11.
#
# Change Date: Four years from the date the Licensed Work is published
#
# Change License: Mozilla Public License 2.0
"""
xtallography.lattices package

It provides implementations of lattice-specific types and methods.
"""
# --- Imports

# External packages
import juliacall


# --- Initialize JuliaCall

jl = juliacall.newmodule("PyXtallography")
jl.seval("using XtallographyUtils")
