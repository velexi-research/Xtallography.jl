# Copyright (c) 2025 Velexi Corporation
#
# Use of this software is governed by the Business Source License included
# in the LICENSE file and at www.mariadb.com/bsl11.
#
# Change Date: Four years from the date the Licensed Work is published
#
# Change License: Mozilla Public License 2.0
"""
The `xtallography` Python package implements a basic Python interface to the
<a href="/Xtallography.jl/">`Xtallography.jl`</a>.
"""
# --- Imports

# External packages
import juliacall


# --- Initialize JuliaCall

_JL = juliacall.newmodule("xtallographyJuliaCall")
_JL.seval("using Xtallography")
