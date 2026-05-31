#   Copyright 2025 Velexi Corporation
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
Unit tests for `SymmetryElement` class
"""
# --- Imports

# Standard library
import unittest

# External packages
import pytest
from xtallography.symmetry import SymmetryElement

# Local packages/modules


# --- Test Suites


class test_xtallography_symmetry_symmetry_elements_SymmetryElement(unittest.TestCase):
    """
    Test suite for the `SymmetryElement` class
    """

    def test_from_julia_invalid_args(self):
        """
        Test argument checks for `from_julia()`
        """
        symmetry_element_jl = 10
        with pytest.raises(ValueError) as exception_info:
            SymmetryElement.from_julia(symmetry_element_jl)

        expected_error = (
            "`symmetry_element_jl` must be a Julia `SymmetryElement` object. "
            f"(symmetry_element_jl={symmetry_element_jl})."
        )
        assert expected_error in str(exception_info)
