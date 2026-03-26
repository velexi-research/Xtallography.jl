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
Unit tests for `xtallography.symmetry.symmetry_elements` module
"""
# --- Imports

# Standard library
import unittest

# External packages
import juliacall
from xtallography.symmetry import GlidePlane, ScrewAxis

# Local packages/modules


# --- Test Suites


class test_xtallography_symmetry_symmetry_elements_GlidePlane(unittest.TestCase):
    """
    Test suite for the `GlidePlane` class
    """

    # --- Fixtures

    def setUp(self):
        """
        Prepare for test.
        """
        # Initialize JuliaCall
        self.jl = juliacall.newmodule("PyXtallographyTest")
        self.jl.seval("using Xtallography")

    def tearDown(self):
        """
        Clean up after test.
        """

    # --- Tests

    @staticmethod
    def test_init():
        """
        Test `__init__()`.
        """
        # --- Tests

        translation = "0,1,0"
        reflection_plane = "1,0,0"
        glide_plane = GlidePlane(translation, reflection_plane)

        assert glide_plane.translation == translation
        assert glide_plane.reflection_plane == reflection_plane

    def test_to_julia(self):
        """
        Test `to_julia()`.
        """
        # --- Preparations

        translation = "0,1,0"
        reflection_plane = "1,0,0"
        glide_plane = GlidePlane(translation, reflection_plane)

        # --- Tests

        glide_plane_jl = glide_plane.to_julia()
        assert self.jl.isa(glide_plane_jl, self.jl.GlidePlane)
        assert glide_plane_jl.translation == translation
        assert glide_plane_jl.reflection_plane == reflection_plane

    def test_from_julia(self):
        """
        Test `from_julia()`.
        """
        # --- Tests

        translation = "0,1,0"
        reflection_plane = "1,0,0"
        glide_plane_jl = self.jl.GlidePlane(translation, reflection_plane)

        glide_plane = GlidePlane.from_julia(glide_plane_jl)

        assert glide_plane.translation == translation
        assert glide_plane.reflection_plane == reflection_plane

    def test_eq(self):
        """
        Test `__eq__()`.
        """
        # --- Preparations

        translation = "0,1,0"
        reflection_plane = "1,0,0"
        glide_plane = GlidePlane(translation, reflection_plane)

        # --- Tests

        # ------ same glide planes

        other_glide_plane = GlidePlane(translation, reflection_plane)
        assert glide_plane == other_glide_plane

        # ------ different glide planes

        # translation differs
        other_glide_plane = GlidePlane("0,0,1", reflection_plane)
        assert glide_plane != other_glide_plane

        # reflection plane differs
        other_glide_plane = GlidePlane(translation, "0,0,1")
        assert glide_plane != other_glide_plane

        # ------ comparison with non-GlidePlane object

        other_glide_plane = ScrewAxis("1,0,0", 3, 2)
        assert glide_plane != other_glide_plane


class test_xtallography_symmetry_symmetry_elements_ScrewAxis(unittest.TestCase):
    """
    Test suite for the `ScrewAxis` class
    """

    # --- Fixtures

    def setUp(self):
        """
        Prepare for test.
        """
        # Initialize JuliaCall
        self.jl = juliacall.newmodule("PyXtallographyTest")
        self.jl.seval("using Xtallography")

    def tearDown(self):
        """
        Clean up after test.
        """

    # --- Tests

    @staticmethod
    def test_init():
        """
        Test `__init__()`.
        """
        # --- Tests

        axis = "1,0,0"
        n = 6
        m = 4
        screw_axis = ScrewAxis(axis, n, m)

        assert screw_axis.axis == axis
        assert screw_axis.n == n
        assert screw_axis.m == m

    def test_to_julia(self):
        """
        Test `to_julia()`.
        """
        # --- Preparations

        axis = "1,0,0"
        n = 6
        m = 4
        screw_axis = ScrewAxis(axis, n, m)

        # --- Tests

        screw_axis_jl = screw_axis.to_julia()
        assert self.jl.isa(screw_axis_jl, self.jl.ScrewAxis)
        assert screw_axis_jl.axis == axis
        assert screw_axis_jl.n == n
        assert screw_axis_jl.m == m

    def test_from_julia(self):
        """
        Test `from_julia()`.
        """
        # --- Tests

        axis = "1,0,0"
        n = 6
        m = 4
        screw_axis_jl = self.jl.ScrewAxis(axis, n, m)

        screw_axis = ScrewAxis.from_julia(screw_axis_jl)

        assert screw_axis.axis == axis
        assert screw_axis.n == n
        assert screw_axis.m == m

    def test_eq(self):
        """
        Test `__eq__()`.
        """
        # --- Preparations

        axis = "1,0,0"
        n = 6
        m = 4
        screw_axis = ScrewAxis(axis, n, m)

        # --- Tests

        # ------ same screw axes

        other_screw_axis = ScrewAxis(axis, n, m)
        assert screw_axis == other_screw_axis

        # ------ different screw axes

        # axis differs
        other_screw_axis = ScrewAxis("0,0,1", n, m)
        assert screw_axis != other_screw_axis

        # n differs
        other_screw_axis = ScrewAxis(axis, n + 1, m)
        assert screw_axis != other_screw_axis

        # m differs
        other_screw_axis = ScrewAxis(axis, n, m - 1)
        assert screw_axis != other_screw_axis

        # ------ comparison with non-ScrewAxis object

        # reflection plane differs
        other_screw_axis = GlidePlane("1,0,1", "0,0,1")
        assert screw_axis != other_screw_axis
