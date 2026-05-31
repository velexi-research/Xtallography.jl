#   Copyright 2026 Velexi Corporation
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
The `xtallography.symmetry.symmetry_elements` package defines classes and methods related
symmetry.
"""
# Package-level types

from .symmetry_element import SymmetryElement  # noqa
from .rotation_axis import RotationAxis  # noqa
from .rotoinversion_axis import RotoinversionAxis  # noqa
from .glide_plane import GlidePlane  # noqa
from .screw_axis import ScrewAxis  # noqa

from .symmetry_element_dynamic import symmetry_element_from_julia


# --- Dynamically update SymmetryElement class


setattr(SymmetryElement, "from_julia", classmethod(symmetry_element_from_julia))
