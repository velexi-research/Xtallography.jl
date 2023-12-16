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
Constants used throughout the PXRD module 
"""
# --- Exports

# Constants
export ALMOST_ZERO

# --- Constants

# Numerical comparison constants
const ALMOST_ZERO = 1e-9

const COS_APPROX_ZERO = 1e-9

# Common trigonometric computations
const SIN_PI_OVER_THREE = sqrt(3) / 2
const SIN_PI_OVER_FOUR = 1 / sqrt(2)
const ACOS_MINUS_ONE_THIRD = acos(-1 / 3)
