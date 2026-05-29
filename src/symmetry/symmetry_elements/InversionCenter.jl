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
InversionCenter type and functions
"""
# --- Exports

# Types
export InversionCenter

# --- Types

"""
    InversionCenter

Type representing an inversion through a point

Supertype: [`SymmetryElement`](@ref)
"""
struct InversionCenter <: SymmetryElement
    # --- Fields

    # location of inversion center
    center::Tuple{Rational,Rational,Rational}
end
