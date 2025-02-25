#!/usr/bin/env bash
#------------------------------------------------------------------------------
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
#------------------------------------------------------------------------------

# Check that juliapkg.json is not the development version
SCRIPT_DIR="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
JULIAPKG_JSON_PATH="$SCRIPT_DIR/../xtallography/juliapkg.json"
RESULT=`grep '"dev":' $JULIAPKG_JSON_PATH | grep "true"`
if [ ! -z "$RESULT" ]; then
    echo "Development version of juliapkg.json should not be committed to the repository."
    exit 1
fi

exit 0
