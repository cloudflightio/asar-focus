# ENVISAT and ERS ASAR instrument focusser for QA4EO activity (c)
# by Cloudflight Austria GmbH
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

SCRIPTDIR=$(dirname $(readlink -f ${BASH_SOURCE[0]}))
cd "$SCRIPTDIR"
docker build .. -f Dockerfile -t asar-focus \
    --build-arg USER_ID=$(id -u) \
    --build-arg GROUP_ID=$(id -g) \
    --progress=plain 2>&1
) | tee build.log
