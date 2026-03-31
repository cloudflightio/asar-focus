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

ARG CUDA_VERSION=12.4.1
ARG CUDA_BUILD_IMAGE=nvidia/cuda:$CUDA_VERSION-devel-ubuntu22.04
ARG CUDA_RUN_IMAGE=nvidia/cuda:$CUDA_VERSION-runtime-ubuntu22.04


FROM $CUDA_BUILD_IMAGE AS build

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
        cmake git clang clang-format clang-tidy \
        gdal-bin libgdal-dev libboost-program-options-dev libboost-date-time-dev libboost-iostreams-dev libboost-log-dev

COPY . /src
RUN mkdir /build

RUN cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CUDA_ARCHITECTURES=75 -DALUS_ENABLE_TESTS=ON -B /build /src

WORKDIR /build
RUN make -j12 VERBOSE=1


FROM build AS test
CMD ["find", "unit-test/", "-type", "f", "-exec", "{}", ";"] 


FROM $CUDA_RUN_IMAGE AS app

ARG BOOST_VERSION=1.74.0
ARG GDAL_VERSION=30
ARG USER_ID=1001
ARG GROUP_ID=1001

RUN groupadd -g ${GROUP_ID} asar-focus && \
    useradd -l -m -s /bin/bash -u ${USER_ID} -g asar-focus asar-focus

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
        libgdal$GDAL_VERSION \
        libboost-program-options$BOOST_VERSION libboost-date-time$BOOST_VERSION libboost-iostreams$BOOST_VERSION libboost-log$BOOST_VERSION

RUN mkdir /dist
WORKDIR /dist

USER asar-focus
COPY --from=build /build/asar_focus .

ENTRYPOINT ["./asar_focus"]