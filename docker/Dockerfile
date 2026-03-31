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

FROM nvidia/cuda:12.4.1-devel-ubuntu22.04

ARG USER_ID=1001
ARG GROUP_ID=1001

RUN groupadd -g ${GROUP_ID} asar-focus && \
    useradd -m -s /bin/bash -u ${USER_ID} -g asar-focus asar-focus

# for tzdata installation - https://askubuntu.com/questions/909277/avoiding-user-interaction-with-tzdata-when-installing-certbot-in-a-docker-contai
ENV DEBIAN_FRONTEND=noninteractive
RUN rm /etc/apt/sources.list.d/cuda-ubuntu2204-x86_64.list
RUN apt update && apt install -y --no-install-recommends \
	gdal-bin libgdal-dev libboost-program-options-dev libboost-date-time-dev libboost-iostreams-dev libboost-log-dev zlib1g-dev

RUN apt install -y cmake git
RUN apt install -y clang clang-format clang-tidy libeigen3-dev unzip

USER asar-focus

COPY --chown=asar-focus:asar-focus . /asar-focus
WORKDIR /asar-focus

RUN cmake -DCMAKE_CUDA_ARCHITECTURES=75 -B build .
WORKDIR build
RUN make VERBOSE=1

USER root
RUN ln -s /asar-focus/build/asar_focus /usr/local/bin
