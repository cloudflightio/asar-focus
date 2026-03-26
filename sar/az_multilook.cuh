/*
* ENVISAT and ERS ASAR instrument focusser for QA4EO activity (c)
* by CGI Estonia AS
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*/
#pragma once


#include "device_padded_image.cuh"
#include "sar_metadata.h"

class AzimuthMultiLook {
    float ratio;
    int in_az_size;
    int out_az_size;
    int range_size;

public:
    AzimuthMultiLook(const DevicePaddedImage& d_slc, SARMetadata& sar_meta);
    void apply(const DevicePaddedImage& d_slc, CudaWorkspace& out);
};
