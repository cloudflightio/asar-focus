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

#include "sar_metadata.h"

// void *PrepareSlantRangeToGroundRangeInterpolation(SARMetadata& sar_meta);

// void DoSlantRangeToGroundRangeInterpolation(void *pargs, const float* d_slant_in, float* d_gr_out);

class SlantRangeToGroundRangeInterpolation {
public:
    struct KernelArgs {
        float grsr_poly[5];
        int out_x_size;
        int in_x_size;
        int y_size;
        float slant_range_spacing;
        float ground_spacing;
    } args;
    SlantRangeToGroundRangeInterpolation(const SARMetadata& sar_meta);
    void apply(const float* d_slant_in, float* d_gr_out);
    void adjust_meta(SARMetadata& sar_meta);
};