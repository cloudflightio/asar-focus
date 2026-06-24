/*
 * ENVISAT and ERS ASAR instrument focusser for QA4EO activity (c)
 * by CGI Estonia AS, Cloudflight Austria GmbH
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
#include "sar_chirp.h"

#include "cufft_plan.h"
#include "cuda_util.h"
#include "cuda_cleanup.h"
#include "cufft_checks.h"


namespace
{
void ApplyHammingWindow(std::vector<std::complex<float>>& data, double coeff) {
    const size_t N = data.size();
    // const double coeff = 0.75; // original 0.54; TODO: take from CON file
    for (size_t i = 0; i < N; i++) {
        size_t n = i;
        double term = (2 * M_PI * n) / N;
        double m = coeff - (1 - coeff) * cos(term);
        data[i] *= m;
    }
}

void ScaleWindowed(std::vector<std::complex<float>>& chirp_data)
{
    cufftHandle plan = {};
    CufftPlanCleanup fft_cleanup(plan);
    size_t fft_size = chirp_data.size();
    // waste to use cuda for 1 small FFT, however no CPU FFT library has been integrated to the project
    CHECK_CUFFT_ERR(cufftPlan1d(&plan, fft_size, CUFFT_C2C, 1));
    cufftComplex* d_chirp = nullptr;
    size_t byte_sz = fft_size * sizeof(cufftComplex);
    CHECK_CUDA_ERR(cudaMalloc(&d_chirp, byte_sz));
    CudaMallocCleanup mem_cleanup(d_chirp);
    CHECK_CUDA_ERR(cudaMemcpy(d_chirp, chirp_data.data(), byte_sz, cudaMemcpyHostToDevice));
    CHECK_CUFFT_ERR(cufftExecC2C(plan, d_chirp, d_chirp, CUFFT_FORWARD));

    std::vector<std::complex<float>> h_fft(fft_size);
    CHECK_CUDA_ERR(cudaMemcpy(h_fft.data(), d_chirp, byte_sz, cudaMemcpyDeviceToHost));

    // total power in frequency bins, TODO(priit) equivalent in time domain?
    double scale = 0.0;
    for (const auto& el : h_fft) {
        float mag = std::abs(el);
        scale += mag * mag;
    }

    scale /= fft_size;

    // scale time td chirp
    for (auto& el : chirp_data) {
        el /= sqrt(scale);
    }
}

void KrFromReplica(ChirpInfo& chirp) {
    if(chirp.replica_pulses.size() == 0)
        return;
    std::vector<double> idx;
    std::vector<double> phase;
    double phi0 = 0.0;
    double i0 = 0.0;
    double dphi = 0.0;
    float maxrep = 0.0f;
    for(auto& rep: chirp.replica_pulses) {
        maxrep = std::max(maxrep, abs(rep));
    }
    int crows=0, rowbeg=0, rowend=0;
    for(int row=0; row < 24; row++) {
        bool ok = true;
        for(int j=0; j < 36; j++) {
            float absrep = abs(chirp.replica_pulses[36*row+j]);
            ok &= absrep == 0.0f || 2*absrep > maxrep;
        }
        if(!ok) {
            if(crows > rowend-rowbeg) {
                rowbeg = row-crows;
                rowend = row;
            }
            crows=0;
        } else {
            crows++;
        }
    }
    if(rowend-rowbeg < 12)
        return;
    for(int i=36*rowbeg; i < 36*rowend; i++) {
        auto rep = chirp.replica_pulses[i];
        if(std::abs(rep) == 0)
            continue;
        double phi = std::arg(rep);
        if(phase.size() > 0) {
            double phi1 = phase.size() > 1 ? phi0 + (i-i0) * dphi : phi0;
            while(phi-phi1 < -M_PI)
                phi += 2*M_PI;
            while(phi-phi1 > +M_PI)
                phi -= 2*M_PI;
            dphi = (phi-phi0)/(i-i0);
        }
        i0 = i;
        phi0 = phi;
        idx.push_back(i0);
        phase.push_back(phi0);
    }
    std::vector<double> poly = Polyfit(idx, phase, 2);
    chirp.Kr = -poly[0] / M_PI * (chirp.range_sampling_rate * chirp.range_sampling_rate);
}

}

void GenerateChirpData(ChirpInfo& chirp, size_t padding_size) {
    KrFromReplica(chirp);

    double Kr = chirp.Kr;

    const int64_t n_samples = chirp.n_samples;
    const double dt = 1.0 / chirp.range_sampling_rate;

    chirp.time_domain.clear();
    chirp.time_domain.resize(n_samples);

    for (int64_t i = 0; i < n_samples; i++) {
        const double t = (i - n_samples / 2 + 0.5) * dt;
        const double phase = M_PI * Kr * t * t;
        std::complex<float> iq(cos(phase), sin(phase));

        chirp.time_domain[i] = iq;
    }

    chirp.padded_windowed_data = chirp.time_domain;

    if(chirp.apply_window)
    {
        ApplyHammingWindow(chirp.padded_windowed_data, chirp.alpha);
    }

    chirp.padded_windowed_data.resize(padding_size);

    // ApplyWindowPad(chirp_data, padding_size);
    if (chirp.apply_window) {
        ScaleWindowed(chirp.padded_windowed_data);
    } else {
        for (auto& el : chirp.padded_windowed_data) {
            el /= sqrt(n_samples);
        }
    }

}