#include "gain_correction.cuh"

#include <vector>
#include <string>
#include <cufft.h>
#include "cuda_util/cuda_cleanup.h"
#include "cuda_util/cuda_util.h"
#include "cuda_util/cufft_plan.h"
#include "cufft_checks.h"


#define GAIN_CORRECTION_GROUNDPOINTS 20
#define RSL_REF_RANGE 800e3

__global__ void GainCorrectionApply(cuComplex* data, const double* correction, int range_size, int azimuth_size, int correction_size) {
    const int x = threadIdx.x + blockDim.x * blockIdx.x;
    const int y = threadIdx.y + blockDim.y * blockIdx.y;
    const int data_idx = y * range_size + x;

    if (x >= correction_size) return;

    data[data_idx].x *= correction[x];
    data[data_idx].y *= correction[x];
}

double gain_corr_point(GeoPos3D spos, GeoPos3D gpos, double* eap_corr, double* angles) {
    GeoPos3D vpos = Difference(gpos, spos);
    double rsl = pow(Distance({0,0,0}, vpos) / RSL_REF_RANGE, 1.5);
    double angle = Degree(spos, vpos);
    double last_dist = 100, dist;
    for (int i = 0; i < 201; i++) {
        dist = abs(angle - angles[i]);
        if (dist > last_dist) {
            return eap_corr[i-1] * rsl;
        }
        last_dist = dist;
    }
    return eap_corr[200] * rsl;
}

void GainCorrection(DevicePaddedImage& img, SARMetadata &sar_meta, ASARMetadata &asar_meta, InstrumentFile &ins, alus::asar::envformat::aux::ExternalCalibration &xca) {

    // get eap correction
    const int swath_idx = alus::asar::envformat::SwathIdx(asar_meta.swath);
    const float angle = ins.fbp.ref_elev_angle[swath_idx];
    std::vector<std::string> pol_options = {"H/H","V/V","H/V","V/H"};
    auto it = std::find(pol_options.begin(), pol_options.end(), asar_meta.polarization.substr(0,3));
    int polaridx = it - pol_options.begin();

    float lut[201];
    switch(swath_idx) {
        case 0:
            memcpy(lut, xca.two_way_ant_elev_pat_gain_lut_is1 + polaridx * 201, 201 * sizeof(float));
            break;
        case 1:
            memcpy(lut, xca.two_way_ant_elev_pat_gain_lut_is2 + polaridx * 201, 201 * sizeof(float));
            break;
        case 2:
            memcpy(lut, xca.two_way_ant_elev_pat_gain_lut_is3_ss2 + polaridx * 201, 201 * sizeof(float));
            break;
        case 3:
            memcpy(lut, xca.two_way_ant_elev_pat_gain_lut_is4_ss3 + polaridx * 201, 201 * sizeof(float));
            break;
        case 4:
            memcpy(lut, xca.two_way_ant_elev_pat_gain_lut_is5_ss4 + polaridx * 201, 201 * sizeof(float));
            break;
        case 5:
            memcpy(lut, xca.two_way_ant_elev_pat_gain_lut_is6_ss5 + polaridx * 201, 201 * sizeof(float));
            break;
        case 6:
            memcpy(lut, xca.two_way_ant_elev_pat_gain_lut_is7 + polaridx * 201, 201 * sizeof(float));
            break;
        default:
            LOGE << "Invalid swath index for eap correction: " << std::to_string(swath_idx);
            break;
    }
    double eap_corr[201];
    double angles[201];
    for (int i = 0; i < 201; i++) {
        eap_corr[i] = pow(10, -lut[i]/20);
        angles[i] = angle - 5 + i * (10.0/201);
    }

    // get center orbit
    const boost::posix_time::ptime center_time = asar_meta.sensing_start + (asar_meta.sensing_stop-asar_meta.sensing_start) / 2;
    std::string center_time_str = boost::posix_time::to_simple_string(center_time);
    OrbitStateVector osv = InterpolateOrbit(sar_meta.osv, center_time);
    GeoPos3D spos = {osv.x_pos, osv.y_pos, osv.z_pos};
    Velocity3D svel = {osv.x_vel, osv.y_vel, osv.z_vel};

    // calc gain correction
    std::vector<double> gc;
    std::vector<double> gc_interp;
    int range_samp, last_range_samp = 0;
    double range_interval = ((float)sar_meta.img.range_size) / (GAIN_CORRECTION_GROUNDPOINTS - 1);
    for (int i = 0; i < GAIN_CORRECTION_GROUNDPOINTS; i++) {
        range_samp = round(i * range_interval);
        const int R0 = CalcSlantRange(sar_meta, range_samp);
        GeoPos3D ground_point = RangeDopplerGeoLocate(svel, spos, sar_meta.center_point, R0);
        gc.push_back(gain_corr_point(spos, ground_point, eap_corr, angles));

        if (i == 0) continue;
        int last_range_interval = range_samp - last_range_samp;
        for (float j = 0; j < last_range_interval; j+=1.0) {
            gc_interp.push_back((j / last_range_interval) * gc.rbegin()[0] + (1 - (j / last_range_interval)) * gc.rbegin()[1]);
        }
        last_range_samp = range_samp;
    }

    // prepare cuda memory
    cuComplex* d_data = img.Data();
    double* d_gc;
    const size_t gc_interp_size = sizeof(double) * gc_interp.size();
    CHECK_CUDA_ERR(cudaMalloc(&d_gc, gc_interp_size));
    CHECK_CUDA_ERR(cudaMemcpy(d_gc, gc_interp.data(), gc_interp_size, cudaMemcpyHostToDevice));

    int x_size = img.XStride();
    int y_size = img.YStride();

    // apply correction
    dim3 block_size(16, 16);
    dim3 grid_size((x_size + 15) / 16, (y_size + 15) / 16);
    GainCorrectionApply<<<grid_size, block_size>>>(d_data, d_gc, x_size, y_size, gc_interp.size());

    CHECK_CUDA_ERR(cudaDeviceSynchronize());
    CHECK_CUDA_ERR(cudaGetLastError());
}