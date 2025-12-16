#include "sar_metadata.h"
#include "envisat_lvl0_parser.h"
#include "orbit_state_vector.h"
#include "geo_tools.h"
#include "cuda_util/device_padded_image.cuh"


double gain_corr_point(GeoPos3D spos, GeoPos3D gpos, double* eap_corr);

void GainCorrection(DevicePaddedImage& img, SARMetadata &sar_meta, ASARMetadata &asar_meta, InstrumentFile &ins, alus::asar::envformat::aux::ExternalCalibration &xca);