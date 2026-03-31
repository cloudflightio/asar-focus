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

#include "ers_im_parse.h"

#include <stdexcept>
#include <string>

#include "fmt/format.h"

#include "alus_log.h"
#include "cuda_algorithm.h"
#include "cuda_util.h"
#include "envisat_format_kernels.h"
#include "envisat_utils.h"
#include "ers_env_format.h"
#include "ers_sbt.h"
#include "parse_util.h"

#define DEBUG_PACKETS 0
#define HANDLE_DIAGNOSTICS 1

namespace {

template <class T>
[[nodiscard]] const uint8_t* CopyBSwapPOD(T& dest, const uint8_t* src) {
    memcpy(&dest, src, sizeof(T));
    dest = bswap(dest);
    return src + sizeof(T);
}

struct EchoMeta {
    mjd isp_sensing_time;

    uint16_t isp_length;  // Part of the FEP annotation, but this is the only useful field, others are static fillers.

    uint32_t data_record_number;
    // IDHT header
    uint8_t packet_counter;
    uint8_t subcommutation_counter;

    // Auxiliary data and replica/calibration pulses
    uint64_t onboard_time;
    uint8_t activity_task;
    uint8_t sample_flags;
    uint32_t image_format_counter;
    uint16_t pri_code;
    uint16_t swst_code;
};

#if DEBUG_PACKETS

struct ErsFepAndPacketMetadata {
    EchoMeta echo_meta;
};

[[maybe_unused]] void DebugErsPackets(const std::vector<ErsFepAndPacketMetadata> dm, size_t limit,
                                      std::ostream& out_stream) {
    out_stream << "isp_sense,isp_sense_micros,isp_length,data_rec_no,packet_counter,subcomm_counter,onboard_time,"
               << "activity_task,image_format_counter,pri_code,swst_code" << std::endl;
    for (size_t i{}; i < limit; i++) {
        const auto& dmi = dm.at(i).echo_meta;

        const auto isp_sensing_time_ptime = MjdToPtime(dmi.isp_sensing_time);
        // clang-format off
        out_stream << boost::posix_time::to_simple_string(isp_sensing_time_ptime) << ","
                   << isp_sensing_time_ptime.time_of_day().total_microseconds() << ","
                   << dmi.isp_length << "," << dmi.data_record_number << "," << (uint32_t)dmi.packet_counter  << ","
                   << (uint32_t)dmi.subcommutation_counter << "," << dmi.onboard_time << ","
                   << dmi.activity_task << "," << dmi.image_format_counter << "," << dmi.pri_code << ","
                   << dmi.swst_code << std::endl;
        // clang-format on
    }
}

#endif

inline void FetchUint32(const uint8_t* array_start, uint32_t& var) {
    [[maybe_unused]] const auto discard = CopyBSwapPOD(var, array_start);
}

#if HANDLE_DIAGNOSTICS
constexpr auto DR_NO_MAX{alus::asar::envformat::parseutil::MaxValueForBits<uint32_t, 32>()};
constexpr auto PACKET_COUNTER_MAX{alus::asar::envformat::parseutil::MaxValueForBits<uint8_t, 8>()};
constexpr uint8_t SUBCOMMUTATION_COUNTER_MAX{48};
constexpr auto IMAGE_FORMAT_COUNTER_MAX{alus::asar::envformat::parseutil::MaxValueForBits<uint32_t, 24>()};
constexpr auto SBT_MAX{alus::asar::envformat::parseutil::MaxValueForBits<uint32_t, 32>()};

void InitializeDataRecordNo(const uint8_t* packet_start, uint32_t& dr_no) {
    FetchUint32(packet_start + alus::asar::envformat::ers::highrate::PROCESSOR_ANNOTATION_ISP_SIZE_BYTES +
                    alus::asar::envformat::ers::highrate::FEP_ANNOTATION_SIZE_BYTES,
                dr_no);
    if (dr_no != 0) {
        dr_no--;
    } else {
        dr_no = DR_NO_MAX;
    }
}

void FetchSatelliteBinaryTime(const uint8_t* packet_start, uint32_t& sbt) {
    /*
         * The update of that binary counter shall occur every 4th PRI.
         * The time relation to the echo data in the format is as following:
         * the transfer of the ICU on-board time to the auxiliary memory
         * occurs t2 before the RF transmit pulse as depicted in
         * Fig. 4.4.2.4.6-3. The least significant bit is equal to 1/256 sec.
         *
         * ER-IS-ESA-GS-0002  - pg. 4.4 - 24
     */
    FetchUint32(packet_start + alus::asar::envformat::ers::highrate::SBT_OFFSET_BYTES, sbt);
}

void InitializeCounters(const uint8_t* packets_start, uint32_t& dr_no, uint8_t& packet_counter,
                        uint8_t& subcommutation_counter, uint32_t& image_format_counter) {
    InitializeDataRecordNo(packets_start, dr_no);
    // 36 for packet counter
    packet_counter = packets_start[alus::asar::envformat::ers::highrate::PROCESSOR_ANNOTATION_ISP_SIZE_BYTES +
                                   alus::asar::envformat::ers::highrate::FEP_ANNOTATION_SIZE_BYTES +
                                   alus::asar::envformat::ers::highrate::DATA_RECORD_NUMBER_SIZE_BYTES];
    if (packet_counter == 0) {
        packet_counter = PACKET_COUNTER_MAX;
    } else {
        packet_counter--;
    }

    // 37 for subcommutation counter.
    subcommutation_counter = packets_start[alus::asar::envformat::ers::highrate::PROCESSOR_ANNOTATION_ISP_SIZE_BYTES +
                                           alus::asar::envformat::ers::highrate::FEP_ANNOTATION_SIZE_BYTES +
                                           alus::asar::envformat::ers::highrate::DATA_RECORD_NUMBER_SIZE_BYTES + 1];
    if (subcommutation_counter == 0) {
        subcommutation_counter = SUBCOMMUTATION_COUNTER_MAX;
    } else {
        subcommutation_counter--;
    };

    FetchUint32(packets_start + alus::asar::envformat::ers::highrate::PROCESSOR_ANNOTATION_ISP_SIZE_BYTES +
                    alus::asar::envformat::ers::highrate::FEP_ANNOTATION_SIZE_BYTES +
                    alus::asar::envformat::ers::highrate::DATA_RECORD_NUMBER_SIZE_BYTES + 18,
                image_format_counter);
    if (image_format_counter == 0) {
        image_format_counter = IMAGE_FORMAT_COUNTER_MAX;
    } else {
        image_format_counter--;
    }
}

// oos - out of sequence
struct PacketParseDiagnostics {
    size_t dr_no_oos;            // No. of times gap happened.
    size_t dr_no_total_missing;  // Total sum of all gaps' value.
    size_t packet_counter_oos;
    size_t packer_counter_total_missing;
    size_t packet_counter_changes;
    size_t subcommutation_counter_oos;
    size_t subcommutation_counter_total_missing;
    size_t image_format_counter_oos;
    size_t image_format_counter_total_missing;
};

#endif

inline void FetchPriCode(const uint8_t* start_array, uint16_t& var) {
    var |= static_cast<uint16_t>(start_array[0]) << 8;
    var |= start_array[1];
}

// This is deprecated, will be removed along with the initial parsing flow.
inline void CalculatePrf(uint16_t pri_code, double& pri, double& prf) {
    // ISCE2 uses 210.943006e-9, but
    // ERS-1-Satellite-to-Ground-Segment-Interface-Specification_01.09.1993_v2D_ER-IS-ESA-GS-0001_EECF05
    // note 5 in 4.4:28 specifies 210.94 ns, which draws more similar results to PF-ERS results.
    pri = (pri_code + 2.0) * 210.94e-9;
    prf = 1 / pri;

    if (prf < 1640 || prf > 1720) {
        LOGW << "PRF value '" << prf << "' is out of the specification (ERS-1-Satellite-to-Ground-Segment-Interface-"
             << "Specification_01.09.1993_v2D_ER-IS-ESA-GS-0001_EECF05) range [1640, 1720] Hz" << std::endl;
    }
}

inline const uint8_t* FetchPacketFrom(const uint8_t* start, alus::asar::envformat::ForecastMeta& meta, uint32_t& dr_no,
                                      uint32_t& sbt) {
    // https://earth.esa.int/eogateway/documents/20142/37627/ERS-products-specification-with-Envisat-format.pdf
    // https://asf.alaska.edu/wp-content/uploads/2019/03/ers_ceos.pdf and ER-IS-ESA-GS-0002 mixed between three.
    auto it = CopyBSwapPOD(meta.isp_sensing_time, start);
    it += 12;  // No MJD for ERS in FEP.
    uint16_t isp_length;
    it = CopyBSwapPOD(isp_length, it);
    InitializeDataRecordNo(start, dr_no);
    FetchSatelliteBinaryTime(start, sbt);

    // ISP size - 1.
    if (isp_length != 11465) {
        throw std::runtime_error("FEP annotation's ISP length shall be 11465. Check the dataset, current one has " +
                                 std::to_string(isp_length));
    }

    return start + alus::asar::envformat::ers::highrate::MDSR_SIZE_BYTES;
}

inline void TransferMetadata(const EchoMeta& echo_meta, alus::asar::envformat::CommonPacketMetadata& common_meta) {
    common_meta.sensing_time = echo_meta.isp_sensing_time;
    common_meta.onboard_time = echo_meta.onboard_time;
    common_meta.pri_code = echo_meta.pri_code;
    common_meta.swst_code = echo_meta.swst_code;
    common_meta.sample_count = alus::asar::envformat::ers::highrate::MEASUREMENT_DATA_SAMPLE_COUNT;
    common_meta.measurement_array_length_bytes = alus::asar::envformat::ers::highrate::MEASUREMENT_DATA_SIZE_BYTES;
}

}  // namespace

namespace alus::asar::envformat {

std::vector<ForecastMeta> FetchErsL0ImForecastMeta(const std::vector<char>& file_data, const DSD_lvl0& mdsr,
                                                   boost::posix_time::ptime packets_start_filter,
                                                   boost::posix_time::ptime packets_stop_filter,
                                                   size_t& packets_before_start, size_t& packets_after_stop) {
    if (mdsr.dsr_size != ers::highrate::MDSR_SIZE_BYTES) {
        throw std::runtime_error("Expected DSR to be " + std::to_string(ers::highrate::MDSR_SIZE_BYTES) +
                                 " bytes, actual - " + std::to_string(mdsr.dsr_size) +
                                 ". This is a failsafe abort, check the dataset.");
    }
    if (file_data.size() - mdsr.ds_offset != mdsr.num_dsr * mdsr.dsr_size) {
        throw std::runtime_error("Actual packets BLOB array (" + std::to_string(file_data.size() - mdsr.ds_offset) +
                                 " bytes) is not equal to dataset defined packets (" +
                                 std::to_string(mdsr.num_dsr * mdsr.dsr_size) + " bytes).");
    }

    std::vector<ForecastMeta> forecast_meta{};
    size_t packets_before_requested_start{};
    size_t packets_after_requested_stop{};
    const auto* const packets_start = reinterpret_cast<const uint8_t*>(file_data.data() + mdsr.ds_offset);
    const uint8_t* next_packet = reinterpret_cast<const uint8_t*>(file_data.data() + mdsr.ds_offset);
    const auto start_mjd = PtimeToMjd(packets_start_filter);
    const auto stop_mjd = PtimeToMjd(packets_stop_filter);

    size_t dr_no_oos{};
    size_t dr_no_total_missing{};
    uint32_t last_dr_no;
    InitializeDataRecordNo(next_packet, last_dr_no);
    uint32_t last_sbt;
    uint32_t sbt_repeat{0};
    FetchSatelliteBinaryTime(next_packet, last_sbt);

    size_t i{};
    for (; i < mdsr.num_dsr; i++) {
        ForecastMeta current_meta;
        uint32_t dr_no;
        current_meta.packet_start_offset_bytes = next_packet - packets_start;
        uint32_t sbt;
        next_packet = FetchPacketFrom(next_packet, current_meta, dr_no, sbt);
        const auto sbt_delta = parseutil::CounterGap<uint32_t, SBT_MAX>(last_sbt, sbt);
        if (sbt_delta == 0) {
            sbt_repeat++;
        } else {
            sbt_repeat = 0;
        }
        bool overflow = (sbt < last_sbt ? true : false);
        ers::AdjustIspSensingTime(current_meta.isp_sensing_time, sbt, sbt_repeat, overflow);

        const auto dr_no_gap = parseutil::CounterGap<uint32_t, DR_NO_MAX>(last_dr_no, dr_no);
        // When gap is 0, this is not handled - could be massive rollover or duplicate or the first packet.
        if (dr_no_gap > 1) {
            if (forecast_meta.size() < 1) {
                throw std::runtime_error(
                    "Implementation error detected - data record no. discontinuity detected, but no packets have been "
                    "fetched yet in order to fill packets");
            }

            const auto packets_to_fill = dr_no_gap - 1;
            dr_no_oos++;
            dr_no_total_missing += packets_to_fill;
            // 900 taken from
            // https://github.com/gmtsar/gmtsar/blob/master/preproc/ERS_preproc/ers_line_fixer/ers_line_fixer.c#L434
            // For ISCE2 this is 30000 which seems to be wayyyy too biig.
            if (packets_to_fill > 900) {
                throw std::runtime_error(
                    fmt::format("There are {} packets missing, which is more than threshold 900 in order to continue. "
                                "It was detected between {}th and {}th packets.",
                                packets_to_fill, i - 1, i));
            }
            const auto& previous_meta = forecast_meta.back();
            // TODO - sensing times should be recalculated as well.
            for (size_t pf{}; pf < packets_to_fill; pf++) {
                forecast_meta.push_back(previous_meta);
            }

            if (previous_meta.isp_sensing_time < start_mjd) {
                packets_before_requested_start += packets_to_fill;
            }

            if (previous_meta.isp_sensing_time > stop_mjd) {
                packets_after_requested_stop += packets_to_fill;
            }
        } else {
            if (current_meta.isp_sensing_time < start_mjd) {
                packets_before_requested_start++;
            }

            if (current_meta.isp_sensing_time > stop_mjd) {
                packets_after_requested_stop++;
            }
        }

        // In case the packets_after_stop is equal to 0. Delete one at the back.
        if (packets_after_requested_stop > packets_after_stop) {
            packets_after_requested_stop--;
            break;
        }

        last_dr_no = dr_no;
        forecast_meta.push_back(current_meta);
    }

    if (packets_before_requested_start > packets_before_start) {
        forecast_meta.erase(forecast_meta.begin(),
                            forecast_meta.begin() + (packets_before_requested_start - packets_before_start));
        packets_before_requested_start -= (packets_before_requested_start - packets_before_start);
        if (packets_before_start == 0) {
            auto it = forecast_meta.begin();
            while (it != forecast_meta.end()) {
                const auto sensing_time = it->isp_sensing_time;
                if (sensing_time < start_mjd) {
                    forecast_meta.erase(it);
                    it = forecast_meta.begin();
                } else {
                    break;
                }
            }
        }
    }

    packets_before_start = packets_before_requested_start;
    packets_after_stop = packets_after_requested_stop;

    LOGD << "Data record number discontinuities/total missing - " << dr_no_oos << "/" << dr_no_total_missing
         << " these were inserted as to be parsed/duplicated during prefetch";

    return forecast_meta;
}

RawSampleMeasurements ParseErsLevel0ImPackets(const std::vector<char>& file_data, size_t mdsr_offset_bytes,
                                              const std::vector<ForecastMeta>& entries_to_be_parsed, InstrumentFile&,
                                              SARMetadata& sar_meta,
                                              std::vector<CommonPacketMetadata>& common_metadata) {
    const uint8_t* const packets_start = reinterpret_cast<const uint8_t*>(file_data.data()) + mdsr_offset_bytes;
    const uint8_t* it = packets_start + entries_to_be_parsed.front().packet_start_offset_bytes;
    const auto num_dsr = entries_to_be_parsed.size();

    // std::vector is chosen IF actual packets should exceed the num_dsr.
    const auto alloc_bytes_for_raw_samples = num_dsr * ers::highrate::MEASUREMENT_DATA_SIZE_BYTES;
    LOGD << "Reserving " << alloc_bytes_for_raw_samples / (1 << 20) << "MiB for raw samples ("
         << ers::highrate::MEASUREMENT_DATA_SIZE_BYTES << "x" << num_dsr << ")";

    RawSampleMeasurements raw_measurements{};
    raw_measurements.max_samples = ers::highrate::MEASUREMENT_DATA_SAMPLE_COUNT;
    raw_measurements.total_samples = ers::highrate::MEASUREMENT_DATA_SAMPLE_COUNT * num_dsr;
    raw_measurements.single_entry_length = ers::highrate::MEASUREMENT_DATA_SIZE_BYTES;
    raw_measurements.entries_total = num_dsr;
    raw_measurements.raw_samples = std::unique_ptr<uint8_t[]>(new uint8_t[alloc_bytes_for_raw_samples]);
    auto echoes_raw = raw_measurements.raw_samples.get();

#if DEBUG_PACKETS
    std::vector<ErsFepAndPacketMetadata> ers_dbg_meta;
    std::ofstream debug_stream("/tmp/" + asar_meta.product_name + ".debug.csv");
#endif
    // Data rec no. is used to detect missing lines as is in ISCE2. Image format counter seems to be suitable as well.
    uint32_t last_dr_no{};
    uint8_t last_packet_counter{};
    uint8_t last_subcommutation_counter{};
    uint32_t last_image_format_counter{};
    InitializeCounters(it, last_dr_no, last_packet_counter, last_subcommutation_counter, last_image_format_counter);
    uint32_t last_sbt{};
    FetchSatelliteBinaryTime(it, last_sbt);
#if HANDLE_DIAGNOSTICS
    PacketParseDiagnostics diagnostics{};
#endif

    struct {
        int count = 0;
        float I = 0.0f;
        float Q = 0.0f;
    } repsample[2][24*36];
    int reprow = 24;

    size_t packet_index{};
    // https://earth.esa.int/eogateway/documents/20142/37627/ERS-products-specification-with-Envisat-format.pdf
    // https://asf.alaska.edu/wp-content/uploads/2019/03/ers_ceos.pdf and ER-IS-ESA-GS-0002 mixed between three.
    for (const auto& entry : entries_to_be_parsed) {
        it = packets_start + entry.packet_start_offset_bytes;
        EchoMeta echo_meta = {};

        uint32_t sbt;
        FetchSatelliteBinaryTime(it, sbt);
        it = CopyBSwapPOD(echo_meta.isp_sensing_time, it);
        bool overflow = (sbt < last_sbt ? true : false);
        // sbt_repeat is not handled, because no mechanism could be developed.
        // SBT stays constant for some period of packets due to being updated after every 4th PRI.
        // When dividing SBT update time period with PRI no sensible period did develop between packets,
        // hence this is not used right now. There would be minimal range lines skipped/included due to this anyway.
        ers::AdjustIspSensingTime(echo_meta.isp_sensing_time, sbt, 0, overflow);
        last_sbt = sbt;

        if (echo_meta.isp_sensing_time != entry.isp_sensing_time) {
            throw std::runtime_error(fmt::format(
                "Discrepancy between packet list '{}' and actual parsed packet sensing time '{}' at packet index {}",
                to_simple_string(MjdToPtime(entry.isp_sensing_time)),
                to_simple_string(MjdToPtime(echo_meta.isp_sensing_time)), packet_index));
        }

        // FEP annotation only consists of useful ISP length for ERS.
        it += 12;  // No MJD for ERS.
        it = CopyBSwapPOD(echo_meta.isp_length, it);
        it += 6;  // No CRC or correction count for ERS.

        // ISP size - 1.
        if (echo_meta.isp_length != 11465) {
            throw std::runtime_error("FEP annotation's ISP length shall be 11465. Check the dataset, current one has " +
                                     std::to_string(echo_meta.isp_length) + " for packet no " +
                                     std::to_string(packet_index + 1));
        }
        it = CopyBSwapPOD(echo_meta.data_record_number, it);

        // ER-IS-ESA-GS-0002 4.4.2.4.5
        /*
         * It is a binary counter which is incremented each time a new source packet of general IDHT header data is
         * transmitted (about every 1 second). It is reset to zero with power-on of the IDHT-DCU.
         * The first format transmitted after power-on will show the value "1". Bit O is defined as MSB. bit 7 as LSB.
         */
        echo_meta.packet_counter = it[0];
        /*
         * binary counter counting the 8 byte segments from 1 to 48; it is reset for each new source packet.
         * Bit o is defined as MSB, bit 7 as LSB.
         */
        echo_meta.subcommutation_counter = it[1];
        it += 2;

        /*
         * These segments ot 8 bytes contain the subcommutated IDHT General Header source packet.
         */
        it += 8;
        // Table 4.4.2.6-2 in ER-IS-ESA-GS-0002
        if (it[0] != 0xAA) {
            throw std::runtime_error(fmt::format(
                "ERS data packet's auxiliary section shall start with 0xAA, instead {:#x} was found for packet no {}",
                it[0], packet_index + 1));
        }
        it += 1;
        // OBRC - On-Board Range Compressed | OGRC - On-Ground Range Compressed - See ER-IS-EPO-GS-0201
        it += 1;  // bit 0 - OGRC/OBRC flag (0/1 respectively). Bits 1-4 Orbit number (0-15 -> 1-16)
        /*
         * The update of that binary counter shall occur every 4th PRI.
         * The time relation to the echo data in the format is as following:
         * the transfer of the ICU on-board time to the auxiliary memory
         * occurs t2 before the RF transmit pulse as depicted in
         * Fig. 4.4.2.4.6-3. The least significant bit is equal to 1/256 sec.
         *
         * ER-IS-ESA-GS-0002  - pg. 4.4 - 24
         */
        echo_meta.onboard_time = sbt;
        it += 4;
        /*
         * This word shall define the activity task within the mode of
         * operation, accompanied by the validity flag bits and the first
         * sample flag bits. The update of the activity task word is controlled by the 4th PRI interrupt.
         * ER-IS-ESA-GS-0002  - pg. 4.4 - 25
         *
         * Pg. 4.4-25
         * Activity means generation of noise, echo, calibration or replica data.
         *
         * MSB  LSB
         * 10001000 - Noise; no calibration
         * 10011001 - No echo; Cal. drift (EM only)
         * 10101001 - Echo; Cal. drift
         * 10101010 - Echo; Replica
         * 10100000 - Echo; no Replica (because of OBRC)
         */
        echo_meta.activity_task = it[0];
        /*
         * ...continued from activity task's pg. 4.4-25
         * Bit
         * 0 - echo invalid/valid (0/1)
         * 1 - Cal. data/Replica data invalid/valid (0/1)
         * 2 - Is noise
         * 3 - Is Cal/Repl.
         * 4 - Is echo
         * ... spare bits
         *
         * The bits 2 to 4 (noise/cal/repl/echo) of byte 23 can only be set or reset with a 4 x PRI interval,
         * therefore the start of a new activity is flagged for the first four formats.
         */
        echo_meta.sample_flags = it[1];
        it += 2;
        // Updated every format. Reset at the beginning of a transmission.
        it = CopyBSwapPOD(echo_meta.image_format_counter, it);
        it = CopyBSwapPOD(echo_meta.swst_code, it);
        it = CopyBSwapPOD(echo_meta.pri_code, it);
        // Next item - Cal. S/S attenuation select - ER-IS-ESA-GS-0002 Table 4.4.2.4.6-2
        // Drift calibration data
        it += 122;
        uint16_t repencoded;
        for(int i=0; i<36; i++) {
            it = CopyBSwapPOD(repencoded, it);
            if(reprow >= 24) {
                continue;
            }
            auto &rep = repsample[1][36*reprow + i];
            if(( repencoded >> 12) == 15) {
                rep = {1, float(repencoded & 0x3F), float((repencoded >> 6) & 0x3F)};
            } else {
                rep = {0, 0.0f, 0.0f};
            }
        }
        reprow++;
        if((repencoded >> 12) != 15) {
            if(reprow == 24) {
                for(int i=0; i<24*36; i++) {
                    auto& rep0 = repsample[0][i];
                    auto& rep1 = repsample[1][i];
                    rep0.count += rep1.count;
                    rep0.I += rep1.I;
                    rep0.Q += rep1.Q;
                }
            }
            reprow = 0;
        }
        // Lets skip rest
        it += 10;
      
        const auto dr_no_gap = parseutil::CounterGap<uint32_t, DR_NO_MAX>(last_dr_no, echo_meta.data_record_number);
        // When gap is 0, it has been corrected by the previously fetched metadata.
        if (dr_no_gap == 0) {
            diagnostics.dr_no_oos++;
            diagnostics.dr_no_total_missing++;
        }
        last_dr_no = echo_meta.data_record_number;
#if HANDLE_DIAGNOSTICS
        const auto packet_counter_gap =
            parseutil::CounterGap<uint8_t, PACKET_COUNTER_MAX>(last_packet_counter, echo_meta.packet_counter);
        // Packet counter increments something like after 1679 packets?? based on IM L0 observation. When this happens
        // subcommutation counter also resets to 1.
        if (packet_counter_gap > 1) {
            diagnostics.packet_counter_oos++;
            diagnostics.packer_counter_total_missing += (packet_counter_gap - 1);
        }
        if (packet_counter_gap == 1) {
            diagnostics.packet_counter_changes++;
        }
        last_packet_counter = echo_meta.packet_counter;

        const auto subcommutation_gap = parseutil::CounterGap<uint8_t, SUBCOMMUTATION_COUNTER_MAX>(
            last_subcommutation_counter, echo_meta.subcommutation_counter);
        // When gap is 0, this is not handled - could be massive rollover or duplicate.
        // Also the counter does not start from 0.
        if (subcommutation_gap > 1 && echo_meta.subcommutation_counter != 1 && packet_counter_gap != 1) {
            diagnostics.subcommutation_counter_oos++;
            diagnostics.subcommutation_counter_total_missing += (subcommutation_gap - 1);
        }
        last_subcommutation_counter = echo_meta.subcommutation_counter;

        const auto image_format_counter_gap = parseutil::CounterGap<uint32_t, IMAGE_FORMAT_COUNTER_MAX>(
            last_image_format_counter, echo_meta.image_format_counter);
        if (image_format_counter_gap > 1) {
            diagnostics.image_format_counter_oos++;
            diagnostics.image_format_counter_total_missing += (image_format_counter_gap - 1);
        }
        last_image_format_counter = echo_meta.image_format_counter;
#endif

        std::memcpy(echoes_raw + packet_index * raw_measurements.single_entry_length, it,
                    raw_measurements.single_entry_length);
        TransferMetadata(echo_meta, common_metadata.at(packet_index));

        packet_index++;
#if DEBUG_PACKETS
        ErsFepAndPacketMetadata meta{};
        meta.echo_meta = echo_meta;
        ers_dbg_meta.push_back(meta);
#endif
    }
    for (auto rep: repsample[0]) {
        if(rep.count) {
            sar_meta.chirp.replica_pulses.push_back({rep.I/rep.count-31, rep.Q/rep.count-31});
        } else {
            sar_meta.chirp.replica_pulses.push_back({0.0f, 0.0f});
        }
    }

#if DEBUG_PACKETS
    DebugErsPackets(ers_dbg_meta, ers_dbg_meta.size(), debug_stream);
#endif

    if (diagnostics.dr_no_oos > 0) {
        raw_measurements.no_of_product_errors_compensated = diagnostics.dr_no_oos;
    }

    LOGD << "Data record number discontinuities/total missing - " << diagnostics.dr_no_oos << "/"
         << diagnostics.dr_no_total_missing << " - were detected during parsing.";

#if HANDLE_DIAGNOSTICS
    LOGD << "Packet counter field diagnostics - " << diagnostics.packet_counter_oos << " "
         << diagnostics.packer_counter_total_missing << " changes " << diagnostics.packet_counter_changes;
    LOGD << "Subcommutation counter field diagnostics - " << diagnostics.subcommutation_counter_oos << " "
         << diagnostics.subcommutation_counter_total_missing;
    LOGD << "Image format counter diagnostics - " << diagnostics.image_format_counter_oos << " "
         << diagnostics.image_format_counter_total_missing;
#endif

    return raw_measurements;
}

}  // namespace alus::asar::envformat
