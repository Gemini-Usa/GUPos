#include <iostream>
#include <fstream>
#include <chrono>
#include <Eigen/Dense>
#include "Algebra.h"
#include "Geodesy.h"
#include "ImuData.h"
#include "InsMechanizer.h"
#include "GnssData.h"
#include "GnssInsFilter.h"

using namespace Utility;
using namespace std::chrono;

constexpr double max_gap{ 0.011 };
double std_pos[3]{ 0.009, 0.008, -0.022 };
double std_vel[3]{ 0.000, 0.000, 0.000 };
double std_att[3]{ D2R(0.05), D2R(0.05), D2R(0.05) };
constexpr double lever[3]{ 0.2350,-0.1000,0.8900 };
double ARW{ D2R(0.2) / 60.0 };
double VRW{ 0.4 / 60.0 };
double std_gb[3]{ D2R(24) / 3600.0, D2R(24) / 3600.0, D2R(24) / 3600.0 };
double std_ab[3]{ 400 * 1e-5, 400 * 1e-5, 400 * 1e-5 };
double std_gs[3]{ 1000 * 1e-6, 1000 * 1e-6, 1000 * 1e-6 };
double std_as[3]{ 1000 * 1e-6, 1000 * 1e-6, 1000 * 1e-6 };

void ReadIMUFile(const std::string &filename, std::vector<ImuData>& imu_data) {
    std::ifstream ifs{ filename, std::ios::in };
    if (!ifs.is_open()) return;
    std::string line;
    ImuData curr_imu, prev_imu;
    while (getline(ifs, line)) {
        prev_imu = curr_imu;
        if (!curr_imu.ParseASC(line)) continue;
        if (curr_imu.isDuplicated(prev_imu)) continue;
        double dt = curr_imu.getSecond() - prev_imu.getSecond();
        if (dt > max_gap && dt < 1.0)
            imu_data.push_back(ImuData::Interpolate(prev_imu, curr_imu, prev_imu.getSecond() + 0.01));
        imu_data.emplace_back(curr_imu);
    }
    ifs.close();
}

EulerAngle initialAlignment(const std::vector<ImuData>& imu_data, const double *initial_pos, double start_time = -1, double end_time = -1) {
    double tos = start_time < 0 ? imu_data.front().getSecond() : start_time;
    double toe = end_time < 0 ? imu_data.back().getSecond() : end_time;
    int epoch{ 0 };
    EulerAngle euler;
    ImuData avg_imu;
    for (auto& imu : imu_data) {
        if (imu.getSecond() < tos) continue;
        if (imu.getSecond() > toe) break;
        avg_imu += imu;
        ++epoch;
    }
    avg_imu /= epoch;
    avg_imu.StaticAlignment(initial_pos, euler);
    return euler;
}

void INSMechanize(const std::vector<ImuData>& imu_data,
                  const double *init_pos, const double *init_vel, const Eigen::Quaterniond &init_att,
                  double start_time = -1, double end_time = -1) {
    double tos = start_time < 0 ? imu_data.front().getSecond() : start_time;
    double toe = end_time < 0 ? imu_data.back().getSecond() : end_time;
    InsMechanizer ins;
    for (auto& imu : imu_data) {
        if (fabs(imu.getSecond() - tos) < 1.0e-3) {
            ins.Initialize(init_pos, init_vel, init_att, imu);
            continue;
        }
        if (imu.getSecond() < tos || imu.getSecond() > toe) continue;
        ins.INSUpdate(imu);
        ins.PrintState();
    }
}

void ReadPosFile(const std::string &filename, std::vector<GnssData>& pos_results) {
    std::ifstream ifs{ filename, std::ios::in };
    if (!ifs.is_open()) return;
    std::string line;
    GnssData pos_result;
    for (int i = 0; i < 2; ++i) getline(ifs, line);
    while (getline(ifs, line)) {
        pos_result.Parse(line);
        pos_results.push_back(pos_result);
    }
    ifs.close();
}

void LooselyCouple(const std::vector<ImuData> &imu_data, const std::vector<GnssData> &gnss_data, const double *init_pos) {
    GnssInsFilter filter;
    // Initial Alignment, get initial position, velocity and attitude
    double tos = imu_data.at(0).getSecond();
    double toe = imu_data.at(0).getSecond() + 300.0;
    auto init_euler = initialAlignment(imu_data, init_pos, tos, toe);
    double init_vel[3]{ 0.0, 0.0, 0.0 };
    auto init_att = EulerAngleToQuaternion(init_euler);
    const auto pos = std::make_pair(init_pos, std_pos);
    const auto vel = std::make_pair(init_vel, std_vel);
    const auto att = std::make_pair(init_att, std_att);
    double init_gb[3]{ 0.0, 0.0, 0.0 };
    double init_ab[3]{ 0.0, 0.0, 0.0 };
    double init_gs[3]{ 0.0, 0.0, 0.0 };
    double init_as[3]{ 0.0, 0.0, 0.0 };
    const auto gb = std::make_pair(init_gb, std_gb);
    const auto ab = std::make_pair(init_ab, std_ab);
    const auto gs = std::make_pair(init_gs, std_gs);
    const auto as = std::make_pair(init_as, std_as);
    // Loosely Couple by EKF
    filter.setLeverArm(lever);
    filter.setMarkovTime(3600.0, 3600.0, 3600.0, 3600.0);
    filter.setRandomWalk(ARW, VRW);
    // Time Synchronize
    double i_time, g_time;
    int i = 0, k = 0;
    while(fabs(imu_data.at(i).getSecond() - toe) > 1.0E-3) ++i;
    filter.Initialize(pos, vel, att, gb, ab, gs, as, imu_data.at(i));
    for (++i; i < imu_data.size(); ++i) {
        i_time = imu_data.at(i).getSecond();
        g_time = gnss_data.at(k).getSecond();
        if (k == gnss_data.size() - 1) break;
        while (g_time < i_time) g_time = gnss_data.at(++k).getSecond();
        if (fabs(g_time - i_time) < 1.0E-9) filter.ProcessData(imu_data.at(i), gnss_data.at(k++));
        else filter.ProcessData(imu_data.at(i));
    }
    filter.PrintState();
}

int main() {
    std::string public_path = "..//Rawdata//";
    std::string static_file = "3-2-static.ASC";
    std::string open_file = "3-2-open.ASC";
    double static_init_pos[3] = { D2R(30.5278045116), D2R(114.3557164314), 22.119 }; // static
    double open_init_pos[3]{ D2R(30.52844163), D2R(114.35697688), 22.334 }; // open
    std::vector<ImuData> imu_data;
    imu_data.reserve(150000);
    std::vector<GnssData> pos_results;
    ReadPosFile("..//Data//open.pos", pos_results);
    auto start = system_clock::now();
    ReadIMUFile(public_path + open_file, imu_data);
    LooselyCouple(imu_data, pos_results, open_init_pos);
    auto end = system_clock::now();
    auto duration = duration_cast<seconds>(end - start);
    std::cout << "duration: " << duration.count() << " seconds" << std::endl;
}
