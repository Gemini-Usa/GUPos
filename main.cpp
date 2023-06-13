#include <iostream>
#include <iomanip>
#include <fstream>
#include <chrono>
#include <Eigen/Dense>
#include "Algebra.h"
#include "Geodesy.h"
#include "ImuData.h"
#include "InsMechanizer.h"
#include "GnssData.h"
#include "GnssInsFilter.h"
#include "GuConfig.h"
#include "String.h"

using namespace Utility;
using namespace std::chrono;

constexpr char conf_file[]{ "..//conf.ini" };
constexpr double max_gap{ 0.011 };

struct TrueResult {
    double second{ 0.0 };
    double blh[3]{ 0.0 };
    double vel[3]{ 0.0 };
    EulerAngle att{};
};

void ReadTrueResult(std::vector<TrueResult> &res, const std::string &file_name) {
    std::ifstream ifs{ file_name, std::ios::in };
    if (!ifs.is_open()) return;
    std::string line;
    std::vector<std::string> substr;
    TrueResult r;
    while (getline(ifs, line)) {
        RemoveSpace(line);
        SplitString(line, substr, " ");
        r.second = stod(substr[1]);
        r.blh[0] = D2R(stod(substr[2]));
        r.blh[1] = D2R(stod(substr[3]));
        r.blh[2] = stod(substr[4]);
        r.vel[0] = stod(substr[6]);
        r.vel[1] = stod(substr[5]);
        r.vel[2] = -stod(substr[7]);
        r.att.yaw = D2R(stod(substr[8]));
        r.att.roll = D2R(stod(substr[10]));
        r.att.pitch = D2R(stod(substr[9]));
        res.push_back(r);
        substr.clear();
    }
    ifs.close();
}

void ReadIMUFile(const std::string &filename, std::vector<ImuData>& imu_data) {
    if (filename.find(".ASC") != std::string::npos) {
        std::ifstream ifs{filename, std::ios::in};
        if (!ifs.is_open()) return;
        std::string line;
        ImuData curr_imu, prev_imu;
        while (getline(ifs, line)) {
            prev_imu = curr_imu;
            if (!curr_imu.ParseAsc(line)) continue;
            if (curr_imu.isDuplicated(prev_imu)) continue;
            double dt = curr_imu.getSecond() - prev_imu.getSecond();
            if (dt > max_gap && dt < 1.0)
                imu_data.push_back(ImuData::Interpolate(prev_imu, curr_imu, prev_imu.getSecond() + 0.01));
            imu_data.emplace_back(curr_imu);
        }
        ifs.close();
    } else if (filename.find(".imr") != std::string::npos) {
        std::ifstream ifs{ filename, std::ios::in | std::ios::binary };
        if (!ifs.is_open()) return;
        char header[512];
        ifs.read(header, 512);
        ImuData::ParseImrHeader(header);
        char buffer[32];
        ImuData curr_imu, prev_imu;
        while (ifs.read(buffer, 32)) {
            prev_imu = curr_imu;
            curr_imu.ParseImr(buffer);
            if (curr_imu.isDuplicated(prev_imu)) continue;
            double dt = curr_imu.getSecond() - prev_imu.getSecond();
            if (dt > max_gap && dt < 1.0)
                imu_data.push_back(ImuData::Interpolate(prev_imu, curr_imu, prev_imu.getSecond() + 0.01));
            imu_data.emplace_back(curr_imu);
        }
        ifs.close();
    } else {
        return;
    }
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

void LooselyCouple(const std::vector<ImuData> &imu_data, const std::vector<GnssData> &gnss_data, const std::vector<TrueResult> &true_results, const GuConfig &conf) {
    GnssInsFilter filter;
    std::ofstream ofs{ conf.output, std::ios::out };
    // Initial Alignment, get initial position, velocity and attitude
    double tos, toe;
    EulerAngle init_euler = conf.init_euler;
    if (conf.is_align) {
        tos = imu_data.at(0).getSecond();
        toe = tos + conf.align_time;
        init_euler = initialAlignment(imu_data, conf.init_pos, tos, toe);
    }
    else toe = imu_data.at(0).getSecond();
    auto init_att = EulerAngleToQuaternion(init_euler);
    const auto pos = std::make_pair(conf.init_pos, conf.init_std_pos);
    const auto vel = std::make_pair(conf.init_vel, conf.init_std_vel);
    const auto att = std::make_pair(init_att, conf.init_std_att);
    const auto gb = std::make_pair(conf.init_gb, conf.init_std_gb);
    const auto ab = std::make_pair(conf.init_ab, conf.init_std_ab);
    const auto gs = std::make_pair(conf.init_gs, conf.init_std_gs);
    const auto as = std::make_pair(conf.init_as, conf.init_std_as);
    // Loosely Couple by EKF
    filter.setLeverArm(conf.lever);
    filter.setMarkovTime(conf.markov_time, conf.markov_time, conf.markov_time, conf.markov_time);
    filter.setRandomWalk(conf.ARW, conf.VRW);
    // Time Synchronize
    double i_time, g_time, r_time;
    int i = 0, k = 0, n = 0;
    while(fabs(imu_data.at(i).getSecond() - toe) > 1.0E-3) ++i;
    while(gnss_data.at(k).getSecond() < toe) ++k;
    filter.Initialize(pos, vel, att, gb, ab, gs, as, imu_data.at(i));
    for (++i; i < imu_data.size(); ++i) {
        if (k == gnss_data.size() || n == true_results.size()) break;
        i_time = imu_data.at(i).getSecond();
        g_time = gnss_data.at(k).getSecond();
        r_time = true_results.at(n).second;
        if (g_time < i_time && fabs(g_time - i_time) > 1.0E-3) {
            ImuData inter_imu = ImuData::Interpolate(imu_data.at(i - 1), imu_data.at(i), g_time);
            filter.ProcessData(inter_imu, gnss_data.at(k++));
        } else if (fabs(g_time - i_time) < 1.0E-3) {
            filter.ProcessData(imu_data.at(i), gnss_data.at(k++));
        } else {
            filter.ProcessData(imu_data.at(i));
        }
        auto [p, v, a] = filter.getState();
        auto euler = QuaternionToEulerAngle(a);
        double ned[3];
        BlhToNed(p, conf.init_pos, ned);
        if (ofs.is_open() && fabs(r_time - i_time) < 1.0E-3) {
            const auto ref = true_results.at(n);
            double ref_ned[3];
            BlhToNed(ref.blh, conf.init_pos, ref_ned);
            double d_yaw = R2D(ref.att.yaw - euler.yaw);
            if (d_yaw < -350) d_yaw = M_PI * 2;
            else if (d_yaw > 350) d_yaw = -M_PI * 2;
            else d_yaw = 0;
            ofs << std::fixed << std::setprecision(3) << i_time << ","
                << std::setprecision(9) << ref_ned[0] - ned[0] << "," << ref_ned[1] - ned[1] << "," << ref_ned[2] - ned[2] << ","
                << std::setprecision(6) << ref.vel[0] - v[0] << "," << ref.vel[1] - v[1] << "," << ref.vel[2] - v[2] << ","
                << std::setprecision(3) << R2D(ref.att.roll - euler.roll) << "," << R2D(ref.att.pitch - euler.pitch) << "," << R2D(d_yaw + ref.att.yaw - euler.yaw)
                << ","
                << filter.isZeroUpdate() << "\n";
            ++n;
        }
    }
    filter.PrintState();
}

int main() {
    GuConfig conf;
    conf.ParseConfig(conf_file);
    std::vector<ImuData> imu_data;
    imu_data.reserve(150000);
    std::vector<GnssData> pos_results;
    ReadPosFile(conf.input_gnss, pos_results);
    std::vector<TrueResult> true_results;
    ReadTrueResult(true_results, "..//..//Data//true_open.pos");
    auto start = system_clock::now();                                 // timing
    ReadIMUFile(conf.input_imu, imu_data);
    LooselyCouple(imu_data, pos_results, true_results, conf);
    auto end = system_clock::now();                                   // timing
    auto duration = duration_cast<seconds>(end - start);    // timing
    std::cout << "duration: " << duration.count() << " seconds" << std::endl;       // timing
}
