#include <iostream>
#include <iomanip>
#include <fstream>
#include <chrono>
#include <map>
#include <Eigen/Dense>
#include "Utility/Algebra.h"
#include "Utility/Geodesy.h"
#include "Utility/String.h"
#include "GUCore/ImuData.h"
#include "GUCore/InsMechanizer.h"
#include "GUCore/GnssData.h"
#include "GUCore/GnssInsFilter.h"
#include "Configure/GuConfig.h"
#ifdef WIN32
#include <corecrt_math_defines.h>
#endif

using namespace Utility;
using namespace std::chrono;

constexpr char conf_file[]{ "Configure//conf_open.ini" };
constexpr double max_gap{ 0.011 };

struct TrueResult {
    double second{ 0.0 };
    double blh[3]{ 0.0 };
    double xyz[3]{ 0.0 };
    double vel[3]{ 0.0 };
    EulerAngle att{};
};

void ReadTrueResult(std::vector<TrueResult> &res, const std::string &file_name) {
    std::ifstream ifs{ file_name, std::ios::in };
    if (!ifs.is_open()) {
        std::cout << "Result File Not Found!" << std::endl;
        return;
    }
    std::string line;
    std::vector<std::string> substr;
    TrueResult r;
    // Read Header
    std::map<std::string, double> kv;
    std::vector<std::string> header;
    getline(ifs, line);
    RemoveSpace(line);
    SplitString(line, substr, " ");
    for (const auto &s: substr) header.push_back(s);
    substr.clear();
    // Read Unit
    getline(ifs, line);
    RemoveSpace(line);
    SplitString(line, substr, " ");
    for (int i = 0; i < header.size(); ++i) header.at(i) += substr.at(i);
    substr.clear();
    // Store In Dictionary
    for (const auto & s : header) kv.insert({ s, 0.0 });
    // Parse File Body
    while (getline(ifs, line)) {
        RemoveSpace(line);
        SplitString(line, substr, " ");
        for (int i = 0; i < substr.size(); ++i) kv[header.at(i)] = stod(substr.at(i));
        r.second = kv["GPSTime(s)"];
        r.blh[0] = D2R(kv["Latitude(deg)"]);
        r.blh[1] = D2R(kv["Longitude(deg)"]);
        r.blh[2] = kv["Ellip-Hgt(m)"];
        r.xyz[0] = kv["ECEF-X(m)"];
        r.xyz[1] = kv["ECEF-Y(m)"];
        r.xyz[2] = kv["ECEF-Z(m)"];
        r.vel[0] = kv["Local-VN(m/s)"];
        r.vel[1] = kv["Local-VE(m/s)"];
        r.vel[2] = -kv["Local-VU(m/s)"];
        r.att.yaw = D2R(kv["Heading(deg)"]);
        r.att.roll = D2R(kv["Roll(deg)"]);
        r.att.pitch = D2R(kv["Pitch(deg)"]);
        res.push_back(r);
        substr.clear();
    }
    ifs.close();
}

void ReadIMUFile(const std::string &filename, std::vector<ImuData>& imu_data) {
    if (filename.find(".ASC") != std::string::npos) { // Try Read ASC
        std::ifstream ifs{filename, std::ios::in};
        if (!ifs.is_open()) return;
        std::string line;
        ImuData curr_imu, prev_imu;
        while (getline(ifs, line)) {
            prev_imu = curr_imu;
            if (!curr_imu.ParseAsc(line)) continue;
            if (curr_imu.IsDuplicated(prev_imu)) continue;
            double dt = curr_imu.GetSecond() - prev_imu.GetSecond();
            if (dt > max_gap && dt < 1.0)
                imu_data.push_back(ImuData::Interpolate(prev_imu, curr_imu, prev_imu.GetSecond() + 0.01));
            imu_data.emplace_back(curr_imu);
        }
        ifs.close();
    } else if (filename.find(".imr") != std::string::npos) { // Try Read imr
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
            if (curr_imu.IsDuplicated(prev_imu)) continue;
            double dt = curr_imu.GetSecond() - prev_imu.GetSecond();
            if (dt > max_gap && dt < 1.0)
                imu_data.push_back(ImuData::Interpolate(prev_imu, curr_imu, prev_imu.GetSecond() + 0.01));
            imu_data.emplace_back(curr_imu);
        }
        ifs.close();
    } else {
        std::cout << "IMU File Not Found!" << std::endl;
        return;
    }
}

EulerAngle initialAlignment(const std::vector<ImuData>& imu_data, const double *initial_pos, double tos, double toe) {
    int epoch{ 0 };
    EulerAngle euler;
    ImuData avg_imu;
    for (auto& imu : imu_data) {
        if (imu.GetSecond() < tos) continue;
        if (imu.GetSecond() > toe) break;
        avg_imu += imu;
        ++epoch;
    }
    avg_imu /= epoch;
    avg_imu.StaticAlignment(initial_pos, euler);
    return euler;
}

void ReadPosFile(const std::string &filename, std::vector<GnssData>& pos_results) {
    std::ifstream ifs{ filename, std::ios::in };
    if (!ifs.is_open()) {
        std::cout << "GNSS File Not Found!" << std::endl;
        return;
    }
    std::string line;
    GnssData pos_result;
    for (int i = 0; i < 2; ++i) getline(ifs, line); // Skip The First Two Line For Convenience
    while (getline(ifs, line)) {
        pos_result.Parse(line);
        pos_results.push_back(pos_result);
    }
    ifs.close();
}

void OutputResult(double time, const GnssInsFilter &filter, const GuConfig &conf, std::ofstream &ofs) {
    auto [pos, vel, att] = filter.GetState();
    auto euler = QuaternionToEulerAngle(att);
    const auto *Gb = filter.GetGb();
    const auto *Ab = filter.GetAb();
    const auto *Gs = filter.GetGs();
    const auto *As = filter.GetAs();
    double ned[3];
    BlhToNed(pos, conf.init_pos, ned);
    ofs << std::fixed << std::setprecision(3) << time << ","
        << std::setprecision(9) << ned[0] << "," << ned[1] << "," << ned[2] << ","
        << std::setprecision(6) << vel[0] << "," << vel[1] << "," << vel[2] << ","
        << std::setprecision(3) << R2D(euler.roll) << "," << R2D(euler.pitch) << "," << R2D(euler.yaw) << ","
        << std::scientific << std::setprecision(9)
        << Eigen::Vector3d(Gb[0], Gb[1], Gb[2]).norm() << ","
        << Eigen::Vector3d(Ab[0], Ab[1], Ab[2]).norm() << ","
        << Eigen::Vector3d(Gs[0], Gs[1], Gs[2]).norm() << ","
        << Eigen::Vector3d(As[0], As[1], As[2]).norm() << ","
        << filter.IsZeroUpdate() << "\n";
}

void OutputResult(double time, const GnssInsFilter &filter, const TrueResult &ref, const GuConfig &conf, std::ofstream &ofs) {
    auto [pos, vel, att] = filter.GetState();
    auto euler = QuaternionToEulerAngle(att);
    const auto *Gb = filter.GetGb();
    const auto *Ab = filter.GetAb();
    const auto *Gs = filter.GetGs();
    const auto *As = filter.GetAs();
    double xyz[3];
    BlhToXyz(pos, xyz);
    double d_yaw = R2D(ref.att.yaw - euler.yaw);
    if (d_yaw < -350) d_yaw = M_PI * 2;
    else if (d_yaw > 350) d_yaw = -M_PI * 2;
    else d_yaw = 0;
    ofs << std::fixed << std::setprecision(3) << time << ","
        << std::setprecision(9) << ref.xyz[0] - xyz[0] << "," << ref.xyz[1] - xyz[1] << "," << ref.xyz[2] - xyz[2] << ","
        << std::setprecision(6) << ref.vel[0] - vel[0] << "," << ref.vel[1] - vel[1] << "," << ref.vel[2] - vel[2] << ","
        << std::setprecision(3) << R2D(ref.att.roll - euler.roll) << "," << R2D(ref.att.pitch - euler.pitch) << "," << R2D(d_yaw + ref.att.yaw - euler.yaw) << ","
        << std::scientific << std::setprecision(9)
        << Eigen::Vector3d(Gb[0], Gb[1], Gb[2]).norm() << ","
        << Eigen::Vector3d(Ab[0], Ab[1], Ab[2]).norm() << ","
        << Eigen::Vector3d(Gs[0], Gs[1], Gs[2]).norm() << ","
        << Eigen::Vector3d(As[0], As[1], As[2]).norm() << ","
        << filter.IsZeroUpdate() << "\n";
}

void LooselyCouple(const std::vector<ImuData> &imu_data, const std::vector<GnssData> &gnss_data, const std::vector<TrueResult> &true_results, const GuConfig &conf) {
    GnssInsFilter filter;
    std::ofstream ofs{ conf.output, std::ios::out };
    // Initial Alignment, get initial position, velocity and attitude
    double tos, toe;
    EulerAngle init_euler = conf.init_euler;
    if (conf.is_align) {
        tos = imu_data.at(0).GetSecond();
        toe = tos + conf.align_time;
        init_euler = initialAlignment(imu_data, conf.init_pos, tos, toe);
    }
    else toe = imu_data.at(0).GetSecond();
    std::cout << "Initial Alignment Result: " <<
    "Yaw = " << R2D(init_euler.yaw) << " " <<
    "Roll = " << R2D(init_euler.roll) << " " <<
    "Pitch = " << R2D(init_euler.pitch) << std::endl;
    auto init_att = EulerAngleToQuaternion(init_euler);
    const auto pos = std::make_pair(conf.init_pos, conf.init_std_pos);
    const auto vel = std::make_pair(conf.init_vel, conf.init_std_vel);
    const auto att = std::make_pair(init_att, conf.init_std_att);
    const auto gb = std::make_pair(conf.init_gb, conf.init_std_gb);
    const auto ab = std::make_pair(conf.init_ab, conf.init_std_ab);
    const auto gs = std::make_pair(conf.init_gs, conf.init_std_gs);
    const auto as = std::make_pair(conf.init_as, conf.init_std_as);
    // Loosely Couple by EKF
    filter.SetLeverArm(conf.lever);
    filter.SetMarkovTime(conf.markov_time, conf.markov_time, conf.markov_time, conf.markov_time);
    filter.SetRandomWalk(conf.ARW, conf.VRW);
    // Time Synchronize And Solve
    double i_time, g_time, r_time;
    int i = 0, k = 0, n = 0;
    while(fabs(imu_data.at(i).GetSecond() - toe) > 1.0E-3) ++i;
    while(gnss_data.at(k).GetSecond() < toe) ++k;
    filter.Initialize(pos, vel, att, gb, ab, gs, as, imu_data.at(i));
    for (++i; i < imu_data.size(); ++i) {
        if (k == gnss_data.size() || n == true_results.size()) break;
        i_time = imu_data.at(i).GetSecond();
        g_time = gnss_data.at(k).GetSecond();
        r_time = true_results.at(n).second;
        if (g_time < i_time && fabs(g_time - i_time) > 1.0E-3) {
            ImuData inter_imu = ImuData::Interpolate(imu_data.at(i - 1), imu_data.at(i), g_time);
            filter.ProcessData(inter_imu, gnss_data.at(k++));
        } else if (fabs(g_time - i_time) < 1.0E-3) {
            filter.ProcessData(imu_data.at(i), gnss_data.at(k++));
        } else {
            filter.ProcessData(imu_data.at(i));
        }
        // OutputResult(i_time, filter, conf, ofs); // Output Result
        if (ofs.is_open() && fabs(r_time - i_time) < 1.0E-3) { // Output Result Different With Reference Result
            const auto ref = true_results.at(n++);
            OutputResult(i_time, filter, ref, conf, ofs);
        }
    }
    filter.PrintState();
}

int main() {
    // Read Configure File
    GuConfig conf;
    conf.ParseConfig(conf_file);
    // Read IMU File
    std::vector<ImuData> imu_data;
    imu_data.reserve(150000);
    ReadIMUFile(conf.input_imu, imu_data);
    // Read GNSS File
    std::vector<GnssData> pos_results;
    ReadPosFile(conf.input_gnss, pos_results);
    // Read Reference Result
    std::vector<TrueResult> true_results;
    ReadTrueResult(true_results, "..//Data//ref_open.pos");
    // Start To Solve
    auto start = system_clock::now();                                               // timing
    LooselyCouple(imu_data, pos_results, true_results, conf);
    auto end = system_clock::now();                                                 // timing
    auto duration = duration_cast<seconds>(end - start);                            // timing
    std::cout << "duration: " << duration.count() << " seconds" << std::endl;       // timing
}
