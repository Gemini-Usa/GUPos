//
// Created by 崔宇璐 on 2023/6/12.
//

#include "GuConfig.h"
#include "Algebra.h"
#include "Geodesy.h"
#include "minIni.h"

using namespace Utility;

void GuConfig::ParseConfig(const std::string &conf_file) {
    minIni ini(conf_file);
    // Position BLH
    init_pos[0] = D2R(ini.getf("Public", "init_pos_B"));
    init_pos[1] = D2R(ini.getf("Public", "init_pos_L"));
    init_pos[2] = ini.getf("Public", "init_pos_H");
    // Velocity v_ned
    init_vel[0] = ini.getf("Public", "init_vel_N");
    init_vel[1] = ini.getf("Public", "init_vel_E");
    init_vel[2] = ini.getf("Public", "init_vel_D");
    // Attitude euler angle
    init_euler.yaw = D2R(ini.getf("Public", "init_att_Y"));
    init_euler.roll = D2R(ini.getf("Public", "init_att_R"));
    init_euler.pitch = D2R(ini.getf("Public", "init_att_P"));
    // Gyro, Accelerator bias & scale factor
    init_gb[0] = init_gb[1] = init_gb[2] = ini.getf("Public", "init_gb");
    init_ab[0] = init_ab[1] = init_ab[2] = ini.getf("Public", "init_ab");
    init_gs[0] = init_gs[1] = init_gs[2] = ini.getf("Public", "init_gs");
    init_as[0] = init_as[1] = init_as[2] = ini.getf("Public", "init_as");
    // Standard variance
    init_std_pos[0] = ini.getf("Public", "init_std_N");
    init_std_pos[1] = ini.getf("Public", "init_std_E");
    init_std_pos[2] = ini.getf("Public", "init_std_D");
    init_std_vel[0] = ini.getf("Public", "init_std_VN");
    init_std_vel[1] = ini.getf("Public", "init_std_VE");
    init_std_vel[2] = ini.getf("Public", "init_std_VD");
    init_std_att[0] = D2R(ini.getf("Public", "init_std_R"));
    init_std_att[1] = D2R(ini.getf("Public", "init_std_Y"));
    init_std_att[2] = D2R(ini.getf("Public", "init_std_P"));
    init_std_gb[0] = init_std_gb[1] = init_std_gb[2] = D2R(ini.getf("Public", "init_std_gb"));
    init_std_ab[0] = init_std_ab[1] = init_std_ab[2] = ini.getf("Public", "init_std_ab");
    init_std_gs[0] = init_std_gs[1] = init_std_gs[2] = ini.getf("Public", "init_std_gs");
    init_std_as[0] = init_std_as[1] = init_std_as[2] = ini.getf("Public", "init_std_as");
    // Other important parameter
    lever[0] = ini.getf("Public", "lever_x");
    lever[1] = ini.getf("Public", "lever_y");
    lever[2] = ini.getf("Public", "lever_z");
    ARW = D2R(ini.getf("Public", "ARW"));
    VRW = ini.getf("Public", "VRW");
    markov_time = ini.getf("Public", "markov_time");
    // Initial alignment
    is_align = ini.getbool("Initial Alignment", "is_align", false);
    align_time = ini.getf("Initial Alignment", "align_time", 0.0);
    // File
    input_imu = ini.gets("File", "input_imu");
    input_gnss = ini.gets("File", "input_gnss");
    output = ini.gets("File", "output", "");
}
