//
// Created by GeminiUsa on 2023/6/12.
//

#ifndef GUPOS_GUCONFIG_H
#define GUPOS_GUCONFIG_H

#include <Eigen/Dense>
#include "../Utility/Algebra.h"

class GuConfig {
public:
    static double ParseValue(const std::string &qty) ;
    void ParseConfig(const std::string &conf_file);
public:
    double init_pos[3];
    double init_vel[3];
    Utility::EulerAngle init_euler;
    double init_gb[3];
    double init_ab[3];
    double init_gs[3];
    double init_as[3];
    double init_std_pos[3];
    double init_std_vel[3];
    double init_std_att[3];
    double init_std_gb[3];
    double init_std_ab[3];
    double init_std_gs[3];
    double init_std_as[3];
    bool is_ZUPT;
    double lever[3];
    double ARW;
    double VRW;
    double markov_time;
    bool is_align;
    double align_time;
    std::string input_imu;
    std::string input_gnss;
    std::string output;
};


#endif //GUPOS_GUCONFIG_H
