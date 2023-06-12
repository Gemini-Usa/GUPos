//
// Created by 崔宇璐 on 2023/5/16.
//

#ifndef GUPOS_INSMECHANIZER_H
#define GUPOS_INSMECHANIZER_H

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Core>
#include "ImuData.h"

using namespace Eigen;

class InsMechanizer {
public:
    void Initialize(const double* pos, const double* vel, const Eigen::Quaterniond& att, const ImuData& imu);
    void INSUpdate(const ImuData& data);
    void PrintState() const;
    std::tuple<const double *, const double *, Eigen::Quaterniond> getState() const;
    double getTimeInterval(double curr_time) const;
    void Correct(const Eigen::Vector3d &dp, const Eigen::Vector3d &dv, const Eigen::Quaterniond &dq);
private:
    void positionUpdate(const ImuData &curr_imu, const double *curr_vel, double *pos);
    void velocityUpdate(const ImuData &curr_imu, double *vel);
    void attitudeUpdate(const ImuData &curr_imu, Eigen::Quaterniond& att);
public:
    const double *getVel() const;
    const double *getPos() const;
    const Quaterniond &getAtt() const;
private:
    ImuData _imu_data;
    std::deque<std::tuple<double*, double*>> _state_queue;
    double _pos[3]{ 0.0, 0.0, 0.0 };
    double _vel[3]{ 0.0, 0.0, 0.0 };
    Eigen::Quaterniond _att;
};


#endif //GUPOS_INSMECHANIZER_H
