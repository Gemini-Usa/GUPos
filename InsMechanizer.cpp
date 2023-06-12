//
// Created by 崔宇璐 on 2023/5/16.
//

#include "InsMechanizer.h"
#include "Algebra.h"
#include "Geodesy.h"
using namespace Utility;
using v3 = Eigen::Vector3d;
using qt = Eigen::Quaterniond;

void InsMechanizer::Initialize(const double *pos, const double *vel, const qt& att, const ImuData& imu) {
    for (int i = 0; i < 3; ++i) {
        _pos[i] = pos[i];
        _vel[i] = vel[i];
    }
    _att = att;
    _imu_data = imu;
}

void InsMechanizer::positionUpdate(const ImuData &curr_imu, const double *curr_vel, double *pos) {
    double t = curr_imu.getSecond() - _imu_data.getSecond();
    double R_M = getRM(_pos[0]);
    double R_N = getRN(_pos[0]);
    double avg_veln = (curr_vel[0] + _vel[0]) / 2.0;
    double avg_vele = (curr_vel[1] + _vel[1]) / 2.0;
    double avg_veld = (curr_vel[2] + _vel[2]) / 2.0;
    pos[2] = _pos[2] - avg_veld * t;
    double avg_h = (pos[2] - _pos[2]) / 2.0;
    pos[0] = _pos[0] + avg_veln * t / (R_M + avg_h);
    double avg_b = (pos[0] + _pos[0]) / 2.0;
    pos[1] = _pos[1] + avg_vele * t / ((R_N + avg_h) * cos(avg_b));
}

void InsMechanizer::velocityUpdate(const ImuData &curr_imu, double *vel) {
    double t = curr_imu.getSecond() - _imu_data.getSecond();
    double p[3]{ _pos[0], _pos[1], _pos[2] }, v[3]{ _vel[0], _vel[1], _vel[2] };
    if (_state_queue.size() == 2) { // extrapolate
        const auto& [pprev_p, pprev_v] = _state_queue.front();
        const auto& [prev_p, prev_v] = _state_queue.back();
        for (int i = 0; i < 3; ++i) {
            p[i] = prev_p[i] * (3.0 / 2.0) - pprev_p[i] * 0.5;
            v[i] = prev_v[i] * (3.0 / 2.0) - pprev_v[i] * 0.5;
        }
    }
    auto curr_gyro = curr_imu.getGyro();
    auto prev_gyro = _imu_data.getGyro();
    auto curr_accl = curr_imu.getAccl();
    auto prev_accl = _imu_data.getAccl();
    v3 g_pn{ 0.0, 0.0, getLocalGravity(p[0], p[2]) };
    v3 omg_ien = getAng_ienVec(p[0]);
    v3 omg_enn = getAng_ennVec(p[0], p[2], v[0], v[1]);
    v3 a_gc = g_pn - (2 * omg_ien + omg_enn).cross(v3{ v[0], v[1], v[2] });
    v3 dv_gn = a_gc * t;
    v3 zeta = (omg_ien + omg_enn) * t;
    v3 rotate = curr_gyro.cross(curr_accl) / 2.0;
    v3 scull = (prev_gyro.cross(curr_accl) + prev_accl.cross(curr_gyro)) / 12.0;
    v3 dv_fb = curr_accl + rotate + scull;
    v3 dv_fn = (Eigen::Matrix3d::Identity() - SkewSymmetric(zeta) * 0.5) * _att.toRotationMatrix() * dv_fb;
    for (int i = 0; i < 3; ++i) vel[i] = _vel[i] + dv_fn(i) + dv_gn(i);
}

void InsMechanizer::attitudeUpdate(const ImuData &curr_imu, Eigen::Quaterniond& att) {
    double t = curr_imu.getSecond() - _imu_data.getSecond();
    auto curr_gyro = curr_imu.getGyro();
    auto prev_gyro = _imu_data.getGyro();
    v3 omg_ien = getAng_ienVec(_pos[0]);
    v3 omg_enn = getAng_ennVec(_pos[0], _pos[2], _vel[0], _vel[1]);
    v3 phi_k = curr_gyro + prev_gyro.cross(curr_gyro) / 12.0;
    qt q_b = RotateVectorToQuaternion(phi_k);
    v3 zeta = (omg_ien + omg_enn) * t;
    qt q_n = RotateVectorToQuaternion(zeta).conjugate();
    att = (q_n * _att * q_b).normalized();
}

void InsMechanizer::INSUpdate(const ImuData &data) {
    double curr_pos[3]{ 0.0, 0.0, 0.0 }, curr_vel[3]{ 0.0, 0.0, 0.0 };
    qt curr_att{};
    attitudeUpdate(data, curr_att);
    velocityUpdate(data, curr_vel);
    positionUpdate(data, curr_vel, curr_pos);
    for (int i = 0; i < 3; ++i) {
        _vel[i] = curr_vel[i];
        _pos[i] = curr_pos[i];
    }
    _att = curr_att;
    _imu_data = data;
    if (_state_queue.size() == 2) _state_queue.pop_front();
    _state_queue.emplace_back(_pos, _vel);
}

void InsMechanizer::PrintState() const {
    printf("time %.3f\n", _imu_data.getSecond());
    printf("\t B: %.9f, L: %.9f, H: %.6f\n", R2D(_pos[0]), R2D(_pos[1]), _pos[2]);
    printf("\t VN: %.6f, VE: %.6f, VD: %.6f\n", _vel[0], _vel[1], _vel[2]);
    EulerAngle euler_angle = QuaternionToEulerAngle(_att);
    printf("\t roll: %.6f, pitch: %.6f, yaw: %.6f\n", R2D(euler_angle.roll), R2D(euler_angle.pitch), R2D(euler_angle.yaw));
}

std::string InsMechanizer::getStateInfo() const {
    char str[200];
    EulerAngle euler_angle = QuaternionToEulerAngle(_att);
    std::sprintf(str, "%10.3f %13.9f %13.9f %13.6f %13.9f %13.9f %13.9f %10.3f %10.3f %10.3f",
         _imu_data.getSecond(), R2D(_pos[0]), R2D(_pos[1]), _pos[2],
         _vel[0], _vel[1], _vel[2],
         R2D(euler_angle.roll), R2D(euler_angle.pitch), R2D(euler_angle.yaw));
    return {str};
}

double InsMechanizer::getTimeInterval(double curr_time) const {
    return curr_time - _imu_data.getSecond();
}

const Quaterniond &InsMechanizer::getAtt() const {
    return _att;
}

const double *InsMechanizer::getPos() const {
    return _pos;
}

const double *InsMechanizer::getVel() const {
    return _vel;
}

void InsMechanizer::Correct(const Vector3d &dp, const Vector3d &dv, const Quaterniond &dq) {
    for (int i = 0; i < 3; ++i) {
        _pos[i] -= dp(i);
        _vel[i] -= dv(i);
    }
    _att = dq * _att;
}
