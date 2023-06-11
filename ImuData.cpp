//
// Created by 崔宇璐 on 2023/5/16.
//

#include <string>
#include <vector>
#include <Eigen/Dense>
#include "ImuData.h"
#include "Algebra.h"
#include "Geodesy.h"

using namespace Eigen;
using namespace Utility;

void SplitString(const std::string& str, std::vector<std::string>& substr, const std::string& split) {
    std::string str0;
    for (auto c : str) {
        if (split.find(c) == std::string::npos) str0.push_back(c);
        else { substr.push_back(str0); str0.clear(); }
    }
}

bool ImuData::ParseASC(const std::string &str) {
    std::string split{ "," };
    split += ";";
    split += "*";
    std::vector<std::string> substr;
    SplitString(str, substr, split);
    if (substr[0] != "%RAWIMUSA") return false;
    this->_second = std::stod(substr[2]);
    this->_accl[0] = -std::stod(substr[7]) * _accl_scale * _frequency;
    this->_accl[1] = std::stod(substr[8]) * _accl_scale * _frequency;
    this->_accl[2] = -std::stod(substr[6]) * _accl_scale * _frequency;
    this->_gyro[0] = -std::stod(substr[10]) * _gyro_scale * _frequency;
    this->_gyro[1] = std::stod(substr[11]) * _gyro_scale * _frequency;
    this->_gyro[2] = -std::stod(substr[9]) * _gyro_scale * _frequency;
    return true;
}

bool ImuData::isDuplicated(const ImuData &other) const {
    if (fabs(this->_second - other._second) < 1.0e-3) return true;
    else return false;
}

double ImuData::getSecond() const {
    return _second;
}

void ImuData::SmoothBy(const std::deque<ImuData> &dataset) {
    auto size = dataset.size();
    for (int i = 0; i < 3; ++i) {
        this->_accl[i] = 0.0;
        this->_gyro[i] = 0.0;
    }
    for (auto& imu : dataset) {
        for (int i = 0; i < 3; ++i) {
            this->_accl[i] += imu._accl[i];
            this->_gyro[i] += imu._gyro[i];
        }
    }
    for (int i = 0; i < 3; ++i) {
        this->_accl[i] /= static_cast<double>(size);
        this->_gyro[i] /= static_cast<double>(size);
    }
}

void ImuData::StaticAlignment(const double *blh, EulerAngle &att) const {
    Vector3d g_n{ 0.0, 0.0, getLocalGravity(blh[0], blh[2]) };
    Vector3d omg_ien{ omg_e * cos(blh[0]), 0.0, -omg_e * sin(blh[0]) };
    Vector3d g_b{ -this->_accl[0], -this->_accl[1], -this->_accl[2] };
    Vector3d omg_ieb{ this->_gyro[0], this->_gyro[1], this->_gyro[2] };

    Vector3d v_g = g_n.normalized();
    Vector3d v_omg = g_n.cross(omg_ien).normalized();
    Vector3d v_gomg = g_n.cross(omg_ien).cross(g_n).normalized();
    Vector3d w_g = g_b.normalized();
    Vector3d w_omg = g_b.cross(omg_ieb).normalized();
    Vector3d w_gomg = g_b.cross(omg_ieb).cross(g_b).normalized();

    Matrix3d V, W, C_bn;
    V.block<1, 3>(0, 0) = v_g;
    V.block<1, 3>(1, 0) = v_omg;
    V.block<1, 3>(2, 0) = v_gomg;
    W.block<1, 3>(0, 0) = w_g;
    W.block<1, 3>(1, 0) = w_omg;
    W.block<1, 3>(2, 0) = w_gomg;
    C_bn = V.inverse() * W;

    att = MatrixToEulerAngle(C_bn);
}

ImuData& ImuData::operator=(const ImuData& other) {
    if (&other == this) return *this;
    this->_second = other._second;
    for (int i = 0; i < 3; ++i) {
        this->_accl[i] = other._accl[i];
        this->_gyro[i] = other._gyro[i];
    }
    return *this;
}

ImuData &ImuData::operator+=(const ImuData &other) {
    if (&other == this) return *this;
    for (int i = 0; i < 3; ++i) this->_accl[i] += other._accl[i];
    for (int i = 0; i < 3; ++i) this->_gyro[i] += other._gyro[i];
    return *this;
}

ImuData &ImuData::operator/=(int num) {
    for (int i = 0; i < 3; ++i) {
        this->_accl[i] /= num;
        this->_gyro[i] /= num;
    }
    return *this;
}

Eigen::Vector3d ImuData::getAccl() const {
    return Eigen::Vector3d{ _accl[0], _accl[1], _accl[2] } / _frequency;
}

Eigen::Vector3d ImuData::getGyro() const {
    return Eigen::Vector3d{ _gyro[0], _gyro[1], _gyro[2] } / _frequency;
}

ImuData ImuData::Interpolate(const ImuData &prev_imu, const ImuData &curr_imu, double time) {
    double dt = time - prev_imu._second;
    ImuData midd_imu;
    midd_imu = prev_imu + (curr_imu - prev_imu) * dt;
    midd_imu._second = time;
    return midd_imu;
}

ImuData operator+(const ImuData &imu1, const ImuData &imu2) {
    ImuData res;
    for (int i = 0; i < 3; ++i) {
        res._accl[i] = imu1._accl[i] + imu2._accl[i];
        res._gyro[i] = imu1._gyro[i] + imu2._gyro[i];
    }
    return res;
}

ImuData operator-(const ImuData &imu1, const ImuData &imu2) {
    ImuData res;
    for (int i = 0; i < 3; ++i) {
        res._accl[i] = imu1._accl[i] - imu2._accl[i];
        res._gyro[i] = imu1._gyro[i] - imu2._gyro[i];
    }
    return res;
}

ImuData operator*(const ImuData &imu, double num) {
    ImuData res;
    for (int i = 0; i < 3; ++i) {
        res._accl[i] = imu._accl[i] * num;
        res._gyro[i] = imu._gyro[i] * num;
    }
    return res;
}

void ImuData::Compensate(const double *gb, const double *ab, const double *gs, const double *as) {
    for (int i = 0; i < 3; ++i) {
        _gyro[i] = (_gyro[i] - gb[i]) / (1.0 + gs[i]);
        _accl[i] = (_accl[i] - ab[i]) / (1.0 + as[i]);
    }
}
