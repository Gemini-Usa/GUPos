//
// Created by 崔宇璐 on 2023/5/16.
//

#ifndef GUPOS_IMUDATA_H
#define GUPOS_IMUDATA_H

#include <cmath>
#include <deque>
#include <Eigen/Dense>
#include "Algebra.h"

class ImuData {
public:
    friend ImuData operator+(const ImuData& imu1, const ImuData& imu2);
    friend ImuData operator-(const ImuData &imu1, const ImuData &imu2);
    friend ImuData operator*(const ImuData& imu, double num);
    ImuData& operator=(const ImuData& other);
    ImuData& operator+=(const ImuData& other);
    ImuData& operator/=(int num);
    bool ParseASC(const std::string& str);
    bool isDuplicated(const ImuData& other) const;
    void SmoothBy(const std::deque<ImuData> &dataset);
    void StaticAlignment(const double *blh, Utility::EulerAngle &att) const;
    void Compensate(const double *gb, const double *ab, const double *gs, const double *as);
    static ImuData Interpolate(const ImuData& prev_imu, const ImuData& curr_imu, double time);
public: // Getter
    double getSecond() const;
    Eigen::Vector3d getAccl() const;
    Eigen::Vector3d getGyro() const;
    static const int getFrequency();
private:
    static constexpr int _frequency{ 100 };
    static constexpr double _gyro_scale{ 0.1 / (3600.0 * 256.0) };
    static constexpr double _accl_scale{ 0.05 / 32768.0 };
    double _second{ 0.0 };
    double _accl[3]{ 0.0, 0.0, 0.0 };
    double _gyro[3]{ 0.0, 0.0, 0.0 };
};


#endif //GUPOS_IMUDATA_H
