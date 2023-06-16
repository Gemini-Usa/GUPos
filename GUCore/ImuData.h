//
// Created by 崔宇璐 on 2023/5/16.
//

#ifndef GUPOS_IMUDATA_H
#define GUPOS_IMUDATA_H

#include <cmath>
#include <deque>
#include <Eigen/Dense>
#include "../Utility/Algebra.h"

class ImuData {
public:
    friend ImuData operator+(const ImuData& imu1, const ImuData& imu2);
    friend ImuData operator-(const ImuData &imu1, const ImuData &imu2);
    friend ImuData operator*(const ImuData& imu, double num);
    ImuData& operator=(const ImuData& other);
    ImuData& operator+=(const ImuData& other);
    ImuData& operator/=(int num);
    bool ParseAsc(const std::string &str);
    static void ParseImrHeader(const char *header);
    void ParseImr(const char *str);
    bool IsDuplicated(const ImuData& other) const;
    void SmoothBy(const std::deque<ImuData> &dataset);
    void StaticAlignment(const double *blh, Utility::EulerAngle &att) const;
    void Compensate(const double *gb, const double *ab, const double *gs, const double *as);
    static ImuData Interpolate(const ImuData& prev_imu, const ImuData& curr_imu, double time);
public: // Getter
    double GetSecond() const;
    Eigen::Vector3d GetAccl() const;
    Eigen::Vector3d GetGyro() const;
    static int GetRate();
private:
    static int _rate;
    static double _gyro_scale;
    static double _accl_scale;
    double _second{ 0.0 };
    double _accl[3]{ 0.0, 0.0, 0.0 };
    double _gyro[3]{ 0.0, 0.0, 0.0 };
};


#endif //GUPOS_IMUDATA_H
