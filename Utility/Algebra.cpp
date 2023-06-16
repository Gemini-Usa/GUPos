//
// Created by GeminiUsa on 2023/5/16.
//
#ifdef _WIN32
#include <corecrt_math_defines.h>
#endif
#include "Algebra.h"

Utility::EulerAngle Utility::MatrixToEulerAngle(const Eigen::Matrix3d &DCM) {
    EulerAngle euler_angle;
    double pitch = atan2(-DCM(2, 0), sqrt(DCM(2, 1) * DCM(2, 1) + DCM(2, 2) * DCM(2, 2)));
    double yaw = atan2(DCM(1, 0), DCM(0, 0));
    double roll = atan2(DCM(2, 1), DCM(2, 2));
    euler_angle.pitch = pitch;
    euler_angle.yaw = yaw;
    euler_angle.roll = roll;
    return euler_angle;
}

Utility::EulerAngle Utility::QuaternionToEulerAngle(const Eigen::Quaterniond& q) {
    Utility::EulerAngle angles;

    // roll (x-axis rotation)
    double sinr_cosp = 2 * (q.w() * q.x() + q.y() * q.z());
    double cosr_cosp = 1 - 2 * (q.x() * q.x() + q.y() * q.y());
    angles.roll = std::atan2(sinr_cosp, cosr_cosp);

    // pitch (y-axis rotation)
    double sinp = std::sqrt(1 + 2 * (q.w() * q.y() - q.x() * q.z()));
    double cosp = std::sqrt(1 - 2 * (q.w() * q.y() - q.x() * q.z()));
    angles.pitch = 2 * std::atan2(sinp, cosp) - M_PI / 2;

    // yaw (z-axis rotation)
    double siny_cosp = 2 * (q.w() * q.z() + q.x() * q.y());
    double cosy_cosp = 1 - 2 * (q.y() * q.y() + q.z() * q.z());
    angles.yaw = std::atan2(siny_cosp, cosy_cosp);
    if (angles.yaw < 0) angles.yaw += 2 * M_PI;
    if (angles.yaw > 2 * M_PI) angles.yaw -= 2 * M_PI;

    return angles;
}

Eigen::Matrix3d Utility::SkewSymmetric(const Eigen::Vector3d& vector) {
    Eigen::Matrix3d matrix;
    matrix << 0.0, -vector(2), vector(1),
            vector(2), 0.0, -vector(0),
            -vector(1), vector(0), 0.0;
    return matrix;
}

Eigen::Quaterniond Utility::RotateVectorToQuaternion(const Eigen::Vector3d& vector) {
    double real = cos((0.5 * vector).norm());
    Eigen::Vector3d image = (0.5 * vector).normalized() * sin((0.5 * vector).norm());
    return Eigen::Quaterniond{ real, image(0), image(1), image(2) };
}

Eigen::Quaterniond Utility::EulerAngleToQuaternion(const Utility::EulerAngle& euler_angle) {
    double yaw = euler_angle.yaw;
    double pitch = euler_angle.pitch;
    double roll = euler_angle.roll;
    double sinpsi = sin(yaw / 2), cospsi = cos(yaw / 2);
    double sinthe = sin(pitch / 2), costhe = cos(pitch / 2);
    double sinphi = sin(roll / 2), cosphi = cos(roll / 2);
    Eigen::Quaterniond q(cosphi * costhe * cospsi + sinphi * sinthe * sinpsi,
                         sinphi * costhe * cospsi - cosphi * sinthe * sinpsi,
                         cosphi * sinthe * cospsi + sinphi * costhe * sinpsi,
                         cosphi * costhe * sinpsi - sinphi * sinthe * cospsi);
    return q;
}

Utility::EulerAngle &Utility::EulerAngle::operator=(const Utility::EulerAngle &other) {
    if (&other == this) return *this;
    this->pitch = other.pitch;
    this->roll = other.roll;
    this->yaw = other.yaw;
    return *this;
}

Utility::EulerAngle &Utility::EulerAngle::operator+=(const Utility::EulerAngle &other) {
    if (&other == this) return *this;
    this->pitch += other.pitch;
    this->roll += other.roll;
    this->yaw += other.yaw;
    return *this;
}

Utility::EulerAngle &Utility::EulerAngle::operator/=(int num) {
    this->pitch /= num;
    this->roll /= num;
    this->yaw /= num;
    return *this;
}
