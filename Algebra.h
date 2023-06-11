//
// Created by 崔宇璐 on 2023/5/16.
//

#ifndef GUPOS_ALGEBRA_H
#define GUPOS_ALGEBRA_H

#include <Eigen/Dense>

namespace Utility {
    struct EulerAngle{
        double roll{0.0};
        double yaw{0.0};
        double pitch{0.0};
        EulerAngle& operator+=(const EulerAngle& other);
        EulerAngle& operator/=(int num);
    };

    EulerAngle MatrixToEulerAngle(const Eigen::Matrix3d &DCM);
    EulerAngle QuaternionToEulerAngle(const Eigen::Quaterniond& q);
    Eigen::Matrix3d skewSymmetric(const Eigen::Vector3d& vector);
    Eigen::Quaterniond RotateVectorToQuaternion(const Eigen::Vector3d& vector);
    Eigen::Quaterniond EulerAngleToQuaternion(const EulerAngle& euler_angle);
};


#endif //GUPOS_ALGEBRA_H
