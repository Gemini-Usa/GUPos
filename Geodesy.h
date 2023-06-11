//
// Created by 崔宇璐 on 2023/5/16.
//

#ifndef GUPOS_GEODESY_H
#define GUPOS_GEODESY_H

#include <Eigen/Dense>

namespace Utility {
    constexpr double omg_e{ 7.2921150 * 1.0e-5 };
    constexpr double major{ 6378137.0 };
    constexpr double flat{ 1.0/298.257223563 };

    double D2R(double deg);
    double R2D(double rad);
    double getRM(double phi);
    double getRN(double phi);
    double getLocalGravity(double phi, double h);
    Eigen::Vector3d getAng_ienVec(double phi);
    Eigen::Vector3d getAng_ennVec(double phi, double h, double vn, double ve);
    Eigen::Vector3d getAng_innVec(double phi, double h, double vn, double ve);
    Eigen::Matrix3d getinv_DR(double phi, double h);
};


#endif //GUPOS_GEODESY_H
