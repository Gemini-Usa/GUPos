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
    double GetRM(double phi);
    double GetRN(double phi);
    double GetLocalGravity(double phi, double h);
    void BlhToNed(const double *blh, const double *blh0, double *ned);
    void BlhToXyz(const double *blh, double *xyz);
    Eigen::Vector3d GetAng_ienVec(double phi);
    Eigen::Vector3d GetAng_ennVec(double phi, double h, double vn, double ve);
    Eigen::Vector3d GetAng_innVec(double phi, double h, double vn, double ve);
    Eigen::Matrix3d GetInv_DRMat(double phi, double h);
};


#endif //GUPOS_GEODESY_H
