//
// Created by 崔宇璐 on 2023/5/16.
//

#include <cmath>
#ifdef _WIN32
#include <corecrt_math_defines.h>
#endif
#include "Geodesy.h"

double Utility::getLocalGravity(double phi, double h) {
    double g0 = (9.7803267715 * (1.0 + 0.0052790414 * sin(phi) * sin(phi) + 0.0000232718 * pow(sin(phi), 4)));
    return g0 - (3.087691089E-6 - 4.397731E-9 * pow(sin(phi), 2)) * h + 0.721E-12 * h * h;
}

double Utility::D2R(double deg) {
    return deg * M_PI / 180.0;
}

double Utility::getRM(double phi) {
    double sqsinphi = sin(phi) * sin(phi);
    double e2 = Utility::flat * (2.0 - Utility::flat);
    return Utility::major * (1 - e2) / sqrt(pow((1 - e2 * sqsinphi), 3));
}

double Utility::getRN(double phi) {
    double sqsinphi = sin(phi) * sin(phi);
    double e2 = Utility::flat * (2.0 - Utility::flat);
    return Utility::major / sqrt(1 - e2 * sqsinphi);
}

double Utility::R2D(double rad) {
    return rad * 180.0 / M_PI;
}

void Utility::BlhToXyz(const double *blh, double *xyz) {
    double R_N = getRN(blh[0]);
    double e2 = 2 * flat - flat * flat;
    xyz[0] = (R_N + blh[2]) * cos(blh[0]) * cos(blh[1]);
    xyz[1] = (R_N + blh[2]) * cos(blh[0]) * sin(blh[1]);
    xyz[2] = (R_N * (1 - e2) + blh[2]) * sin(blh[0]);
}

void Utility::BlhToNed(const double *blh, const double *blh0, double *ned) {
    double sin_l = sin(blh0[1]);
    double cos_l = cos(blh0[1]);
    double sin_b = sin(blh0[0]);
    double cos_b = cos(blh0[0]);
    Eigen::Matrix3d R;
    R << -sin_b * cos_l, -sin_b * sin_l, cos_b,
    -sin_l, cos_l, 0,
    -cos_b * cos_l,  -cos_b * sin_l, -sin_b;
    double xyz0[3], xyz[3];
    BlhToXyz(blh0, xyz0);
    BlhToXyz(blh, xyz);
    Eigen::Vector3d d_xyz{ xyz[0] - xyz0[0], xyz[1] - xyz0[1], xyz[2] - xyz0[2] };
    Eigen::Vector3d res = R * d_xyz;
    for (int i = 0; i < 3; ++i) ned[i] = res(i);
}

Eigen::Vector3d Utility::getAng_ienVec(double phi) {
    Eigen::Vector3d Ang_ien;
    Ang_ien << omg_e * cos(phi), 0, -omg_e * sin(phi);
    return Ang_ien;
}

Eigen::Vector3d Utility::getAng_ennVec(double phi, double h, double vn, double ve) {
    Eigen::Vector3d Ang_enn;
    double R_N = getRN(phi);
    double R_M = getRM(phi);
    Ang_enn << ve / (R_N + h), -vn / (R_M + h), -ve * tan(phi) / (R_N + h);
    return Ang_enn;
}

Eigen::Vector3d Utility::getAng_innVec(double phi, double h, double vn, double ve) {
    Eigen::Vector3d Ang_ien = getAng_ienVec(phi);
    Eigen::Vector3d Ang_enn = getAng_ennVec(phi, h, vn, ve);
    return Ang_ien + Ang_enn;
}

Eigen::Matrix3d Utility::getInv_DR(double phi, double h) {
    double R_M = getRM(phi);
    double R_N = getRN(phi);
    return Eigen::Vector3d(1 / (R_M + h), 1 / ((R_N + h) * cos(phi)), -1).asDiagonal();
}
