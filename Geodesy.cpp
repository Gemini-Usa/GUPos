//
// Created by 崔宇璐 on 2023/5/16.
//

#include <cmath>
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

Eigen::Matrix3d Utility::getinv_DR(double phi, double h) {
    double R_M = getRM(phi);
    double R_N = getRN(phi);
    return Eigen::Vector3d(1 / (R_M + h), 1 / ((R_N + h) * cos(phi)), -1).asDiagonal();
}
