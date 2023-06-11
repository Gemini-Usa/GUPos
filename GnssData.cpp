//
// Created by 崔宇璐 on 2023/5/29.
//

#include "GnssData.h"
#include "Geodesy.h"

void GnssData::Parse(const std::string &str) {
    _second = stod(str.substr(9, 13));
    _blh[0] = Utility::D2R(stod(str.substr(22, 19)));
    _blh[1] = Utility::D2R(stod(str.substr(41, 19)));
    _blh[2] = stod(str.substr(60, 15));
    _blh_std[0] = stod(str.substr(75, 15));
    _blh_std[1] = stod(str.substr(90, 15));
    _blh_std[2] = stod(str.substr(105, 15));
    _vel[0] = stod(str.substr(120, 12));
    _vel[1] = stod(str.substr(132, 12));
    _vel[2] = stod(str.substr(144, 12));
    _vel_std[0] = stod(str.substr(156, 12));
    _vel_std[1] = stod(str.substr(168, 12));
    _vel_std[2] = stod(str.substr(180, 12));
}

double GnssData::getSecond() const {
    return _second;
}

const double *GnssData::getBlh() const {
    return _blh;
}

const double *GnssData::getBlhStd() const {
    return _blh_std;
}

const double *GnssData::getVel() const {
    return _vel;
}

const double *GnssData::getVelStd() const {
    return _vel_std;
}
