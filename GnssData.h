//
// Created by 崔宇璐 on 2023/5/29.
//

#ifndef GUPOS_GNSSDATA_H
#define GUPOS_GNSSDATA_H

#include <iostream>

class GnssData {
public:
    void Parse(const std::string &str);
    double getSecond() const;
    const double *getBlh() const;
    const double *getBlhStd() const;
    const double *getVel() const;
    const double *getVelStd() const;
private:
    double _second{ 0.0 };
    double _blh[3]{ 0.0, 0.0, 0.0 };
    double _blh_std[3]{ 0.0, 0.0, 0.0 };
    double _vel[3]{ 0.0, 0.0, 0.0 };
    double _vel_std[3]{ 0.0, 0.0, 0.0 };
};


#endif //GUPOS_GNSSDATA_H
