//
// Created by 崔宇璐 on 2023/6/8.
//

#ifndef GUPOS_GNSSINSFILTER_H
#define GUPOS_GNSSINSFILTER_H

#include "KalmanFilter.hpp"
#include "InsMechanizer.h"
#include "GnssData.h"

class GnssInsFilter : KalmanFilter<21, 6> {
public:
    void Initialize(const std::pair<const double *, const double *> &pos, const std::pair<const double *, const double *> &vel,
                    const std::pair<const Quaterniond &, const double *> &att,
                    const std::pair<const double *, const double *> &gyro_bias, const std::pair<const double *, const double *> &accl_bias,
                    const std::pair<const double *, const double *> &gyro_scale, const std::pair<const double *, const double *> &accl_scale,
                    const ImuData &imu);
    void setMarkovTime(double Tgb, double Tab, double Tgs, double Tas);
    void setLeverArm(const double *lever);
    void setRandomWalk(double ARW, double VRW);
    void ProcessData(const ImuData &imu, const GnssData &gnss = GnssData());
    void PrintState() const;
private:
    static void setF_rr(Fdim &F, double vn, double ve, double vd, double phi, double h);
    static void setF_vr(Fdim &F, double vn, double ve, double vd, double phi, double h);
    static void setF_vv(Fdim &F, double vn, double ve, double vd, double phi, double h);
    static void setF_phir(Fdim &F, double vn, double ve, double vd, double phi, double h);
    static void setF_phiv(Fdim &F, double phi, double h);
    static void setF_Anginn(Fdim &F, double vn, double ve, double phi, double h);
    Fdim buildF(const ImuData &imu, double dt) const;
    Hdim buildH(const ImuData &imu, double dt) const;
    Qdim buildQ(const Fdim &F, double dt) const;
    zdim buildz(const ImuData &imu, const GnssData &gnss, double dt) const;
    Rdim buildR(const GnssData &gnss) const;
private:
    InsMechanizer _ins;
    double _gb[3]{ 0.0 }, _ab[3]{ 0.0 }, _gs[3]{ 0.0 }, _as[3]{ 0.0 };
    double _std_gb[3]{ 0.0 }, _std_ab[3]{ 0.0 }, _std_gs[3]{ 0.0 }, _std_as[3]{ 0.0 };
    double _Tgb{ 0.0 }, _Tab{ 0.0 }, _Tgs{ 0.0 }, _Tas{ 0.0 };
    double _ARW{ 0.0 }, _VRW{ 0.0 };
    double _lever[3]{ 0.0, 0.0, 0.0 };
};


#endif //GUPOS_GNSSINSFILTER_H
