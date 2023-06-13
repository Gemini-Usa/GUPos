//
// Created by 崔宇璐 on 2023/6/8.
//

#include "GnssInsFilter.h"
#include "Geodesy.h"
using namespace Utility;

void GnssInsFilter::Initialize(const std::pair<const double *, const double *> &pos, const std::pair<const double *, const double *> &vel,
                               const std::pair<const Quaterniond &, const double *> &att,
                               const std::pair<const double *, const double *> &gyro_bias, const std::pair<const double *, const double *> &accl_bias,
                               const std::pair<const double *, const double *> &gyro_scale, const std::pair<const double *, const double *> &accl_scale,
                               const ImuData &imu) {
    const auto& [p, p_std] = pos;
    const auto& [v, v_std] = vel;
    const auto& [a, a_std] = att;
    const auto& [gb, gb_std] = gyro_bias;
    const auto& [ab, ab_std] = accl_bias;
    const auto& [gs, gs_std] = gyro_scale;
    const auto& [as, as_std] = accl_scale;
    _ins.Initialize(p, v, a, imu);
    for (int i = 0; i < 3; ++i) {
        _std_gb[i] = gb_std[i];
        _std_ab[i] = ab_std[i];
        _std_gs[i] = gs_std[i];
        _std_as[i] = as_std[i];
    }
    Vector<double, 21> std0 = Vector<double, 21>::Zero(); // initial_P
    std0.segment(0, 3) = Eigen::Vector3d(p_std[0], p_std[1], p_std[2]);
    std0.segment(3, 3) = Eigen::Vector3d(v_std[0], v_std[1], v_std[2]);
    std0.segment(6, 3) = Eigen::Vector3d(a_std[0], a_std[1], a_std[2]);
    std0.segment(9, 3) = Eigen::Vector3d(gb_std[0], gb_std[1], gb_std[2]);
    std0.segment(12, 3) = Eigen::Vector3d(ab_std[0], ab_std[1], ab_std[2]);
    std0.segment(15, 3) = Eigen::Vector3d(gs_std[0], gs_std[1], gs_std[2]);
    std0.segment(18, 3) = Eigen::Vector3d(as_std[0], as_std[1], as_std[2]);
    _P = std0.asDiagonal();
    _P *= _P;

    _zupt.WindowMoveOn(imu);
}

void GnssInsFilter::setMarkovTime(double Tgb, double Tab, double Tgs, double Tas) {
    _Tgb = Tgb;
    _Tab = Tab;
    _Tgs = Tgs;
    _Tas = Tas;
}

void GnssInsFilter::setLeverArm(const double *lever) {
    for (int i = 0; i < 3; ++i) _lever[i] = lever[i];
}

void GnssInsFilter::setF_rr(GnssInsFilter::Fdim &F, double vn, double ve, double vd, double phi, double h) {
    double R_M = getRM(phi);
    double R_N = getRN(phi);
    Eigen::Matrix3d Frr = Eigen::Matrix3d::Zero();
    Frr(0, 0) = -vd / (R_M + h);
    Frr(0, 2) = vn / (R_M + h);
    Frr(1, 0) = ve * tan(phi) / (R_N + h);
    Frr(1, 1) = -(vd + vn * tan(phi)) / (R_N + h);
    Frr(1, 2) = ve / (R_N + h);
    F.block(0, 0, 3, 3) = Frr;
}

void GnssInsFilter::setF_vr(GnssInsFilter::Fdim &F, double vn, double ve, double vd, double phi, double h) {
    double R_M = getRM(phi);
    double R_N = getRN(phi);
    double gp = getLocalGravity(phi, h);
    double sinphi = sin(phi);
    double cosphi = cos(phi);
    double tanphi = tan(phi);
    Eigen::Matrix3d Fvr = Eigen::Matrix3d::Zero();
    Fvr(0, 0) = -2 * ve * omg_e * cosphi / (R_M + h) - ve * ve / ((R_M + h) * (R_N + h) * cosphi * cosphi);
    Fvr(0, 2) = vn * vd / ((R_M + h) * (R_M + h)) - ve * ve * tanphi / ((R_N + h) * (R_N + h));
    Fvr(1, 0) = 2 * omg_e * (vn * cosphi - vd * sinphi) / (R_M + h) + vn * ve / ((R_M + h) * (R_N + h) * cosphi * cosphi);
    Fvr(1, 2) = (ve * vd + vn * ve * tanphi) / ((R_N + h) * (R_N + h));
    Fvr(2, 0) = 2 * omg_e * ve * sinphi / (R_M + h);
    Fvr(2, 2) = -ve * ve / ((R_N + h) * (R_N + h)) - vn * vn / ((R_M + h) * (R_M + h)) + 2 * gp / (sqrt(R_M * R_N) + h);
    F.block(3, 0, 3, 3) = Fvr;
}

void GnssInsFilter::setF_vv(GnssInsFilter::Fdim &F, double vn, double ve, double vd, double phi, double h) {
    double R_M = getRM(phi);
    double R_N = getRN(phi);
    double sinphi = sin(phi);
    double cosphi = cos(phi);
    double tanphi = tan(phi);
    Eigen::Matrix3d Fvv = Eigen::Matrix3d::Zero();
    Fvv(0, 0) = vd / (R_M + h);
    Fvv(0, 1) = -2 * (omg_e * sinphi + ve * tanphi / (R_N + h));
    Fvv(0, 2) = vn / (R_M + h);
    Fvv(1, 0) = 2 * omg_e * sinphi + ve * tanphi / (R_N + h);
    Fvv(1, 1) = (vd + vn * tanphi) / (R_N + h);
    Fvv(1, 2) = 2 * omg_e * cosphi + ve / (R_N + h);
    Fvv(2, 0) = -2 * vn / (R_M + h);
    Fvv(2, 1) = -2 * (omg_e * cosphi + ve / (R_N + h));
    F.block(3, 3, 3, 3) = Fvv;
}

void GnssInsFilter::setF_phir(GnssInsFilter::Fdim &F, double vn, double ve, double vd, double phi, double h) {
    double R_M = getRM(phi);
    double R_N = getRN(phi);
    double sinphi = sin(phi);
    double cosphi = cos(phi);
    double tanphi = tan(phi);
    Eigen::Matrix3d Fphir = Eigen::Matrix3d::Zero();
    Fphir(0, 0) = -omg_e * sinphi / (R_M + h);
    Fphir(0, 2) = ve / ((R_N + h) * (R_N + h));
    Fphir(1, 2) = -vn / ((R_M + h) * (R_M + h));
    Fphir(2, 0) = -omg_e * cosphi / (R_M + h) - ve / ((R_M + h) * (R_N + h) * cosphi * cosphi);
    Fphir(2, 2) = -ve * tanphi / ((R_N + h) * (R_N + h));
    F.block(6, 0, 3, 3) = Fphir;
}

void GnssInsFilter::setF_phiv(GnssInsFilter::Fdim &F, double phi, double h) {
    double R_M = getRM(phi);
    double R_N = getRN(phi);
    double tanphi = tan(phi);
    Eigen::Matrix3d Fphiv = Eigen::Matrix3d::Zero();
    Fphiv(0, 1) = 1 / (R_N + h);
    Fphiv(1, 0) = -1 / (R_M + h);
    Fphiv(2, 1) = -tanphi / (R_N + h);
    F.block(6, 3, 3, 3) = Fphiv;
}

void GnssInsFilter::setF_Anginn(GnssInsFilter::Fdim &F, double vn, double ve, double phi, double h) {
    Eigen::Matrix3d FAng_inn = Eigen::Matrix3d::Zero();
    double R_M = getRM(phi);
    double R_N = getRN(phi);
    double sinphi = sin(phi);
    double cosphi = cos(phi);
    double tanphi = tan(phi);
    Eigen::Vector3d Ang_inn = Eigen::Vector3d::Zero();
    Ang_inn(0) = omg_e * cosphi + ve / (R_N + h);
    Ang_inn(1) = -vn / (R_M + h);
    Ang_inn(2) = -omg_e * sinphi - ve * tanphi / (R_N + h);
    FAng_inn = SkewSymmetric({Ang_inn(0), Ang_inn(1), Ang_inn(2)});
    F.block(6, 6, 3, 3) = -FAng_inn;
}

GnssInsFilter::Fdim GnssInsFilter::buildF(const ImuData &imu, double dt) const {
    GnssInsFilter::Fdim F{ GnssInsFilter::Fdim::Zero() };
    auto ang = imu.getGyro();
    auto spc = imu.getAccl();
    Eigen::Vector3d fb{ spc[0], spc[1], spc[2] };
    Eigen::Vector3d omg_ibb{ ang[0], ang[1], ang[2] };
    fb /= dt, omg_ibb /= dt;
    auto C_bn = _ins.getAtt().toRotationMatrix();
    auto I_33 = Eigen::Matrix3d::Identity();
    const auto *pos = _ins.getPos();
    const auto *vel = _ins.getVel();

    F.block(0, 3, 3, 3) = I_33;
    F.block(9, 9, 3, 3) = -I_33 / _Tgb;
    F.block(12, 12, 3, 3) = -I_33 / _Tab;
    F.block(15, 15, 3, 3) = -I_33 / _Tgs;
    F.block(18, 18, 3, 3) = -I_33 / _Tas;

    F.block(6, 9, 3, 3) = -C_bn;
    F.block(3, 12, 3, 3) = C_bn;
    F.block(3, 6, 3, 3) = SkewSymmetric(C_bn * fb);
    F.block(6, 15, 3, 3) = -C_bn * omg_ibb.asDiagonal();
    F.block(3, 18, 3, 3) = C_bn * fb.asDiagonal();

    setF_rr(F, vel[0], vel[1], vel[2], pos[0], pos[2]);
    setF_vr(F, vel[0], vel[1], vel[2], pos[0], pos[2]);
    setF_vv(F, vel[0], vel[1], vel[2], pos[0], pos[2]);
    setF_phir(F, vel[0], vel[1], vel[2], pos[0], pos[2]);
    setF_phiv(F, pos[0], pos[2]);
    setF_Anginn(F, vel[0], vel[1], pos[0], pos[2]);
    return F;
}

GnssInsFilter::Hdim<6> GnssInsFilter::buildH(const ImuData &imu, double dt) const {
    using _Hdim = Hdim<6>;
    auto ang = imu.getGyro();
    Eigen::Vector3d omg_ibb{ ang[0], ang[1], ang[2] };
    omg_ibb /= dt;
    const auto &att = _ins.getAtt();
    const auto *pos = _ins.getPos();
    const auto *vel = _ins.getVel();
    _Hdim H{ _Hdim::Zero() };

    Eigen::Matrix<double, 3, 21> Hr = Eigen::Matrix<double, 3, 21>::Zero();
    Eigen::Vector3d lb{ _lever[0], _lever[1], _lever[2] };
    Hr.block(0, 0, 3, 3) = Eigen::Matrix<double, 3, 3>::Identity();
    Hr.block(0, 6, 3, 3) = SkewSymmetric(att.toRotationMatrix() * lb);
    H.block(0, 0, 3, 21) = Hr;

    Eigen::Matrix<double, 3, 21> Hv = Eigen::Matrix<double, 3, 21>::Zero();
    double R_M = getRM(pos[0]);
    double R_N = getRN(pos[0]);
    double sinphi = sin(pos[0]);
    double cosphi = cos(pos[0]);
    double tanphi = tan(pos[0]);
    Eigen::Matrix3d C_bn = att.toRotationMatrix();
    Eigen::Vector3d Ang_inn = getAng_innVec(pos[0], pos[2], vel[0], vel[1]);
    Eigen::Matrix3d H_v3 = -(SkewSymmetric(Ang_inn) * SkewSymmetric(C_bn * lb)) -
                           SkewSymmetric(C_bn * lb.cross(omg_ibb));
    Eigen::Matrix3d H_v6 = -C_bn * SkewSymmetric(lb) * omg_ibb.asDiagonal();
    Hv.block(0, 3, 3, 3) = Eigen::Matrix3d::Identity();
    Hv.block(0, 6, 3, 3) = H_v3;
    Hv.block(0, 9, 3, 3) = -C_bn * SkewSymmetric(lb);
    Hv.block(0, 15, 3, 3) = H_v6;
    H.block(3, 0, 3, 21) = Hv;

    return H;
}

GnssInsFilter::Qdim GnssInsFilter::buildQ(const GnssInsFilter::Fdim &F, double dt) const {
    using qdim = Eigen::Matrix<double, 18, 18>;
    using Gdim = Eigen::Matrix<double, 21, 18>;
    GnssInsFilter::Fdim Phi = Fdim::Identity() + F * dt;
    qdim q = qdim::Zero();

    Eigen::Matrix3d I_33 = Eigen::Matrix3d::Identity();
    q.block(0, 0, 3, 3) = _VRW * _VRW * I_33;
    q.block(3, 3, 3, 3) = _ARW * _ARW * I_33;
    Eigen::Matrix3d std_gb = Eigen::Vector3d(_std_gb[0], _std_gb[1], _std_gb[2]).asDiagonal();
    Eigen::Matrix3d std_ab = Eigen::Vector3d(_std_ab[0], _std_ab[1], _std_ab[2]).asDiagonal();
    Eigen::Matrix3d std_gs = Eigen::Vector3d(_std_gs[0], _std_gs[1], _std_gs[2]).asDiagonal();
    Eigen::Matrix3d std_as = Eigen::Vector3d(_std_as[0], _std_as[1], _std_as[2]).asDiagonal();
    q.block(6, 6, 3, 3) = 2 * std_gb * std_gb / _Tgb;
    q.block(9, 9, 3, 3) = 2 * std_ab * std_ab / _Tab;
    q.block(12, 12, 3, 3) = 2 * std_gs * std_gs / _Tgs;
    q.block(15, 15, 3, 3) = 2 * std_as * std_as / _Tas;

    Gdim G = Gdim::Zero();
    Eigen::Matrix3d C_bn = _ins.getAtt().toRotationMatrix();
    G.block(3, 0, 3, 3) = C_bn;
    G.block(6, 3, 3, 3) = C_bn;
    G.block(9, 6, 3, 3) = I_33;
    G.block(12, 9, 3, 3) = I_33;
    G.block(15, 12, 3, 3) = I_33;
    G.block(18, 15, 3, 3) = I_33;

    GnssInsFilter::Qdim Q = (Phi * G * q * G.transpose() * Phi.transpose() + G * q * G.transpose()) * 0.5 * dt;
    return Q;
}

void GnssInsFilter::setRandomWalk(double ARW, double VRW) {
    _ARW = ARW;
    _VRW = VRW;
}

GnssInsFilter::zdim<6> GnssInsFilter::buildz(const ImuData &imu, const GnssData &gnss, double dt) const {
    using _zdim = zdim<6>;
    Eigen::Vector3d lb{ _lever[0], _lever[1], _lever[2] };
    _zdim z{ _zdim::Zero() };
    auto ang = imu.getGyro();
    const auto *i_pos = _ins.getPos();
    const auto *g_pos = gnss.getBlh();
    const auto *i_vel = _ins.getVel();
    const auto *g_vel = gnss.getVel();
    const auto &Cb_n = _ins.getAtt();
    Eigen::Vector3d omg_ibb{ ang[0], ang[1], ang[2] };
    omg_ibb /= dt;
    // position observation
    Eigen::Vector3d r_G(g_pos[0], g_pos[1], g_pos[2]);
    Eigen::Matrix3d invD_R = getInv_DR(i_pos[0], i_pos[2]);
    Eigen::Vector3d hatr_I(i_pos[0], i_pos[1], i_pos[2]);
    Eigen::Vector3d hatr_G = hatr_I + invD_R * Cb_n * lb;
    Eigen::Vector3d z_r = invD_R.inverse() * (hatr_G - r_G);
    // velocity observation
    Eigen::Vector3d v_G(g_vel[0], g_vel[1], g_vel[2]);
    Eigen::Vector3d hatv_I(i_vel[0], i_vel[1], i_vel[2]);
    Eigen::Vector3d omg_inn = getAng_ienVec(i_pos[0]) + getAng_ennVec(i_pos[0], i_pos[2], i_vel[0], i_vel[1]);
    Eigen::Vector3d hatv_G = hatv_I - SkewSymmetric(omg_inn) * Cb_n * lb - Cb_n * lb.cross(omg_ibb);
    Eigen::Vector3d z_v = hatv_G - v_G;
    z << z_r, z_v;
    return z;
}

GnssInsFilter::Rdim<6> GnssInsFilter::buildR(const GnssData &gnss) const {
    using _Rdim = Rdim<6>;
    const auto *blh_std = gnss.getBlhStd();
    const auto *vel_std = gnss.getVelStd();
    Eigen::Vector<double, 6> R1{ pow(blh_std[0], 2), pow(blh_std[1], 2), pow(blh_std[2], 2),
                                 pow(vel_std[0], 2), pow(vel_std[1], 2), pow(vel_std[2], 2) };
    _Rdim R = R1.asDiagonal();
    return R;
}

void GnssInsFilter::ProcessData(const ImuData &imu, const GnssData &gnss) {
    auto curr_imu = imu;
    double dt = _ins.getTimeInterval(curr_imu.getSecond());
    auto prev_euler = QuaternionToEulerAngle(_ins.getAtt());
    curr_imu.Compensate(_gb, _ab, _gs, _as);
    _ins.INSUpdate(curr_imu);
    auto F = this->buildF(curr_imu, dt);
    auto Q = this->buildQ(F, dt);
    this->KFPredict(F, Q, dt);
    this->setX(GnssInsFilter::xdim::Zero());
    if (fabs(imu.getSecond() - gnss.getSecond()) < 1.0E-9) { // GNSS
        auto z = this->buildz(imu, gnss, dt);
        auto H = this->buildH(imu, dt);
        auto R = this->buildR(gnss);
        this->KFUpdate<6>(z, H, R, dt);
    }
    if (_zupt.WindowMoveOn(curr_imu), _zupt.isZeroUpdate()) { // ZUPT
        auto curr_euler = QuaternionToEulerAngle(_ins.getAtt());
        auto z = this->ZeroBuildz(curr_imu, prev_euler, curr_euler, dt);
        auto H = this->ZeroBuildH(curr_imu, curr_euler, dt);
        auto R = this->ZeroBuildR();
        this->KFUpdate<7>(z, H, R, dt);
    }
    // Correct position, velocity, attitude and device parameter
    const auto *pos = _ins.getPos();
    double R_M = getRM(pos[0]);
    double R_N = getRN(pos[0]);
    Eigen::Vector3d dr{ _x.segment(0, 3) }; // delta ned
    Eigen::Vector3d dv{ _x.segment(3, 3) };
    Eigen::Vector3d dphi{ _x.segment(6, 3) }; // delta rotate-vector
    Eigen::Matrix3d invD_R = getInv_DR(pos[0], pos[2]);
    Eigen::Vector3d dp = invD_R * dr; // delta blh
    Eigen::Quaterniond dq = RotateVectorToQuaternion(dphi); // delta quaternion
    _ins.Correct(dp, dv, dq);
    Eigen::Vector3d dgb = _x.segment(9, 3);
    Eigen::Vector3d dab = _x.segment(12, 3);
    Eigen::Vector3d dgs = _x.segment(15, 3);
    Eigen::Vector3d das = _x.segment(18, 3);
    for (int i = 0; i < 3; ++i) {
        _gb[i] += dgb(i) * dt;
        _ab[i] += dab(i) * dt;
        _gs[i] += dgs(i);
        _as[i] += das(i);
    }
    this->setX(GnssInsFilter::xdim::Zero());
}

void GnssInsFilter::PrintState() const {
    _ins.PrintState();
}

std::tuple<const double *, const double *, Eigen::Quaterniond> GnssInsFilter::getState() const {
    return _ins.getState();
}

bool GnssInsFilter::isZeroUpdate() const {
    return _zupt.isZeroUpdate();
}

GnssInsFilter::zdim<7>
GnssInsFilter::ZeroBuildz(const ImuData &imu, const EulerAngle &prev_euler, const EulerAngle &curr_euler, double dt) const {
    using _zdim = GnssInsFilter::zdim<7>;
    const auto *vel = _ins.getVel();
    const auto &gyro = imu.getGyro();
    double yaw_rate = (curr_euler.yaw - prev_euler.yaw) / dt;
    _zdim z;
    z << vel[0], vel[1], vel[2], yaw_rate, gyro(0), gyro(1), gyro(2);
    return z;
}

GnssInsFilter::Hdim<7> GnssInsFilter::ZeroBuildH(const ImuData &imu, const EulerAngle &euler, double dt) const {
    using _Hdim = GnssInsFilter::Hdim<7>;
    _Hdim H = _Hdim::Zero();
    Eigen::Matrix<double, 1, 3> H_psi;
    H_psi << 0, sin(euler.roll) / cos(euler.pitch), cos(euler.roll) / cos(euler.pitch);
    H.block(0, 3, 3, 3) = Eigen::Matrix3d::Identity();
    H.block(3, 9, 1, 3) = H_psi;
    H.block(4, 9, 3, 3) = Eigen::Matrix3d::Identity();
    H.block(4, 15, 3, 3) = Eigen::Matrix3d::Identity();
    return H;
}

GnssInsFilter::Rdim<7> GnssInsFilter::ZeroBuildR() const {
    return GnssInsFilter::Rdim<7>::Identity();
}
