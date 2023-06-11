//
// Created by 崔宇璐 on 2023/6/1.
//

#ifndef GUPOS_KALMANFILTER_HPP
#define GUPOS_KALMANFILTER_HPP

#include <Eigen/Dense>

template<int State, int Measure>
class KalmanFilter {
public:
    using xdim = Eigen::Vector<double, State>;
    using Pdim = Eigen::Matrix<double, State, State>;
    using Fdim = Eigen::Matrix<double, State, State>;
    using Qdim = Eigen::Matrix<double, State, State>;
    using zdim = Eigen::Vector<double, Measure>;
    using Hdim = Eigen::Matrix<double, Measure, State>;
    using Rdim = Eigen::Matrix<double, Measure, Measure>;
public:
    void Initialize(const xdim& x, const Pdim& P);
    void KFPredict(const Fdim &F, const Qdim &Q, double dt);
    void KFUpdate(const zdim &z, const Hdim &H, const Rdim &R, double dt);
public: // Getter and Setter
    const xdim &getX() const {
        return _x;
    }

    void setX(const xdim &x) {
        _x = x;
    }

protected:
    xdim _x{ xdim::Zero() };
    Pdim _P{ Pdim::Zero() };
};

template<int State, int Measure>
void KalmanFilter<State, Measure>::KFUpdate(const KalmanFilter::zdim &z, const KalmanFilter::Hdim &H,
                                            const KalmanFilter::Rdim &R, double dt) {
    using Kdim = Eigen::Matrix<double, State, Measure>;
    Kdim K = _P * H.transpose() * (H * _P * H.transpose() + R).inverse();
    Pdim I = Pdim::Identity();
    _x = _x + K * (z - H * _x);
    _P = (I - K * H) * _P * (I - K * H).transpose() + K * R * K.transpose();
}

template<int State, int Measure>
void KalmanFilter<State, Measure>::KFPredict(const KalmanFilter::Fdim &F, const KalmanFilter::Qdim &Q, double dt) {
    Fdim Phi = Fdim::Identity() + F * dt;
    _x = Phi * _x;
    _P = Phi * _P * Phi.transpose() + Q;
}

template<int State, int Measure>
void KalmanFilter<State, Measure>::Initialize(const KalmanFilter::xdim &x, const KalmanFilter::Pdim &P) {
    _x = x;
    _P = P;
}

#endif //GUPOS_KALMANFILTER_HPP
