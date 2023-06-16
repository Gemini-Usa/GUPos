//
// Created by GeminiUsa on 2023/6/1.
//

#ifndef GUPOS_KALMANFILTER_HPP
#define GUPOS_KALMANFILTER_HPP

#include <Eigen/Dense>

template<int State>
class KalmanFilter {
public:
    using xdim = Eigen::Vector<double, State>;
    using Pdim = Eigen::Matrix<double, State, State>;
    using Fdim = Eigen::Matrix<double, State, State>;
    using Qdim = Eigen::Matrix<double, State, State>;
    template<int Measure> using zdim = Eigen::Vector<double, Measure>;
    template<int Measure> using Hdim = Eigen::Matrix<double, Measure, State>;
    template<int Measure> using Rdim = Eigen::Matrix<double, Measure, Measure>;
public:
    void Initialize(const xdim& x, const Pdim& P);
    void KFPredict(const Fdim &F, const Qdim &Q, double dt);
    template<int Measure>
    void KFUpdate(const zdim<Measure> &z, const Hdim<Measure> &H, const Rdim<Measure> &R, double dt);
public: // Getter and Setter
    const xdim &GetX() const {
        return _x;
    }
    void SetX(const xdim &x) {
        _x = x;
    }
protected:
    xdim _x{ xdim::Zero() };
    Pdim _P{ Pdim::Zero() };
};

template<int State>
template<int Measure>
void KalmanFilter<State>::KFUpdate(const zdim<Measure> &z, const Hdim<Measure> &H, const Rdim<Measure> &R, double dt) {
    using Kdim = Eigen::Matrix<double, State, Measure>;
    Kdim K = _P * H.transpose() * (H * _P * H.transpose() + R).inverse();
    Pdim I = Pdim::Identity();
    _x = _x + K * (z - H * _x);
    _P = (I - K * H) * _P * (I - K * H).transpose() + K * R * K.transpose();
}

template<int State>
void KalmanFilter<State>::KFPredict(const KalmanFilter::Fdim &F, const KalmanFilter::Qdim &Q, double dt) {
    Fdim Phi = Fdim::Identity() + F * dt;
    _x = Phi * _x;
    _P = Phi * _P * Phi.transpose() + Q;
}

template<int State>
void KalmanFilter<State>::Initialize(const KalmanFilter::xdim &x, const KalmanFilter::Pdim &P) {
    _x = x;
    _P = P;
}

#endif //GUPOS_KALMANFILTER_HPP
