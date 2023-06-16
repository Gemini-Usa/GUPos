//
// Created by GeminiUsa on 2023/6/12.
//

#include "ZeroDetector.h"

void ZeroDetector::WindowMoveOn(const ImuData &imu) {
    if (_detector.size() < _window) {
        _detector.push_back(imu);
    } else {
        _detector.pop_front();
        _detector.push_back(imu);
    }
}

bool ZeroDetector::IsZeroUpdate() const {
    double sum = 0.0;
    for (const auto &imu : _detector) {
        sum += pow(imu.GetGyro().norm() * ImuData::GetRate(), 2);
    }
    if ((sum / static_cast<double>(_detector.size())) < _threshold) return true;
    else return false;
}
