//
// Created by admin on 2023/6/12.
//

#ifndef GUPOS_ZERODETECTOR_H
#define GUPOS_ZERODETECTOR_H

#include <deque>
#include "ImuData.h"

class ZeroDetector {
public:
    ZeroDetector() = delete;
    ZeroDetector(int window, double threshold)
    : _window(window), _threshold(threshold)
    {}
    void WindowMoveOn(const ImuData &imu);
    bool IsZeroUpdate() const;
private:
    int _window{ 0 };
    double _threshold{ 0.0 };
    std::deque<ImuData> _detector;
};


#endif //GUPOS_ZERODETECTOR_H
