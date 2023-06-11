cmake_minimum_required(VERSION 3.24)
project(GUPos)

set(CMAKE_CXX_STANDARD 17)

include_directories("/usr/local/Cellar/eigen/3.4.0_1/include/eigen3/")

add_executable(GUPos main.cpp ImuData.cpp ImuData.h Geodesy.cpp Geodesy.h InsMechanizer.cpp InsMechanizer.h Algebra.cpp Algebra.h GnssData.cpp GnssData.h KalmanFilter.hpp GnssInsFilter.cpp GnssInsFilter.h)
