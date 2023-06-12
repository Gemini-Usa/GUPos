cmake_minimum_required(VERSION 3.13)  # CMake version check
project(GUPos)                        # Create project "simple_example"
set(CMAKE_CXX_STANDARD 17)            # Enable c++14 standard

include_directories("/usr/local/Cellar/eigen/3.4.0_1/include/eigen3/")

# Add main.cpp file of project root directory as source file
set(SOURCE_FILES
        main.cpp
        Algebra.cpp
        Algebra.h
        Geodesy.cpp
        Geodesy.h
        GnssData.cpp
        GnssData.h
        GnssInsFilter.cpp
        GnssInsFilter.h
        ImuData.cpp
        ImuData.h
        InsMechanizer.cpp
        InsMechanizer.h
        KalmanFilter.hpp
        minIni.c
        minIni.h
        minGlue.h GuConfig.cpp GuConfig.h)

# Add executable target with source files listed in SOURCE_FILES variable
add_executable(GUPos ${SOURCE_FILES})
