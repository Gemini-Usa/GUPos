cmake_minimum_required(VERSION 3.13)  # CMake version check
project(GUPos)                        # Create project "simple_example"
set(CMAKE_CXX_STANDARD 17)            # Enable c++14 standard

if (WIN32)
    add_compile_options(-bigobj)
    include_directories("D:/eigen-master/")
elseif (APPLE)
    include_directories("/usr/local/Cellar/eigen/3.4.0_1/include/eigen3/")
endif()

# Add main.cpp file of project root directory as source file
set(SOURCE_FILES
        main.cpp
        Utility/Algebra.cpp
        Utility/Algebra.h
        Utility/Geodesy.cpp
        Utility/Geodesy.h
        GUCore/GnssData.cpp
        GUCore/GnssData.h
        GUCore/GnssInsFilter.cpp
        GUCore/GnssInsFilter.h
        GUCore/ImuData.cpp
        GUCore/ImuData.h
        GUCore/InsMechanizer.cpp
        GUCore/InsMechanizer.h
        GUCore/KalmanFilter.hpp
        Configure/minIni.c
        Configure/minIni.h
        Configure/minGlue.h Configure/GuConfig.cpp Configure/GuConfig.h GUCore/ZeroDetector.cpp GUCore/ZeroDetector.h Utility/String.cpp Utility/String.h)

# Add executable target with source files listed in SOURCE_FILES variable
add_executable(GUPos ${SOURCE_FILES})
