cmake_minimum_required(VERSION 3.7)

add_definitions(-std=c++17)

set(CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS, "${CXX_FLAGS}")

project(6DOF_AC_Sim)

set(SOURCE_FILES src/lib/matrixLibrary/Mat.cpp src/lib/numericalInegrator/ERK.cpp src/main.cpp src/Aircraft.cpp)

add_executable(6DOF_Sim_test ${SOURCE_FILES})