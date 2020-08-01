cmake_minimum_required(VERSION 3.0.0)
project("Ray Tracing")

set(CMAKE_CXX_FLAGS "--std=c++17 -O3 -Wall -Wextra -fopenmp") 
include_directories(cpplap cpplap/test)

add_executable(dmd_test dmd_test.cpp)
add_executable(raytracing_test raytracing_test.cpp)
add_executable(run_raytracing run_raytracing.cpp)
