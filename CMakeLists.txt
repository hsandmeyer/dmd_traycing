cmake_minimum_required(VERSION 3.0.0)
project("Ray Tracing")

set(CMAKE_CXX_FLAGS "--std=c++17 -O3 -Wall -Wextra -fopenmp") 
include_directories(cpplap cpplap/test src)

add_executable(dmd_test src/dmd_test.cpp)
add_executable(raytracing_test src/raytracing_test.cpp)
add_executable(run_raytracing src/run_raytracing.cpp)
