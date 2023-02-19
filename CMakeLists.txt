cmake_minimum_required(VERSION 3.24)
project(tspm_cpp_backend)

set(CMAKE_CXX_STANDARD 14)

# added -fopenmp
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -fopenmp")
include_directories(utils)

add_executable(tspm_cpp_backend
        utils/utils.cpp
        utils/utils.h
        main.cpp)
