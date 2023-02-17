cmake_minimum_required(VERSION 3.24)
project(tspm_cpp_backend)

set(CMAKE_CXX_STANDARD 14)

include_directories(utils)

add_executable(tspm_cpp_backend
        utils/utils.cpp
        utils/utils.h
        main.cpp)
