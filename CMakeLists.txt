cmake_minimum_required(VERSION 3.24)
project(tspm_cpp_backend)

set(CMAKE_CXX_STANDARD 17)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -lpthread -std=c++17 -fopenmp -latomic ")
include_directories(utils)
include_directories(example)
include_directories(lib/ips4o)

add_executable(tspm_cpp_backend
        utils/utils.cpp
        utils/utils.h
        utils/sequencing.cpp
        utils/sequencing.h
        utils/dbMartEntry.h
        utils/temporalSequence.h utils/sorter.cpp utils/sorter.h utils/workflows.cpp utils/workflows.h example/example_tSPM_calls.cpp tspmPlus.h)
