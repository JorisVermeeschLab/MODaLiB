cmake_minimum_required(VERSION 3.26)
project(MODaLiB)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(MODaLiB
        fragment.h
        main.cpp
        MethylDackel.h
        read_info.h
        serialization.h)
