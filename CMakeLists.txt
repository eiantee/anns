cmake_minimum_required(VERSION 3.28.3)
project(work)

set(CMAKE_CXX_STANDARD 17)

add_executable(work main.cpp
        graph.h
        knn.h
        nsg.h
        kmean.h
)
