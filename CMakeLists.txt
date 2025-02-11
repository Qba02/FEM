cmake_minimum_required(VERSION 3.24)
project(Projekt_MES)


set(CMAKE_CXX_STANDARD 17)

set(SOURCE_FILES main.cpp src/grid.cpp include/grid.h
        include/GaussIntegral.h src/GaussIntegral.cpp
        include/UniElement.h src/UniElement.cpp
        src/matrix.cpp include/matrix.h
        src/SoE.cpp
        include/SoE.h)
add_executable(Projekt_MES ${SOURCE_FILES})
