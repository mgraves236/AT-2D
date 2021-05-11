cmake_minimum_required(VERSION 3.17)
project(AT_2D_CLASS)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-std=c++11 -O3")
set(CMAKE_EXE_LINKER_FLAGS "-static")

add_executable(AT_2D_CLASS main.cpp Node.cpp Node.h Lattice.cpp Lattice.h func.cpp func.h)