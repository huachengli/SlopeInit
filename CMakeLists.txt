cmake_minimum_required(VERSION 3.19)
project(SlopeInit)

set(CMAKE_CXX_STANDARD 20)

set(SRCS InputParser.c SlopeGeometry.c SlopeGeometry.h)

add_executable(SlopeInit main.cpp ${SRCS})