cmake_minimum_required(VERSION 3.19)
project(SlopeInit)

set(CMAKE_CXX_STANDARD 20)

set(SRCS src/InputParser.c src/SlopeGeometry.c src/eos_state.c src/eos_cpp.h src/eos_cpp.cpp)

add_executable(SlopeInit main.cpp ${SRCS})