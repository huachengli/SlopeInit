cmake_minimum_required(VERSION 3.19)
project(SlopeInit)

set(CMAKE_CXX_STANDARD 20)

set(SRCS src/InputParser.c src/SlopeGeometry.c src/eos_state.c src/eos_cpp.h src/eos_cpp.cpp src/SlopeSlove.cpp)

set(DEAL_II_CAME_PATH "/home/huacheng/hlib/dealii")
FIND_PACKAGE(deal.II 10.0.0 REQUIRED PATHS ${DEAL_II_CAME_PATH})
DEAL_II_INITIALIZE_CACHED_VARIABLES()

ADD_EXECUTABLE(SlopeInit main.cpp ${SRCS})
DEAL_II_SETUP_TARGET(SlopeInit)
