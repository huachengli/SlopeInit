cmake_minimum_required(VERSION 3.19)
project(SlopeInit)

set(CMAKE_CXX_STANDARD 20)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_FLAGS_DEBUG "-O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

#set dealii path
#if(NOT DEFINED ${DEAL_II_CAME_PATH})
#    message(FATAL_ERROR "DEAL_II_CAME_PATH is unknown")
#endif()
#
##set vtk path
#if(NOT DEFINED ${VTK_DIR})
#    message(FATAL_ERROR "VTK_DIR is unknown")
#endif()

message(STATUS "build INFORMATION:")
message(STATUS "CXX_FLAGS_RE = ${CMAKE_CXX_FLAGS_RELEASE} ")
message(STATUS "DEAL_II_CAME_PATH = ${DEAL_II_CAME_PATH}")
message(STATUS "VTK_DIR = ${VTK_DIR}")

#set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

include_directories(${PROJECT_SOURCE_DIR}/src})
add_subdirectory(src)
add_subdirectory(postprocess)
