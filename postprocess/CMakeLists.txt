

file(GLOB_RECURSE PostProcessFiles *.cpp *.c)
find_package(VTK PATHS ${VTK_DIR})

if(VTK_FOUND)
    message(STATUS "VTK_VERSION: ${VTK_VERSION}")
    message(STATUS "POSTPROCESS/CMAKE_TYPE = ${CMAKE_BUILD_TYPE}")

    # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static")
    # set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -static")
    # set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")

    message(STATUS "POSTPROCESS/CMAKE_C_FLAGES = ${CMAKE_C_FLAGES}")
    include_directories(${PROJECT_SOURCE_DIR}/src)
    add_executable(PostProcess ${PostProcessFiles})
    target_link_libraries(PostProcess ${VTK_LIBRARIES})
    # vtk_module_autoinit is needed
    vtk_module_autoinit(TARGETS PostProcess MODULES ${VTK_LIBRARIES})
    target_link_libraries(PostProcess InputParser)
else()
    message(STATUS "cannot find vtk libraries: ${VTK_NOT_FOUND_MESSAGE}, postprocess will not be compiled.")
endif()


