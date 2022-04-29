

file(GLOB_RECURSE PostProcessFiles *.cpp *.c)
find_package(VTK PATHS ${VTK_DIR})

if (NOT VTK_FOUND)
    message(FATAL_ERROR "cannot find vtk libraries: ${VTK_NOT_FOUND_MESSAGE}")
    return()
endif()
message(STATUS "VTK_VERSION: ${VTK_VERSION}")
#message(STATUS "POSTPROCESS/CMAKE_TYPE = ${CMAKE_BUILD_TYPE}")


#add_compile_options(-fPIC)

include_directories(${PROJECT_SOURCE_DIR}/src)
add_executable(PostProcess ${PostProcessFiles})
target_link_libraries(PostProcess ${VTK_LIBRARIES})
# vtk_module_autoinit is needed
vtk_module_autoinit(TARGETS PostProcess MODULES ${VTK_LIBRARIES})
target_link_libraries(PostProcess InputParser)

