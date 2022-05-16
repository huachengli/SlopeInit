#CMAKELISTS for SlopInit process
file(GLOB_RECURSE SlopInitFiles *.cpp *.c)
link_directories(${PROJECT_BINARY_DIR}/lib)
FIND_PACKAGE(deal.II 10.0.0 PATHS ${DEAL_II_CAME_PATH})

if(deal.II_FOUND)
    DEAL_II_INITIALIZE_CACHED_VARIABLES()
    ADD_EXECUTABLE(SlopeInit main.cpp ${SlopInitFiles})
    DEAL_II_SETUP_TARGET(SlopeInit)
else()
    MESSAGE(STATUS "DEALII not found, SlopeInit wouldn't be compiled")
endif()


add_compile_options(-fPIC)
set(InputParserFile InputParser.c InputParser.h)
add_library(InputParser SHARED ${InputParserFile})