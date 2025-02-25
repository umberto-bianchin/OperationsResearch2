cmake_minimum_required(VERSION 3.10)
project(Homework01)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
set(CMAKE_BINARY_DIR "${CMAKE_SOURCE_DIR}/build")
set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")

set(CPLEXDIR /Applications/CPLEX_Studio2211/cplex)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -I. -I${CPLEXDIR}/include/ilcplex")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L${CPLEXDIR}/lib/arm64_osx/static_pic -L. -lcplex -lpthread -ldl")

include_directories(include)

add_executable(Homework01
    src/main.c
    src/chrono.c
    src/tsp.c
)