cmake_minimum_required(VERSION 3.10)
project(Homework01)


set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
set(CMAKE_BINARY_DIR "${CMAKE_SOURCE_DIR}/build")
set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")
    
if(APPLE)
    set(CPLEXDIR /Applications/CPLEX_Studio2211/cplex)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -I. -I${CPLEXDIR}/include/ilcplex")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L${CPLEXDIR}/lib/arm64_osx/static_pic -L. -lcplex -lpthread -ldl")
else()
    set(CPLEXDIR "C:\\Program Files\\IBM\\ILOG\\CPLEX_Studio2212\\cplex")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /I. /I\"${CPLEXDIR}\\include\\ilcplex\"")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /LIBPATH:\"${CPLEXDIR}\\lib\\x64_windows_msvc14\\stat_mda\" /LIBPATH:\"C:\\Program Files\\IBM\\ILOG\\CPLEX_Studio2212\\concert\\lib\\x64_windows_msvc14\\stat_mda\" cplex2212.lib ilocplex.lib")
endif()

include_directories(include)

add_executable(Homework01
    src/main.c
    src/chrono.c
    src/tsp_utilities.c
    src/parsers.c
)

if(!APPLE)
    target_include_directories(Homework01 PRIVATE "C:/Users/fraf/Documents/gnuplot_c/gnuplot_c/src")
endif()