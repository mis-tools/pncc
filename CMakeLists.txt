cmake_minimum_required(VERSION 2.8.12)
 
project(ParallelNormalizedCrossCorrelation)
 
FIND_PACKAGE( OpenMP )
if(OPENMP_FOUND)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

add_executable(pncc src/pncc.cxx)
INSTALL(TARGETS pncc DESTINATION bin/)
