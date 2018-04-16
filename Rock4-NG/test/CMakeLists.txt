cmake_minimum_required(VERSION 2.6)
project(Odes-NG)
enable_language(CXX)
#
# To use intel compiler:
# CXX=icpc cmake ..
# otherwise, to use g++:
#  cmake ..
#
include_directories(
${CMAKE_SOURCE_DIR}/../include
${CMAKE_SOURCE_DIR}/../../common/include
)

if (${CMAKE_CXX_COMPILER} MATCHES "icpc.*$")
  set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS}  -DICC -DALIGN_BLOCKED -restrict -O3   -g -xavx -fargument-noalias  -ansi-alias -qopt-report-phase=vec -qopt-report-file=stdout -qopt-subscript-in-range -std=c++14 -qopenmp")


else ()
  set (USING_GNU TRUE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -Wall  -O3  -std=c++14 -g -DGCC -DALIGN_BLOCKED  -march=native  -fopenmp ")

endif ()

add_executable(
  run
  ../main.cc
  )
if (${CMAKE_CXX_COMPILER} MATCHES "icpc.*$")
  target_link_libraries(
    run
    mkl_intel_lp64
    mkl_sequential
    mkl_core
    mkl_lapack95_lp64
    mkl_blas95_lp64
    tbb
    tbbmalloc_proxy
    tbbmalloc
    rt
    )
else ()
  target_link_libraries(
    run
    lapack 
    )
endif ()



