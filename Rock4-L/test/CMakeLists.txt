cmake_minimum_required(VERSION 2.8)
project(Radau5-NG)
enable_language(CXX)
#
# To use intel compiler:
# CXX=icpc  cmake ..
# to use clang:
# CXX=clang cmake ..
# otherwise, to use g++:
#  cmake ..
#
if (${CMAKE_CXX_COMPILER} MATCHES "icpc.*$")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -DICC -DALIGN_BLOCKED -restrict -O3   -g -xavx -fargument-noalias  -ansi-alias -qopt-report-phase=vec -qopt-report-file=stdout -qopt-subscript-in-range -std=c++14 ")

else ()
  set (USING_GNU TRUE)
  set(CMAKE_CXX_COMPILER "g++")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}   -g -O3 -std=c++14 -Wall  -DGCC  -march=native -ftree-vectorize  -fstrict-aliasing -ffast-math")


endif ()

include_directories(
${CMAKE_SOURCE_DIR}/../include
${CMAKE_SOURCE_DIR}/../../common/include
${CMAKE_SOURCE_DIR}/../../Radau5-NG/include/
)
add_executable(
  run
  ../main.cpp
   )
if (${CMAKE_CXX_COMPILER} MATCHES "icpc.*$")
  target_link_libraries(
    run
    mkl_intel_lp64
    mkl_sequential
    mkl_core
    mkl_lapack95_lp64
    mkl_blas95_lp64
    rt
    )
else ()
  target_link_libraries(
    run
    lapack

    )
endif ()



