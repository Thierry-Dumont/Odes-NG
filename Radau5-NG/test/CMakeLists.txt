cmake_minimum_required(VERSION 3.5)
project(NewRadau5)
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
  set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS}  -DICC -DALIGN_BLOCKED -restrict -O3   -g -xavx -fargument-noalias  -ansi-alias -qopt-report-phase=vec -qopt-report-file=stdout -qopt-subscript-in-range -std=c++14 ")
elseif(${CMAKE_CXX_COMPILER} MATCHES "clang++.*$")
  set (USING_GNU TRUE)
  set(CMAKE_CXX_COMPILER "clang++")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall  -DCLANG")
else ()
  set (USING_GNU TRUE)
  #set(CMAKE_CXX_COMPILER "g++-4.9")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -O3 -g -Wall -std=c++11  -DALIGN_BLOCKED -DGCC  -march=native -ftree-vectorize  -fstrict-aliasing -ffast-math")
#-fopt-info-vec-missed  -fopt-info-vec -ffast-math
endif ()
include_directories(
${CMAKE_SOURCE_DIR}/../include
${CMAKE_SOURCE_DIR}/../../common/include
)
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
    rt
    )
else ()
  target_link_libraries(
    run
    lapack 
    )
endif ()



