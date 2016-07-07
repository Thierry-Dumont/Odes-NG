cmake_minimum_required(VERSION 2.6)
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
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS}  -DICC -DALIGN_BLOCKED -restrict -O3  -g -xavx -fargument-noalias  -ansi-alias  -opt-subscript-in-range  -vec-report=0")
include_directories(
${CMAKE_SOURCE_DIR}/.
$ENV{HOME}/odes/trunk/common/include
$ENV{HOME}/odes/trunk/Radau5-NG/include
$ENV{HOME}/odes/trunk/Rock4-L/include
)
elseif(${CMAKE_CXX_COMPILER} MATCHES "clang.*$")
  set (USING_GNU TRUE)
  set(CMAKE_CXX_COMPILER "clang++")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall  -DCLANG")
else ()
  set (USING_GNU TRUE)
  set(CMAKE_CXX_COMPILER "g++")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}   -g -O3 -std=c++14 -Wall  -DGCC  -march=native -ftree-vectorize  -fstrict-aliasing -ffast-math")
include_directories(
${CMAKE_SOURCE_DIR}/.
$ENV{HOME}/odes/trunk/common/include
$ENV{HOME}/odes/trunk/Radau5-NG/include
$ENV{HOME}/odes/trunk/Rock4-L/include
)

endif ()

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



