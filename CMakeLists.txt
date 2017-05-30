cmake_minimum_required(VERSION 2.8.1)
project(MyDGtalContrib)

find_package(DGtal REQUIRED)
include_directories(${DGTAL_INCLUDE_DIRS})
link_directories(${DGTAL_LIBRARY_DIRS})

set(EXAMPLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})

FIND_PACKAGE(Boost 1.46.0 REQUIRED COMPONENTS)
if ( Boost_FOUND )
  message(STATUS "Boost and boost_program_options found.")
  include_directories( ${Boost_INCLUDE_DIRS} )
endif( Boost_FOUND )

if (NOT WITH_EIGEN)
  message(FATAL_ERROR "Please compile DGtal with Eigen")
endif()

if (NOT WITH_VISU3D_QGLVIEWER)
  message(FATAL_ERROR "Please compile DGtal with QGLViewer")
endif()

include_directories(${CMAKE_SOURCE_DIR}/src/)
option(BUILD_TESTS "Build test programs" OFF)
option(BUILD_EXAMPLES "Build examples programs" ON)
option(WITH_OMP "Parallelization with OpenMP" OFF)

if (WITH_OMP)
  find_package(OpenMP)
  if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  endif()
endif()



if (BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()


message(STATUS ${CMAKE_CXX_FLAGS})

