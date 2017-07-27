cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(lp VERSION 0.1 LANGUAGES C CXX)

find_package(Eigen3 REQUIRED)
add_library(eigen3 INTERFACE IMPORTED)
set_property(TARGET eigen3 PROPERTY
  INTERFACE_INCLUDE_DIRECTORIES ${EIGEN3_INCLUDE_DIR})

find_package(Threads REQUIRED)
add_library(threads INTERFACE IMPORTED)
set_property(TARGET threads PROPERTY
  INTERFACE_LINK_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})

enable_testing()
add_subdirectory(lib)
add_subdirectory(app)
