cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(baryonyx VERSION 0.1 LANGUAGES C CXX)

find_package(Threads REQUIRED)
add_library(threads INTERFACE IMPORTED)
set_property(TARGET threads PROPERTY
  INTERFACE_LINK_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})

enable_testing()
add_subdirectory(lib)
add_subdirectory(app)
