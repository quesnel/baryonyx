cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(baryonyx VERSION 0.2.99.0 LANGUAGES CXX C)

find_package(Threads REQUIRED)
add_library(threads INTERFACE IMPORTED)
set_property(TARGET threads PROPERTY
  INTERFACE_LINK_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})

add_definitions(-DFMT_HEADER_ONLY)
include_directories(${CMAKE_SOURCE_DIR}/external/fmt/include)

enable_testing()
add_subdirectory(lib)
add_subdirectory(app)
