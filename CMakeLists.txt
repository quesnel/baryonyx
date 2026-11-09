cmake_minimum_required(VERSION 3.20)
project(baryonyx VERSION 0.5.0 LANGUAGES CXX)

# Use the GNU standard installation directories. See
# https://cmake.org/cmake/help/v3.0/module/GNUInstallDirs.html
include(GNUInstallDirs)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

find_package(Threads REQUIRED)
add_library(threads INTERFACE IMPORTED)
set_property(TARGET threads PROPERTY
  INTERFACE_LINK_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})

find_package(NLopt CONFIG QUIET)
if (NLopt_FOUND)
    set(NLOPT_FOUND TRUE)
    set(NLOPT_LIBRARIES NLopt::nlopt)
    set(NLOPT_INCLUDE_DIRS "")
else ()
    find_package(NLOPT)
endif ()
message(STATUS "NLopt support: ${NLOPT_FOUND}")

add_definitions(-DFMT_HEADER_ONLY)
include_directories(${CMAKE_SOURCE_DIR}/external/fmt/include)
include_directories(${CMAKE_SOURCE_DIR}/external/ut/include)

enable_testing()
add_subdirectory(lib)
add_subdirectory(app)
