cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(baryonyx VERSION 0.4.0 LANGUAGES CXX)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

find_package(Threads REQUIRED)
add_library(threads INTERFACE IMPORTED)
set_property(TARGET threads PROPERTY
  INTERFACE_LINK_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})

find_package(NLOPT)
if (NLOPT_FOUND)
  add_library(nlopt INTERFACE IMPORTED)
  set_property(TARGET nlopt PROPERTY
    INTERFACE_LINK_LIBRARIES ${NLOPT_LIBRARIES})
endif ()

add_definitions(-DFMT_HEADER_ONLY)
include_directories(${CMAKE_SOURCE_DIR}/external/fmt/include)

find_program(CLANG_TIDY_EXE
  NAMES "clang-tidy"
  DOC "Path to clang-tidy executable")

if (NOT CLANG_TIDY_EXE)
  message(STATUS "clang-tidy not found")
else ()
  message(STATUS "clang-tidy found: ${CLANG_TIDY_EXE}")
  set (DO_CLANG_TIDY "${CLANG_TIDY_EXE}"
    "-checks=-*,cppcoreguidelines-*,readibility-*,modernize-*,performance-*,bugprone-*")
endif ()

enable_testing()
add_subdirectory(lib)
add_subdirectory(app)
