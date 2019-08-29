cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project(libbaryonyx VERSION 0.5.0 LANGUAGES CXX)

option(WITH_FULL_OPTIMIZATION "Disable all logging facilities and active heavy optimization code. [default: off]" OFF)

option(WITH_DEBUG "enable maximium debug code. [default: ON]" ON)
option(WITH_LOG "enable log message. [default: ON]" ON)

set(public_baryonyx_header
 include/baryonyx/core-compare
 include/baryonyx/core-out
 include/baryonyx/core-test
 include/baryonyx/core-utils)

set(private_baryonyx_source
 src/nlopt-optimizer.cpp
 src/manual-optimizer.cpp
 src/branch-optimizer.cpp
 src/branch-and-bound-solver.hpp
 src/debug.hpp
 src/bit-array.hpp
 src/fixed-2darray.hpp
 src/fixed-array.hpp
 src/itm-common.hpp
 src/itm-optimizer-common.hpp
 src/itm-solver-common.hpp
 src/itm-solver-equalities-01.cpp
 src/itm-solver-equalities-101.cpp
 src/itm-solver-inequalities-01.cpp
 src/itm-solver-inequalities-101-buffered.cpp
 src/itm-solver-inequalities-101.cpp
 src/itm-solver-inequalities-Z.cpp
 src/itm.hpp
 src/knapsack-dp-solver.hpp
 src/lpcore.cpp
 src/memory.cpp
 src/memory.hpp
 src/merged-constraint.cpp
 src/observer.hpp
 src/parser.cpp
 src/pnm.hpp
 src/preprocessor.cpp
 src/private.cpp
 src/private.hpp
 src/problem.cpp
 src/problem.hpp
 src/random-solver.cpp
 src/result.cpp
 src/result.hpp
 src/resume.hpp
 src/sparse-matrix.hpp
 src/utils.hpp)

add_library(libbaryonyx ${public_baryonyx_header} ${private_baryonyx_source})

option(BUILD_SHARED_LIBS "Build shared library" OFF)
include(GenerateExportHeader)
generate_export_header(libbaryonyx
    EXPORT_MACRO_NAME BX_EXPORT
    NO_EXPORT_MACRO_NAME BX_NO_EXPORT
    EXPORT_FILE_NAME ${CMAKE_CURRENT_BINARY_DIR}/baryonyx/export)

target_compile_options(libbaryonyx
  PRIVATE
  $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:
      -pipe -march=native
      $<$<CONFIG:Debug>:-O0 -Wall -Wextra>>
  $<$<CXX_COMPILER_ID:MSVC>:
      $<$<CONFIG:Debug>:/Od /W3 /Zi>>
  $<$<CXX_COMPILER_ID:MSVC>:/bigobj>)

target_include_directories(libbaryonyx PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/>)

target_link_libraries(libbaryonyx
  PRIVATE
  threads
  $<$<BOOL:${NLOPT_FOUND}>:nlopt>)

set_target_properties(libbaryonyx PROPERTIES
  OUTPUT_NAME "baryonyx-${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}"
  POSITION_INDEPENDENT_CODE ON
  CXX_VISIBILITY_PRESET hidden
  VISIBILITY_INLINES_HIDDEN ON
  CXX_STANDARD 17
  CXX_STANDARD_REQUIRED ON)

if (CLANG_TIDY_EXE)
  set_target_properties(libbaryonyx PROPERTIES
    CXX_CLANG_TIDY "${DO_CLANG_TIDY}")
endif ()

target_compile_definitions(libbaryonyx
  PRIVATE
  $<$<BOOL:${WITH_FULL_OPTIMIZATION}>:BARYONYX_FULL_OPTIMIZATION>
  $<$<BOOL:${WITH_DEBUG}>:BARYONYX_ENABLE_DEBUG>
  $<$<BOOL:${WITH_LOG}>:BARYONYX_ENABLE_LOG>
  $<$<BOOL:${NLOPT_FOUND}>:BARYONYX_HAVE_NLOPT>
  $<$<CXX_COMPILER_ID:MSVC>:_CRT_SECURE_NO_WARNINGS>
  $<$<CXX_COMPILER_ID:MSVC>:_SCL_SECURE_NO_WARNINGS>
  VERSION_MAJOR=${PROJECT_VERSION_MAJOR}
  VERSION_MINOR=${PROJECT_VERSION_MINOR}
  VERSION_PATCH=${PROJECT_VERSION_PATCH})

install(TARGETS libbaryonyx
  EXPORT libbaryonyx-targets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

configure_file(baryonyx.pc.in
  "${CMAKE_BINARY_DIR}/baryonyx-${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.pc"
  @ONLY)

install(FILES
  ${CMAKE_BINARY_DIR}/baryonyx-${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.pc
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)

install(DIRECTORY include/baryonyx
  DESTINATION
  ${CMAKE_INSTALL_INCLUDEDIR}/baryonyx-${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR})

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/baryonyx/export
  DESTINATION
  ${CMAKE_INSTALL_INCLUDEDIR}/baryonyx-${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}/baryonyx)

install(FILES test/sudoku.lp test/assignment_problem_1.lp
  test/assignment_problem_2.lp test/assignment_problem_3.lp
  DESTINATION
  ${CMAKE_INSTALL_DATADIR}/doc/baryonyx-${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR})

install(EXPORT libbaryonyx-targets
    FILE libbaryonyx-targets.cmake
    NAMESPACE libbaryonyx::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/baryonyx-${VLE_ABI})

# This makes the project importable from the build directory
export(TARGETS libbaryonyx FILE LibbaryonyxConfig.cmake)

function(baryonyx_add_test test_name)
  add_executable(${test_name} ${ARGN})

  set_target_properties(${test_name} PROPERTIES
    COMPILE_DEFINITIONS EXAMPLES_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/test\"
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON)

  target_compile_definitions(${test_name}
    PRIVATE
    $<$<BOOL:${WITH_DEBUG}>:BARYONYX_ENABLE_DEBUG>
    $<$<BOOL:${WITH_LOG}>:BARYONYX_ENABLE_LOG>
    $<$<BOOL:${NLOPT_FOUND}>:BARYONYX_HAVE_NLOPT>
    $<$<CXX_COMPILER_ID:MSVC>:_CRT_SECURE_NO_WARNINGS>
    $<$<CXX_COMPILER_ID:MSVC>:_SCL_SECURE_NO_WARNINGS>)

  target_include_directories(${test_name} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    PRIVATE
    src
    ${CMAKE_CURRENT_BINARY_DIR})

  target_link_libraries(${test_name} libbaryonyx threads)

  add_test(${test_name} ${test_name})
endfunction()

if (NOT BUILD_SHARED_LIBS)
  baryonyx_add_test(testio test/io.cpp)
  baryonyx_add_test(testlib test/lib.cpp)
  baryonyx_add_test(testpreprocess test/preprocess.cpp)
  baryonyx_add_test(testsolver test/solve.cpp)
  baryonyx_add_test(testoptimizer test/optimize.cpp)
endif ()
