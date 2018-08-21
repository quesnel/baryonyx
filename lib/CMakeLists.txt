project(libbaryonyx VERSION 0.2.99.0 LANGUAGES CXX)

option(WITH_FULL_OPTIMIZATION "Disable all logging facilities and active heavy optimization code. [default: off]" OFF)

option(WITH_DEBUG "enable debug log message. [default: ON]" ON)

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
 src/fixed-2darray.hpp
 src/fixed-array.hpp
 src/itm-common.hpp
 src/itm-optimizer-common.hpp
 src/itm-solver-common.hpp
 src/itm-solver-equalities-01.cpp
 src/itm-solver-equalities-01.hpp
 src/itm-solver-equalities-101.cpp
 src/itm-solver-equalities-101.hpp
 src/itm-solver-inequalities-01.cpp
 src/itm-solver-inequalities-01.hpp
 src/itm-solver-inequalities-101-buffered.cpp
 src/itm-solver-inequalities-101-buffered.hpp
 src/itm-solver-inequalities-101.cpp
 src/itm-solver-inequalities-101.hpp
 src/itm-solver-inequalities-Z.cpp
 src/itm-solver-inequalities-Z.hpp
 src/itm.hpp
 src/knapsack-dp-solver.hpp
 src/lpcore.cpp
 src/memory.cpp
 src/memory.hpp
 src/merged-constraint.cpp
 src/observer.hpp
 src/pnm.hpp
 src/preprocessor-2.cpp
 src/preprocessor.cpp
 src/private.cpp
 src/private.hpp
 src/problem.cpp
 src/problem.hpp
 src/result.cpp
 src/result.hpp
 src/resume.hpp
 src/sparse-matrix.hpp
 src/sparse-vector.hpp
 src/utils.hpp)

add_library(libbaryonyx SHARED ${public_baryonyx_header} ${private_baryonyx_source})

target_include_directories(libbaryonyx PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  PRIVATE src ${CMAKE_BINARY_DIR})

target_link_libraries(libbaryonyx threads
  $<$<BOOL:${NLOPT_FOUND}>:nlopt>)

set_target_properties(libbaryonyx PROPERTIES
  OUTPUT_NAME "baryonyx-${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}"
  POSITION_INDEPENDENT_CODE ON
  COMPILE_DEFINITIONS "BARYONYX_DLL"
  DEFINE_SYMBOL "libbaryonyx_EXPORTS"
  INTERPROCEDURAL_OPTIMIZATION ON
  CXX_VISIBILITY_PRESET hidden
  VISIBILITY_INLINES_HIDDEN ON
  CXX_STANDARD 14
  CXX_STANDARD_REQUIRED ON)

target_compile_definitions(libbaryonyx
  PRIVATE
  $<$<BOOL:${WITH_FULL_OPTIMIZATION}>:BARYONYX_FULL_OPTIMIZATION>
  $<$<BOOL:${WITH_DEBUG}>:BARYONYX_DISABLE_DEBUG>
  $<$<BOOL:${NLOPT_FOUND}>:BARYONYX_HAVE_NLOPT>
  $<$<CXX_COMPILER_ID:MSVC>:_CRT_SECURE_NO_WARNINGS>
  $<$<CXX_COMPILER_ID:MSVC>:_SCL_SECURE_NO_WARNINGS>
  VERSION_MAJOR=${PROJECT_VERSION_MAJOR}
  VERSION_MINOR=${PROJECT_VERSION_MINOR}
  VERSION_PATCH=${PROJECT_VERSION_PATCH}
  VERSION_TWEAK=${PROJECT_VERSION_TWEAK})

target_compile_options(libbaryonyx PRIVATE
  $<$<CXX_COMPILER_ID:MSVC>:/bigobj>)

install(TARGETS libbaryonyx
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin)

add_library(libbaryonyx-static STATIC ${public_baryonyx_header} ${private_baryonyx_source})

target_include_directories(libbaryonyx-static PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  src ${CMAKE_BINARY_DIR})

target_link_libraries(libbaryonyx-static threads
  $<$<BOOL:${NLOPT_FOUND}>:nlopt>)

set_target_properties(libbaryonyx-static PROPERTIES
  POSITION_INDEPENDENT_CODE true
  OUTPUT_NAME "baryonyx-static-${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}"
  INTERPROCEDURAL_OPTIMIZATION ON
  CXX_STANDARD 14
  CXX_STANDARD_REQUIRED ON)

target_compile_definitions(libbaryonyx-static
  PRIVATE
  $<$<BOOL:${WITH_FULL_OPTIMIZATION}>:BARYONYX_FULL_OPTIMIZATION>
  $<$<BOOL:${WITH_DEBUG}>:BARYONYX_DISABLE_DEBUG>
  $<$<BOOL:${NLOPT_FOUND}>:BARYONYX_HAVE_NLOPT>
  $<$<CXX_COMPILER_ID:MSVC>:_CRT_SECURE_NO_WARNINGS>
  $<$<CXX_COMPILER_ID:MSVC>:_SCL_SECURE_NO_WARNINGS>
  VERSION_MAJOR=${PROJECT_VERSION_MAJOR}
  VERSION_MINOR=${PROJECT_VERSION_MINOR}
  VERSION_PATCH=${PROJECT_VERSION_PATCH}
  VERSION_TWEAK=${PROJECT_VERSION_TWEAK})

target_compile_options(libbaryonyx-static PRIVATE
  $<$<CXX_COMPILER_ID:MSVC>:/bigobj>)

install(TARGETS libbaryonyx-static
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin)

configure_file(baryonyx.pc.in
  "${CMAKE_BINARY_DIR}/baryonyx-${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.pc"
  @ONLY)

install(FILES
  "${CMAKE_BINARY_DIR}/baryonyx-${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.pc"
  DESTINATION "lib/pkgconfig")

install(DIRECTORY include/ DESTINATION
  include/baryonyx-${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR})

install(FILES test/sudoku.lp test/assignment_problem_1.lp
  test/assignment_problem_2.lp test/assignment_problem_3.lp
  DESTINATION
  share/doc/baryonyx-${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR})

# This makes the project importable from the build directory
export(TARGETS libbaryonyx FILE LibbaryonyxConfig.cmake)

add_executable(testio test/io.cpp)
set_target_properties(testio PROPERTIES
  COMPILE_DEFINITIONS EXAMPLES_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/test\"
  INTERPROCEDURAL_OPTIMIZATION ON
  CXX_STANDARD 14
  CXX_STANDARD_REQUIRED ON)
target_compile_definitions(testio
  PRIVATE
  $<$<CXX_COMPILER_ID:MSVC>:_CRT_SECURE_NO_WARNINGS>
  $<$<CXX_COMPILER_ID:MSVC>:_SCL_SECURE_NO_WARNINGS>)
target_include_directories(testio PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  PRIVATE src)
target_link_libraries(testio libbaryonyx-static threads)
add_test(testio testio)

add_executable(testlib test/lib.cpp)
set_target_properties(testlib PROPERTIES
  INTERPROCEDURAL_OPTIMIZATION ON
  COMPILE_DEFINITIONS EXAMPLES_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/test\"
  CXX_STANDARD 14
  CXX_STANDARD_REQUIRED ON)
target_compile_definitions(testlib
  PRIVATE
  $<$<CXX_COMPILER_ID:MSVC>:_CRT_SECURE_NO_WARNINGS>
  $<$<CXX_COMPILER_ID:MSVC>:_SCL_SECURE_NO_WARNINGS>)
target_include_directories(testlib PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  PRIVATE src)
target_link_libraries(testlib libbaryonyx-static threads)
add_test(testlib testlib)

add_executable(testpreprocess test/preprocess.cpp)
set_target_properties(testpreprocess PROPERTIES
  INTERPROCEDURAL_OPTIMIZATION ON
  COMPILE_DEFINITIONS EXAMPLES_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/test\"
  CXX_STANDARD 14
  CXX_STANDARD_REQUIRED ON)
target_compile_definitions(testpreprocess
  PRIVATE
  $<$<CXX_COMPILER_ID:MSVC>:_CRT_SECURE_NO_WARNINGS>
  $<$<CXX_COMPILER_ID:MSVC>:_SCL_SECURE_NO_WARNINGS>)
target_include_directories(testpreprocess PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  PRIVATE src)
target_link_libraries(testpreprocess libbaryonyx-static threads)
add_test(testpreprocess testpreprocess)

add_executable(testsolver test/solve.cpp)
set_target_properties(testsolver PROPERTIES
  INTERPROCEDURAL_OPTIMIZATION ON
  COMPILE_DEFINITIONS EXAMPLES_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/test\"
  CXX_STANDARD 14
  CXX_STANDARD_REQUIRED ON)
target_compile_definitions(testsolver
  PRIVATE
  $<$<CXX_COMPILER_ID:MSVC>:_CRT_SECURE_NO_WARNINGS>
  $<$<CXX_COMPILER_ID:MSVC>:_SCL_SECURE_NO_WARNINGS>)
target_include_directories(testsolver PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  PRIVATE src)
target_link_libraries(testsolver libbaryonyx-static threads)
add_test(testsolver testsolver)

add_executable(testoptimizer test/optimize.cpp)
set_target_properties(testoptimizer PROPERTIES
  INTERPROCEDURAL_OPTIMIZATION ON
  COMPILE_DEFINITIONS EXAMPLES_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/test\"
  CXX_STANDARD 14
  CXX_STANDARD_REQUIRED ON)
target_compile_definitions(testoptimizer
  PRIVATE
  $<$<CXX_COMPILER_ID:MSVC>:_CRT_SECURE_NO_WARNINGS>
  $<$<CXX_COMPILER_ID:MSVC>:_SCL_SECURE_NO_WARNINGS>)
target_include_directories(testoptimizer PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  PRIVATE src)
target_link_libraries(testoptimizer libbaryonyx-static threads)
add_test(testoptimizer testoptimizer)
