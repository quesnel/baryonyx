project(libbaryonyx VERSION 0.2.99.0 LANGUAGES CXX)

option(WITH_FULL_OPTIMIZATION "Disable all logging facilities and active heavy optimization code. [default: off]" OFF)

option(WITH_DEBUG "enable debug log message. [default: ON]" ON)

add_library(libbaryonyx SHARED include/baryonyx/core-out include/baryonyx/core-test include/baryonyx/core include/baryonyx/core-compare src/knapsack-dp-solver.hpp src/private.hpp src/itm.cpp src/scoped_array.hpp src/select.cpp src/preprocessor.cpp src/fixed_array.hpp src/consistency.cpp src/lpcore.cpp src/branch-and-bound-solver.hpp src/utils.hpp src/fixed_2darray.hpp src/matrix.hpp src/itm-solver-inequalities-Z.hpp src/itm-solver-common.hpp src/itm-common.hpp src/sol-format.cpp src/itm-solver-equalities-01.hpp src/lpformat-io.cpp src/itm.hpp)

target_include_directories(libbaryonyx PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  PRIVATE src ${CMAKE_BINARY_DIR})
target_link_libraries(libbaryonyx threads)
set_target_properties(libbaryonyx PROPERTIES
  OUTPUT_NAME "baryonyx-${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}"
  POSITION_INDEPENDENT_CODE ON
  COMPILE_DEFINITIONS "BARYONYX_DLL"
  DEFINE_SYMBOL "libbaryonyx_EXPORTS"
  CXX_VISIBILITY_PRESET hidden
  VISIBILITY_INLINES_HIDDEN ON
  CXX_STANDARD 14
  CXX_STANDARD_REQUIRED ON)
target_compile_definitions(libbaryonyx
  PRIVATE
  $<$<BOOL:${WITH_FULL_OPTIMIZATION}>:BARYONYX_FULL_OPTIMIZATION>
  $<$<BOOL:${WITH_DEBUG}>:BARYONYX_DISABLE_DEBUG>
  $<$<CXX_COMPILER_ID:MSVC>:_CRT_SECURE_NO_WARNINGS>
  $<$<CXX_COMPILER_ID:MSVC>:_SCL_SECURE_NO_WARNINGS>
  VERSION_MAJOR=${PROJECT_VERSION_MAJOR}
  VERSION_MINOR=${PROJECT_VERSION_MINOR}
  VERSION_PATCH=${PROJECT_VERSION_PATCH}
  VERSION_TWEAK=${PROJECT_VERSION_TWEAK})
install(TARGETS libbaryonyx
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin)

add_library(libbaryonyx-static STATIC include/baryonyx/core-out include/baryonyx/core-test include/baryonyx/core include/baryonyx/core-compare src/knapsack-dp-solver.hpp src/private.hpp src/itm.cpp src/scoped_array.hpp src/select.cpp src/preprocessor.cpp src/fixed_array.hpp src/consistency.cpp src/lpcore.cpp src/branch-and-bound-solver.hpp src/utils.hpp src/fixed_2darray.hpp src/matrix.hpp src/itm-solver-inequalities-Z.hpp src/itm-solver-common.hpp src/itm-common.hpp src/sol-format.cpp src/itm-solver-equalities-01.hpp src/lpformat-io.cpp src/itm.hpp)

target_include_directories(libbaryonyx-static PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  src ${CMAKE_BINARY_DIR})
target_link_libraries(libbaryonyx-static threads)
set_target_properties(libbaryonyx-static PROPERTIES
  POSITION_INDEPENDENT_CODE true
  OUTPUT_NAME "baryonyx-${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}"
  CXX_STANDARD 14
  CXX_STANDARD_REQUIRED ON)
target_compile_definitions(libbaryonyx-static
  PRIVATE
  $<$<BOOL:${WITH_FULL_OPTIMIZATION}>:BARYONYX_FULL_OPTIMIZATION>
  $<$<BOOL:${WITH_DEBUG}>:BARYONYX_DISABLE_DEBUG>
  $<$<CXX_COMPILER_ID:MSVC>:_CRT_SECURE_NO_WARNINGS>
  $<$<CXX_COMPILER_ID:MSVC>:_SCL_SECURE_NO_WARNINGS>
  VERSION_MAJOR=${PROJECT_VERSION_MAJOR}
  VERSION_MINOR=${PROJECT_VERSION_MINOR}
  VERSION_PATCH=${PROJECT_VERSION_PATCH}
  VERSION_TWEAK=${PROJECT_VERSION_TWEAK})
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
target_link_libraries(testio libbaryonyx threads)
add_test(testio testio)

add_executable(testlib test/lib.cpp)
set_target_properties(testlib PROPERTIES
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
target_link_libraries(testlib libbaryonyx threads)
add_test(testlib testlib)

add_executable(testsolver test/solve.cpp)
set_target_properties(testsolver PROPERTIES
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
target_link_libraries(testsolver libbaryonyx threads)
add_test(testsolver testsolver)

add_executable(testoptimizer test/optimize.cpp)
set_target_properties(testoptimizer PROPERTIES
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
target_link_libraries(testoptimizer libbaryonyx threads)
add_test(testoptimizer testoptimizer)
