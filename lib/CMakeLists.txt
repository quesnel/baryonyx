project(liblp VERSION 0.1.0.0 LANGUAGES CXX C)

option(WITH_FULL_OPTIMIZATION "Disable all logging facilities and active heavy optimization code. [default: off]" OFF)

option(WITH_DEBUG "enable debug log message. [default: ON]" ON)

add_library(liblp SHARED src/lpcore.cpp src/mitm.cpp src/utils.cpp)
target_include_directories(liblp PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  PRIVATE src ${CMAKE_BINARY_DIR})
target_link_libraries(liblp boost eigen3 threads)
set_target_properties(liblp PROPERTIES
  OUTPUT_NAME "lp-${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}"
  COMPILE_DEFINITIONS "LP_DLL"
  COMPILE_FLAGS "-fvisibility=hidden -fvisibility-inlines-hidden")
target_compile_definitions(liblp
  PRIVATE
  $<$<BOOL:${WITH_FULL_OPTIMIZATION}>:LP_FULL_OPTIMIZATION>
  $<$<BOOL:${WITH_DEBUG}>:LP_DISABLE_DEBUG>
  VERSION_MAJOR=${PROJECT_VERSION_MAJOR}
  VERSION_MINOR=${PROJECT_VERSION_MINOR}
  VERSION_PATCH=${PROJECT_VERSION_PATCH}
  VERSION_TWEAK=${PROJECT_VERSION_TWEAK})
install(TARGETS liblp
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin)

add_library(liblp-static STATIC src/lpcore.cpp src/mitm.cpp src/utils.cpp)
target_include_directories(liblp-static PRIVATE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  src ${CMAKE_BINARY_DIR})
target_link_libraries(liblp-static boost eigen3 threads)
set_target_properties(liblp-static PROPERTIES
  POSITION_INDEPENDENT_CODE true
  OUTPUT_NAME "lp-${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}"
  DEFINE_SYMBOL "liblp_EXPORTS"
  COMPILE_FLAGS "-fvisibility=hidden -fvisibility-inlines-hidden")
target_compile_definitions(liblp-static
  PRIVATE
  $<$<BOOL:${WITH_FULL_OPTIMIZATION}>:LP_FULL_OPTIMIZATION>
  $<$<BOOL:${WITH_DEBUG}>:LP_DISABLE_DEBUG>
  VERSION_MAJOR=${PROJECT_VERSION_MAJOR}
  VERSION_MINOR=${PROJECT_VERSION_MINOR}
  VERSION_PATCH=${PROJECT_VERSION_PATCH}
  VERSION_TWEAK=${PROJECT_VERSION_TWEAK})
install(TARGETS liblp-static
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin)


# If we have compiler requirements for this library, list them
# here
# target_compile_features(lib
#   PUBLIC cxx_auto_type
#   PRIVATE cxx_variadic_templates)

include(CheckCXXCompilerFlag)

check_cxx_compiler_flag("-std=c++14" COMPILER_SUPPORTS_CXX14)
check_cxx_compiler_flag("-std=c++1y" COMPILER_SUPPORTS_CXX1Y)

if (COMPILER_SUPPORTS_CXX14)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
elseif(COMPILER_SUPPORTS_CXX1Y)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1y")
else()
  message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++14/1y "
    "spport. Please use a different C++ compiler.")
endif()

configure_file(lp.pc.in
  "${CMAKE_BINARY_DIR}/lp-${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.pc"
  @ONLY)

install(FILES
  "${CMAKE_BINARY_DIR}/lp-${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.pc"
  DESTINATION "lib/pkgconfig")

# target_compile_features(lib
#  PUBLIC cxx_constexpr cxx_auto_type cxx_variadic_templates)


install(DIRECTORY include/ DESTINATION
  include/lp-${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR})

install(FILES test/sudoku.lp test/assignment_problem_1.lp
  test/assignment_problem_2.lp test/assignment_problem_3.lp
  DESTINATION
  share/doc/lp-${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR})

# This makes the project importable from the build directory
export(TARGETS liblp FILE LiblpConfig.cmake)

add_executable(testio test/io.cpp)
set_target_properties(testio PROPERTIES
  COMPILE_DEFINITIONS EXAMPLES_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/test\")
target_include_directories(testio PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  PRIVATE src)
target_link_libraries(testio liblp boost eigen3 threads)
add_test(testio testio)

add_executable(testlib test/lib.cpp)
set_target_properties(testlib PROPERTIES
  COMPILE_DEFINITIONS EXAMPLES_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/test\")
target_include_directories(testlib PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  PRIVATE src)
target_link_libraries(testlib liblp boost eigen3 threads)
add_test(testlib testlib)

add_executable(testsolver test/solve.cpp)
set_target_properties(testsolver PROPERTIES
  COMPILE_DEFINITIONS EXAMPLES_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/test\")
target_include_directories(testsolver PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  PRIVATE src)
target_link_libraries(testsolver liblp boost eigen3 threads)
add_test(testsolver testsolver)

add_executable(testoptimizer test/optimize.cpp)
set_target_properties(testoptimizer PROPERTIES
  COMPILE_DEFINITIONS EXAMPLES_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/test\")
target_include_directories(testoptimizer PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  PRIVATE src)
target_link_libraries(testoptimizer liblp boost eigen3 threads)
add_test(testoptimizer testoptimizer)
