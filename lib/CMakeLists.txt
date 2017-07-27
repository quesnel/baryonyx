project(libbaryonyx VERSION 0.1.0.0 LANGUAGES CXX C)

option(WITH_FULL_OPTIMIZATION "Disable all logging facilities and active heavy optimization code. [default: off]" OFF)

option(WITH_DEBUG "enable debug log message. [default: ON]" ON)

add_library(libbaryonyx SHARED src/lpcore.cpp src/mitm.cpp src/utils.cpp)
target_include_directories(libbaryonyx PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  PRIVATE src ${CMAKE_BINARY_DIR})
target_link_libraries(libbaryonyx threads)
set_target_properties(libbaryonyx PROPERTIES
  OUTPUT_NAME "baryonyx-${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}"
  COMPILE_DEFINITIONS "BARYONYX_DLL"
  COMPILE_FLAGS "-fvisibility=hidden -fvisibility-inlines-hidden")
target_compile_definitions(libbaryonyx
  PRIVATE
  $<$<BOOL:${WITH_FULL_OPTIMIZATION}>:BARYONYX_FULL_OPTIMIZATION>
  $<$<BOOL:${WITH_DEBUG}>:BARYONYX_DISABLE_DEBUG>
  VERSION_MAJOR=${PROJECT_VERSION_MAJOR}
  VERSION_MINOR=${PROJECT_VERSION_MINOR}
  VERSION_PATCH=${PROJECT_VERSION_PATCH}
  VERSION_TWEAK=${PROJECT_VERSION_TWEAK})
install(TARGETS libbaryonyx
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin)

add_library(libbaryonyx-static STATIC src/lpcore.cpp src/mitm.cpp src/utils.cpp)
target_include_directories(libbaryonyx-static PRIVATE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  src ${CMAKE_BINARY_DIR})
target_link_libraries(libbaryonyx-static threads)
set_target_properties(libbaryonyx-static PROPERTIES
  POSITION_INDEPENDENT_CODE true
  OUTPUT_NAME "baryonyx-${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}"
  DEFINE_SYMBOL "libbaryonyx_EXPORTS"
  COMPILE_FLAGS "-fvisibility=hidden -fvisibility-inlines-hidden")
target_compile_definitions(libbaryonyx-static
  PRIVATE
  $<$<BOOL:${WITH_FULL_OPTIMIZATION}>:BARYONYX_FULL_OPTIMIZATION>
  $<$<BOOL:${WITH_DEBUG}>:BARYONYX_DISABLE_DEBUG>
  VERSION_MAJOR=${PROJECT_VERSION_MAJOR}
  VERSION_MINOR=${PROJECT_VERSION_MINOR}
  VERSION_PATCH=${PROJECT_VERSION_PATCH}
  VERSION_TWEAK=${PROJECT_VERSION_TWEAK})
install(TARGETS libbaryonyx-static
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

configure_file(baryonyx.pc.in
  "${CMAKE_BINARY_DIR}/baryonyx-${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.pc"
  @ONLY)

install(FILES
  "${CMAKE_BINARY_DIR}/baryonyx-${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.pc"
  DESTINATION "lib/pkgconfig")

# target_compile_features(lib
#  PUBLIC cxx_constexpr cxx_auto_type cxx_variadic_templates)


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
  COMPILE_DEFINITIONS EXAMPLES_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/test\")
target_include_directories(testio PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  PRIVATE src)
target_link_libraries(testio libbaryonyx threads)
add_test(testio testio)

add_executable(testlib test/lib.cpp)
set_target_properties(testlib PROPERTIES
  COMPILE_DEFINITIONS EXAMPLES_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/test\")
target_include_directories(testlib PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  PRIVATE src)
target_link_libraries(testlib libbaryonyx threads)
add_test(testlib testlib)

add_executable(testsolver test/solve.cpp)
set_target_properties(testsolver PROPERTIES
  COMPILE_DEFINITIONS EXAMPLES_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/test\")
target_include_directories(testsolver PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  PRIVATE src)
target_link_libraries(testsolver libbaryonyx threads)
add_test(testsolver testsolver)

add_executable(testoptimizer test/optimize.cpp)
set_target_properties(testoptimizer PROPERTIES
  COMPILE_DEFINITIONS EXAMPLES_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/test\")
target_include_directories(testoptimizer PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  PRIVATE src)
target_link_libraries(testoptimizer libbaryonyx threads)
add_test(testoptimizer testoptimizer)
