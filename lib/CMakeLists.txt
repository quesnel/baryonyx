project(liblp VERSION 0.1.0.0 LANGUAGES CXX C)

#configure_file(src/config.hpp.in ${CMAKE_BINARY_DIR}/config.hpp)

add_library(liblp SHARED src/lpcore.cpp src/mitm.cpp)

target_include_directories(liblp PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  PRIVATE src ${CMAKE_BINARY_DIR})

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

# target_compile_features(lib
#  PUBLIC cxx_constexpr cxx_auto_type cxx_variadic_templates)

target_link_libraries(liblp boost eigen3 threads)

set_target_properties(liblp PROPERTIES
  OUTPUT_NAME "lp-${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}")

install(TARGETS liblp
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin)

install(DIRECTORY include/ DESTINATION include/lp)

install(FILES test/sudoku.lp test/assignment_problem_1.lp
              test/assignment_problem_2.lp test/assignment_problem_3.lp
        DESTINATION share/doc/lp)

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
