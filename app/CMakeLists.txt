add_executable(app src/main.cpp)

target_link_libraries(app libbaryonyx-static threads)

set_target_properties(app PROPERTIES
  OUTPUT_NAME "baryonyx-${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}"
  CXX_STANDARD 14
  CXX_STANDARD_REQUIRED ON)

install(TARGETS app
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin)
