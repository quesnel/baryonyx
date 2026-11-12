find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(_baryonyx src/module.cpp)
target_link_libraries(_baryonyx PRIVATE libbaryonyx threads)

install(TARGETS _baryonyx
  LIBRARY DESTINATION baryonyx
  COMPONENT python)
