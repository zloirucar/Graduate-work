cmake_minimum_required (VERSION 2.8.3) 
project (vibro-transport-package) 

set (CMAKE_CXX_STANDARD 11)

set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set (BUILD_SHARED_LIBS true)

include_directories(${PROJECT_SOURCE_DIR}/external/ode_num_int/include)

add_subdirectory(external)
add_subdirectory(vibro-transport)

