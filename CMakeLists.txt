cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project(khuri)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_definitions(-Wall -pedantic -O3)


set(SOURCE_DIR "khuri/cpp/")
set(SOURCES "${SOURCE_DIR}/cauchy.cpp" "${SOURCE_DIR}/gsl_interface.cpp" "${SOURCE_DIR}/omnes.cpp")

add_subdirectory(lib/pybind11)
pybind11_add_module(omnes ${SOURCES} "${SOURCE_DIR}/omnes_bindings.cpp")
target_link_libraries(omnes PRIVATE gsl cln gslcblas)
