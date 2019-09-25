cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project(khuri)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_definitions(-Wall -pedantic -O3)

set(SOURCE_DIR "./cpp/src")
set(INCLUDE_DIR "./cpp/include")
set(BINDING_DIR "./cpp/bindings")

include_directories(${INCLUDE_DIR})

set(SOURCES "${SOURCE_DIR}/cauchy.cpp" "${SOURCE_DIR}/gsl_interface.cpp" "${SOURCE_DIR}/omnes.cpp")
add_subdirectory(lib/pybind11)
pybind11_add_module(_khuri_omnes ${SOURCES} "${BINDING_DIR}/omnes_bindings.cpp")
target_link_libraries(_khuri_omnes PRIVATE gsl gslcblas)
