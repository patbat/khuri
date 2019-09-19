cmake_minimum_required(VERSION 3.5 FATAL_ERROR)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS}  ../include)

add_definitions(-Wall -pedantic -O3)

add_executable(gsl_test
    gsl_interface.cpp
    ../src/gsl_interface.cpp)
target_link_libraries(gsl_test ${GTEST_LIBRARIES} pthread gsl gslcblas)

## tests ######################################################################

enable_testing()

add_test(NAME gsl
    COMMAND ./gsl_test)
