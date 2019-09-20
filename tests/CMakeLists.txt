cmake_minimum_required(VERSION 3.5 FATAL_ERROR)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# BEGIN GOOGLETEST SETUP ######################################################

# Download and unpack googletest at configure time
configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif()

# Prevent overriding the parent project's compiler/linker
# settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This defines
# the gtest and gtest_main targets.
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
                 ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
                 EXCLUDE_FROM_ALL)

# END GOOGLETEST SETUP ########################################################

include_directories(../include)

add_definitions(-Wall -pedantic -O3)

add_executable(test_gsl
    test_gsl_interface.cpp
    ../src/gsl_interface.cpp)
target_link_libraries(test_gsl gtest pthread gsl gslcblas)

add_executable(test_cauchy
    test_cauchy.cpp
    ../src/gsl_interface.cpp
    ../src/cauchy.cpp)
target_link_libraries(test_cauchy gmock gtest pthread gsl gslcblas)

enable_testing()

add_test(NAME gsl
    COMMAND ./test_gsl)

add_test(NAME cauchy
    COMMAND ./test_cauchy)
