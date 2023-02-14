cmake_minimum_required(VERSION 3.1...3.25)

project(
  hydra
  VERSION 1.0
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-O3 -Wall -Wextra -Wpedantic")

include(cmake/sources.cmake)
add_library(hydra STATIC ${HYDRA_SOURCES})
target_include_directories(hydra PUBLIC .)

add_subdirectory(test)

include(cmake/packages.cmake)


message( STATUS "-------- Creating config.h -------------")
find_package(Git REQUIRED)

# Get the Git Hash
execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
 WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
 OUTPUT_VARIABLE   HYDRA_GIT_HASH
 OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Compute HYDRA_HOSTNAME and HYDRA_COMPILEDBY
execute_process(COMMAND hostname 
    OUTPUT_VARIABLE HYDRA_HOSTNAME 
    RESULT_VARIABLE returncode 
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
execute_process(COMMAND id -un 
    OUTPUT_VARIABLE HYDRA_COMPILEDBY
    RESULT_VARIABLE returncode
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

configure_file(./cmake/config.h.in ${PROJECT_SOURCE_DIR}/hydra/config.h @ONLY)

