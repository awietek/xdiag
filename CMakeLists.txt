cmake_minimum_required(VERSION 3.15)

project(
  hydra
  VERSION 1.0
  LANGUAGES CXX
)

option(BUILD_TESTS "Build test programs" OFF)

# Set release build as default
if(NOT CMAKE_BUILD_TYPE)
 set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})

add_subdirectory(hydra)

if(BUILD_TESTS)
  message(STATUS "Building test programs")
  enable_testing()
  add_subdirectory(test)
endif()
