cmake_minimum_required(VERSION 3.15)

project(
  hydra
  VERSION 1.0
  LANGUAGES CXX
)

option(BUILD_TESTS "Build test programs" OFF)
option(BUILD_JULIA "Build the Julia wrapper" OFF)

option(OPTIMIZE_FOR_NATIVE "Optimize for native architecture" OFF)
option(OMP_THREADING "Enable OpenMp threading" ON)
option(MKL_THREADING "Enable IntelMKL threading" ON)
option(SANITIZER "Build with -fsanitize=address" OFF)

message(STATUS "Compiler ID     : " ${CMAKE_CXX_COMPILER_ID})
message(STATUS "Compiler version: " ${CMAKE_CXX_COMPILER_VERSION})

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

if(BUILD_JULIA)
  message(STATUS "Building Julia wrapper")
  set_property(TARGET hydra PROPERTY POSITION_INDEPENDENT_CODE ON)
  add_subdirectory(julia)
endif()

if(SANITIZER)
    message(STATUS "Compiling with sanitizers")
    target_compile_options(hydra PUBLIC -fsanitize=address)
    target_link_options(hydra PUBLIC -fsanitize=address)
endif()
