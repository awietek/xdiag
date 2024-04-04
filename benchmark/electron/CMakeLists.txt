cmake_minimum_required(VERSION 3.15)

project(
  electron_no_sym
  VERSION 1.0
  LANGUAGES CXX
)

find_package(xdiag REQUIRED HINTS "../../../")
add_executable(main main.cpp)

target_compile_features(main PUBLIC cxx_std_17)
target_compile_definitions(main PUBLIC ${XDIAG_DEFINITIONS})
target_link_libraries(main PUBLIC ${XDIAG_LIBRARIES})
target_include_directories(main PUBLIC ${XDIAG_INCLUDE_DIRS})
set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Build type" FORCE)
