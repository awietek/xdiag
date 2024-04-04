cmake_minimum_required(VERSION 3.19)

project(
  tj_distributed_time_evolve
)

find_package(xdiag_distributed REQUIRED HINTS ../../../install)
add_executable(main main.cpp)
target_link_libraries(main PUBLIC xdiag::xdiag_distributed)
