cmake_minimum_required(VERSION 3.19)

project(
  tj_distributed_time_evolve
  LANGUAGES CXX
)

find_package(hydra_distributed REQUIRED HINTS ../../../install)
add_executable(main main.cpp)
target_link_libraries(main PUBLIC hydra::hydra_distributed)
