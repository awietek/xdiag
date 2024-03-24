cmake_minimum_required(VERSION 3.15)

project(
  tj_distributed_time_evolve
  VERSION 1.0
  LANGUAGES CXX
)

find_package(hydra_distributed CONFIG REQUIRED HINTS "/home/awietek/Research/Software/hydra")
add_executable(main main.cpp)
target_link_libraries(main hydra_distributed)
