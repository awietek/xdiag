cmake_minimum_required(VERSION 3.19)

project(
  tj_square_time_evolve
  VERSION 1.0
  LANGUAGES CXX
)

find_package(hydra REQUIRED HINTS "../../../install")
add_executable(main main.cpp)
target_link_libraries(main PUBLIC hydra::hydra)
