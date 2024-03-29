cmake_minimum_required(VERSION 3.19)

project(
  hello_world
  LANGUAGES CXX
)

find_package(hydra REQUIRED HINTS "../../install")
add_executable(main main.cpp)
target_link_libraries(main PUBLIC hydra::hydra)
