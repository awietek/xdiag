cmake_minimum_required(VERSION 3.15)

project(
  floquet_time_evo
  VERSION 1.0
  LANGUAGES CXX
)

find_package(xdiag REQUIRED HINTS "../../../")
add_executable(main main.cpp)
