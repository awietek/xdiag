cmake_minimum_required(VERSION 3.19)

project(
  checkerboard_tpq
)

add_executable(main main.cpp)
find_package(xdiag REQUIRED HINTS "../../..")
target_link_libraries(main PUBLIC xdiag::xdiag)
