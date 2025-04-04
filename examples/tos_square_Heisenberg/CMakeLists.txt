cmake_minimum_required(VERSION 3.19)

project(
  square-Heisenberg-TOS
)

add_executable(main main.cpp)
find_package(xdiag REQUIRED HINTS "/home/ssarkar/research/xdiag/install")
target_link_libraries(main PRIVATE xdiag::xdiag)
