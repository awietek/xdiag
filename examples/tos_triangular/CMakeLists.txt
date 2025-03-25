cmake_minimum_required(VERSION 3.19)

project(tos_triangular)

find_package(xdiag REQUIRED HINTS "../../install")
add_executable(main main.cpp)
target_link_libraries(main PUBLIC xdiag::xdiag)
