cmake_minimum_required(VERSION 3.19)

project(hubbard_greens_f)

find_package(xdiag REQUIRED HINTS "../../install")
add_executable(main main.cpp)
target_link_libraries(main PUBLIC xdiag::xdiag)
