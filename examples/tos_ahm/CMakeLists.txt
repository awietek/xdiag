cmake_minimum_required(VERSION 3.19)

project(tos_ahm)

find_package(xdiag REQUIRED HINTS "../xdiag/install")
add_executable(main tos_ahm.cpp)
target_link_libraries(main PRIVATE xdiag::xdiag)
target_compile_options(main PRIVATE -Wall -march=native)
