cmake_minimum_required(VERSION 3.19)

project(spinhalf_chain_e0_syms)

add_executable(main main.cpp)
find_package(xdiag REQUIRED HINTS "../../../install")
target_link_libraries(main PUBLIC xdiag::xdiag)
