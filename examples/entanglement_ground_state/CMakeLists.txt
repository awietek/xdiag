cmake_minimum_required(VERSION 3.19)

project(entanglement_ground_state)

find_package(xdiag REQUIRED HINTS "../../install")
add_executable(main main.cpp)
target_link_libraries(main PUBLIC xdiag::xdiag)