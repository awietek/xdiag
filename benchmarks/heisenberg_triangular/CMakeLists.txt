cmake_minimum_required(VERSION 3.19)

project(Ybmat)

find_package(xdiag REQUIRED HINTS "$ENV{HOME}/Research/Software/xdiag/install/")
add_executable(main main.cpp)
target_link_libraries(main PRIVATE xdiag::xdiag)
target_compile_options(main PRIVATE -O3 -march=native)
