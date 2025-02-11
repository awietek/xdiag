cmake_minimum_required(VERSION 3.19)

project(
  benchmark_heisenberg_chain
)

find_package(xdiag REQUIRED HINTS "~/Research/Software/xdiag/install")
add_executable(bench bench.cpp)
target_link_libraries(bench PRIVATE xdiag::xdiag)
