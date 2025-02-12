cmake_minimum_required(VERSION 3.19)

project(
  tj_random_all_to_all
)

add_executable(main main.cpp)
find_package(xdiag REQUIRED HINTS "../../../install")
target_link_libraries(main PRIVATE xdiag::xdiag)
