cmake_minimum_required(VERSION 3.19)

project(
  hb_chain_dynamical_sf_nosyms
)

add_executable(main main.cpp)
find_package(xdiag REQUIRED HINTS "../../..")
target_link_libraries(xdiag PUBLIC ${XDIAG_LIBRARIES})
