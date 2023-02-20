cmake_minimum_required(VERSION 3.15)

project(
  hb_chain_dynamical_sf_fulled_nosyms
  VERSION 1.0
  LANGUAGES CXX
)

add_executable(main main.cpp)

# set hydra compile options
find_package(hydra REQUIRED HINTS "../../..")
target_compile_features(main PUBLIC cxx_std_17)
target_compile_definitions(main PUBLIC ${HYDRA_DEFINITIONS})
target_link_libraries(main PUBLIC ${HYDRA_LIBRARIES})
target_include_directories(main PUBLIC ${HYDRA_INCLUDE_DIRS})
set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
