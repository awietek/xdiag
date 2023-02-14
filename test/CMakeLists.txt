include(../cmake/test_sources.cmake)
include(../cmake/packages.cmake)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/../cmake")

add_executable(tests ${HYDRA_TEST_SOURCES})
target_link_libraries(tests PUBLIC hydra ${lapack_libraries} ${hdf5_libraries})
