# # get julia cxxwrap path
# execute_process(COMMAND julia -e "using CxxWrap; ; println(CxxWrap.prefix_path())"
#  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
#  OUTPUT_VARIABLE JULIA_CXXWRAP_DIR
#  OUTPUT_STRIP_TRAILING_WHITESPACE
# )

# set(JlCxx_DIR "${JULIA_CXXWRAP_DIR}/lib/cmake/JlCxx")
# message(STATUS "JlCxx_DIR: " ${JlCxx_DIR})

find_package(JlCxx REQUIRED)
get_target_property(JlCxx_location JlCxx::cxxwrap_julia LOCATION)
get_filename_component(JlCxx_location ${JlCxx_location} DIRECTORY)

message(STATUS ${JlCxx_location})
add_library(hydrajl SHARED hydrajl.cpp)

target_link_libraries(hydrajl PRIVATE JlCxx::cxxwrap_julia)
target_link_libraries(hydrajl PRIVATE JlCxx::cxxwrap_julia_stl)
target_link_libraries(hydrajl PRIVATE ${JLCXX_TARGET} ${JLCXX_STL_TARGET} ${Julia_LIBRARY})
target_link_libraries(hydrajl PRIVATE ${HYDRA_LIBRARY})
