find_package(JlCxx REQUIRED)
get_target_property(JlCxx_location JlCxx::cxxwrap_julia LOCATION)
get_filename_component(JlCxx_location ${JlCxx_location} DIRECTORY)

add_library(hydrajl SHARED hydrajl.cpp)

target_link_libraries(hydrajl PRIVATE JlCxx::cxxwrap_julia)
target_link_libraries(hydrajl PRIVATE ${JLCXX_TARGET} ${JLCXX_STL_TARGET} ${Julia_LIBRARY})
target_link_libraries(hydrajl PRIVATE hydra)
