# SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
#
# SPDX-License-Identifier: Apache-2.0

message(STATUS "------   BUILDING JULIA WRAPPER XDIAG LIBRARY   ----------")

find_package(JlCxx REQUIRED)
get_target_property(JlCxx_location JlCxx::cxxwrap_julia LOCATION)
get_filename_component(JlCxx_location ${JlCxx_location} DIRECTORY)

message(STATUS ${JlCxx_location})
list(TRANSFORM XDIAG_SOURCES PREPEND "../xdiag/")

message(STATUS ${JlCxx_location})
list(TRANSFORM XDIAG_JULIA_SOURCES PREPEND "src/")

add_library(${XDIAG_LIBRARY} SHARED ${XDIAG_SOURCES} ${XDIAG_JULIA_SOURCES})

target_include_directories(
  ${XDIAG_LIBRARY} PUBLIC "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>")
target_compile_features(${XDIAG_LIBRARY} PUBLIC cxx_std_17)

# set_target_properties(${XDIAG_LIBRARY} PROPERTIES
#   VERSION ${xdiag_VERSION}
#   SOVERSION ${xdiag_VERSION_MAJOR})

target_include_directories(${XDIAG_LIBRARY} PRIVATE ${JlCxx_DIR}/../../../include ${JlCxx_DIR}/../../../include/julia)
target_link_libraries(${XDIAG_LIBRARY} PRIVATE JlCxx::cxxwrap_julia)
target_link_libraries(${XDIAG_LIBRARY} PRIVATE JlCxx::cxxwrap_julia_stl)
target_link_libraries(${XDIAG_LIBRARY} PRIVATE ${JLCXX_TARGET} ${JLCXX_STL_TARGET} ${Julia_LIBRARY})
