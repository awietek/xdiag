# SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
#
# SPDX-License-Identifier: Apache-2.0

message(STATUS "------   BUILDING JULIA WRAPPER XDIAG LIBRARY   ----------")

find_package(JlCxx REQUIRED)
get_target_property(JlCxx_location JlCxx::cxxwrap_julia LOCATION)
get_filename_component(JlCxx_location ${JlCxx_location} DIRECTORY)

message(STATUS ${JlCxx_location})

# Expand the kernel instantiation groups exactly like xdiag/CMakeLists.txt, so
# the grouped sources are compiled as many small per-group translation units
# rather than one monolithic TU each (which would be very slow / memory-heavy
# and inflate the shared library).
detect_grouped_sources(
  "${XDIAG_SOURCES}"
  XDIAG_GROUPED_SOURCES
  XDIAG_NORMAL_SOURCES)

set(SOURCE_ROOT  "${CMAKE_SOURCE_DIR}/xdiag")
set(INCLUDE_ROOT "${CMAKE_SOURCE_DIR}")
set(INST_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}")

set(XDIAG_INSTANTIATION_SOURCES)
foreach(SRC IN LISTS XDIAG_GROUPED_SOURCES)
  set(ABSOLUTE_SRC "${SOURCE_ROOT}/${SRC}")
  get_instantiation_outputs(
    "${ABSOLUTE_SRC}" "${SOURCE_ROOT}" "${INST_OUTPUT_DIR}" GEN_FILES)
  add_custom_command(
    OUTPUT  ${GEN_FILES}
    COMMAND ${CMAKE_COMMAND}
      -DSOURCE_FILE=${ABSOLUTE_SRC}
      -DSOURCE_ROOT=${SOURCE_ROOT}
      -DINCLUDE_ROOT=${INCLUDE_ROOT}
      -DOUTPUT_DIR=${INST_OUTPUT_DIR}
      -P ${CMAKE_SOURCE_DIR}/cmake/generate_instantiation.cmake
    DEPENDS "${ABSOLUTE_SRC}"
    COMMENT "Generating instantiation files for ${SRC}")
  list(APPEND XDIAG_INSTANTIATION_SOURCES ${GEN_FILES})
endforeach()

# Non-grouped core sources live under ../xdiag/; wrapper sources under src/.
list(TRANSFORM XDIAG_NORMAL_SOURCES PREPEND "../xdiag/")
list(TRANSFORM XDIAG_JULIA_SOURCES PREPEND "src/")

# The generated wrapper is split per subsystem (misc/wrapgen/emit_cpp.py); glob
# whatever it produced so adding/removing a subsystem file needs no edit here.
file(GLOB XDIAG_JULIA_GENERATED CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/src/generated/*.cpp")

add_library(${XDIAG_LIBRARY} SHARED
  ${XDIAG_NORMAL_SOURCES}
  ${XDIAG_INSTANTIATION_SOURCES}
  ${XDIAG_JULIA_SOURCES}
  ${XDIAG_JULIA_GENERATED})

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
