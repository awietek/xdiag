# SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
#
# SPDX-License-Identifier: Apache-2.0

message(STATUS "------   BUILDING JULIA WRAPPER XDIAG LIBRARY   ----------")

find_package(JlCxx REQUIRED)
get_target_property(JlCxx_location JlCxx::cxxwrap_julia LOCATION)
get_filename_component(JlCxx_location ${JlCxx_location} DIRECTORY)
message(STATUS ${JlCxx_location})

# Julia wrapper source list (relative to this julia/ directory).
include(${CMAKE_CURRENT_SOURCE_DIR}/src/xdiag_julia_sources.cmake)

# ---------------------------------------------------------------------------
# Instantiation-group expansion (mirrors xdiag/CMakeLists.txt)
#
# The grouped base sources carry BEGIN_INSTANTIATION_GROUP markers and must be
# split into one .cpp per group at build time, exactly like the main library --
# otherwise each grouped file would be compiled as a single, huge translation
# unit.  detect_grouped_sources() has already run at the top level, populating
# XDIAG_GROUPED_SOURCES and XDIAG_NORMAL_SOURCES (both relative to xdiag/), and
# those variables are inherited into this subdirectory scope.
# ---------------------------------------------------------------------------
set(SOURCE_ROOT     "${CMAKE_SOURCE_DIR}/xdiag")
set(INCLUDE_ROOT    "${CMAKE_SOURCE_DIR}")
set(INST_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}")

# For each grouped source: compute the generated file paths at configure time,
# then register a build-time command to (re)generate one .cpp per group
# whenever the original source changes.
set(XDIAG_INSTANTIATION_SOURCES)
foreach(SRC IN LISTS XDIAG_GROUPED_SOURCES)
  set(ABSOLUTE_SRC "${SOURCE_ROOT}/${SRC}")

  get_instantiation_outputs(
    "${ABSOLUTE_SRC}" "${SOURCE_ROOT}" "${INST_OUTPUT_DIR}"
    GEN_FILES)

  add_custom_command(
    OUTPUT  ${GEN_FILES}
    COMMAND ${CMAKE_COMMAND}
      -DSOURCE_FILE=${ABSOLUTE_SRC}
      -DSOURCE_ROOT=${SOURCE_ROOT}
      -DINCLUDE_ROOT=${INCLUDE_ROOT}
      -DOUTPUT_DIR=${INST_OUTPUT_DIR}
      -P ${CMAKE_SOURCE_DIR}/cmake/generate_instantiation.cmake
    DEPENDS "${ABSOLUTE_SRC}"
    COMMENT "Generating instantiation files for ${SRC}"
  )

  list(APPEND XDIAG_INSTANTIATION_SOURCES ${GEN_FILES})
endforeach()

# The non-grouped base sources are relative to xdiag/; make them absolute so
# they resolve from this julia/ subdirectory.
set(XDIAG_NORMAL_SOURCES_ABS)
foreach(SRC IN LISTS XDIAG_NORMAL_SOURCES)
  list(APPEND XDIAG_NORMAL_SOURCES_ABS "${SOURCE_ROOT}/${SRC}")
endforeach()

add_library(${XDIAG_LIBRARY} SHARED
  ${XDIAG_NORMAL_SOURCES_ABS}
  ${XDIAG_INSTANTIATION_SOURCES}
  ${XDIAG_JULIA_SOURCES})

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
