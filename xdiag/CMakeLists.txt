# SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
#
# SPDX-License-Identifier: Apache-2.0

add_library(${XDIAG_LIBRARY})
add_library(xdiag::${XDIAG_LIBRARY} ALIAS ${XDIAG_LIBRARY})
set_target_properties(${XDIAG_LIBRARY} PROPERTIES
  VERSION ${xdiag_VERSION}
  SOVERSION ${xdiag_VERSION_MAJOR})
target_include_directories(
    ${XDIAG_LIBRARY} PUBLIC "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>")
target_compile_features(${XDIAG_LIBRARY} PUBLIC cxx_std_17)

if (CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
  message(STATUS "switching to precise floating point arithmetic of IntelLLVM")
  target_compile_options(${XDIAG_LIBRARY} PRIVATE -fp-model precise)
endif()



if(XDIAG_DISTRIBUTED)
  message(STATUS "--------- BUILDING DISTRIBUTED XDIAG LIBRARY ------------")
  target_sources(${XDIAG_LIBRARY} PRIVATE ${XDIAG_SOURCES} ${XDIAG_DISTRIBUTED_SOURCES})

  find_package(MPI REQUIRED)
  target_link_libraries(xdiag_distributed PUBLIC ${MPI_CXX_LIBRARIES})
  target_compile_definitions(xdiag_distributed PUBLIC XDIAG_DISTRIBUTED)
else()
  message(STATUS "--------    BUILDING NORMAL XDIAG LIBRARY   -------------")

  # XDIAG_GROUPED_SOURCES and XDIAG_NORMAL_SOURCES are set by the top-level
  # CMakeLists.txt via detect_grouped_sources() before add_subdirectory(xdiag).

  set(SOURCE_ROOT  "${CMAKE_SOURCE_DIR}/xdiag")
  set(INCLUDE_ROOT "${CMAKE_SOURCE_DIR}")
  set(INST_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}")

  # For each grouped source: determine the generated file paths at configure
  # time, then register a build-time command to (re)generate them whenever the
  # original source changes.
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

  target_sources(${XDIAG_LIBRARY}
    PRIVATE
    ${XDIAG_NORMAL_SOURCES}
    ${XDIAG_INSTANTIATION_SOURCES})

endif()
