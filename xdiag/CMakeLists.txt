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
  target_compile_definitions(xdiag_distributed PUBLIC XDIAG_USE_MPI)
else()
  message(STATUS "--------    BUILDING NORMAL XDIAG LIBRARY   -------------")
  target_sources(${XDIAG_LIBRARY} PRIVATE ${XDIAG_SOURCES})
endif()
