cmake_minimum_required(VERSION 3.19)

project(
  xdiag
  VERSION 0.2.0
  LANGUAGES CXX
)

option(BUILD_TESTING "Build the tests" Off)
option(XDIAG_DISTRIBUTED "Build the distibuted parallelization libraries" Off)
option(XDIAG_JULIA_WRAPPER "Build the Julia wrapper" Off)
option(XDIAG_DISABLE_OPENMP "Disables the library being compiled with OpenMP" Off)
option(XDIAG_DISABLE_HDF5 "Disables the library being compiled with HDF5" Off)
option(XDIAG_OPTIMIZE_FOR_NATIVE "Optimize for native architecture" Off)

message(STATUS "Compiler ID     : " ${CMAKE_CXX_COMPILER_ID})
message(STATUS "Compiler version: " ${CMAKE_CXX_COMPILER_VERSION})

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Set release build as default
if(NOT CMAKE_BUILD_TYPE)
 set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})

#########################################################
# Set a different default install prefix than system dirs
if(DEFINED CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    message(
        STATUS
        "CMAKE_INSTALL_PREFIX is not set\n"
        "Default value: ${CMAKE_INSTALL_PREFIX}\n"
        "Will set it to ${CMAKE_SOURCE_DIR}/install"
    )
    set(CMAKE_INSTALL_PREFIX
        "${CMAKE_SOURCE_DIR}/install"
        CACHE PATH "Where the library will be installed to" FORCE
    )
endif()

  
###############################################################
# On windows all symbols should be exported 
if (MSVC)
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()

# The follow lines would be necessary to deal with visibility
# of a shared library,
# https://alexreinking.com/blog/building-a-dual-shared-and-static-library-with-cmake.html

# if (NOT DEFINED CMAKE_CXX_VISIBILITY_PRESET AND
#     NOT DEFINED CMAKE_VISIBILITY_INLINES_HIDDEN)
#   set(CMAKE_CXX_VISIBILITY_PRESET hidden)
#   set(CMAKE_VISIBILITY_INLINES_HIDDEN YES)
# endif ()

## Let SomeLib_SHARED_LIBS override BUILD_SHARED_LIBS
if (DEFINED XDIAG_SHARED_LIBS)
    set(BUILD_SHARED_LIBS "${XDIAG_SHARED_LIBS}")
endif ()

# determine which library is being built
if (XDIAG_JULIA_WRAPPER AND XDIAG_DISTRIBUTED)
  message(FATAL_ERROR "The julia wrapper cannot be build with distributed components")
elseif (XDIAG_JULIA_WRAPPER)
  set(XDIAG_LIBRARY "xdiag")

  #library must be built as a shared library for the wrapper
  if (NOT XDIAG_SHARED_LIBS)
    message(STATUS "Forcing build of a shared library for the julia wrapper")
  endif()
  set(BUILD_SHARED_LIBS On)
elseif (XDIAG_DISTRIBUTED)
  set(XDIAG_LIBRARY "xdiag_distributed")
else()
  set(XDIAG_LIBRARY "xdiag")
endif()

# include main library
add_subdirectory(xdiag)


## would be needed for visibility contron in shared objects ->

# # Generate the export header for SomeLib and attach it to the target
# include(GenerateExportHeader)
# generate_export_header(${XDIAG_LIBRARY} EXPORT_FILE_NAME include/xdiag/export.h)
# target_compile_definitions(
#     ${XDIAG_LIBRARY} PUBLIC "$<$<NOT:$<BOOL:${BUILD_SHARED_LIBS}>>:XDIAG_STATIC_DEFINE>")
# target_include_directories(
#   ${XDIAG_LIBRARY} PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>")

# Add sanitizers to debug mode
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
  message(STATUS "Adding Sanitizers to debug build")
  target_compile_options(${XDIAG_LIBRARY} PUBLIC -fsanitize=address)
  target_link_options(${XDIAG_LIBRARY} PUBLIC -fsanitize=address)
endif()
  
#########################################################
# Testing
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    include(CTest)
endif()

if (BUILD_TESTING)
  enable_testing()
  if((CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME) AND BUILD_TESTING)
    message(STATUS "-----------    Building test programs      --------------")
    add_subdirectory(tests)
  endif()
endif()

if(XDIAG_JULIA_WRAPPER)
  message(STATUS "-----------    Building Julia wrapper       -------------")
  add_subdirectory(julia)
endif()

if(XDIAG_OPTIMIZE_FOR_NATIVE)
message(STATUS "Using native architecture optimizations")
target_compile_options(${XDIAG_LIBRARY} PRIVATE -march=native)
target_compile_options(${XDIAG_LIBRARY} PRIVATE -mtune=native)
endif()

# Install the headers
install(DIRECTORY "${CMAKE_SOURCE_DIR}/xdiag" # source directory
  DESTINATION include
  FILES_MATCHING # install only matched files
  PATTERN "*.hpp" # select header files
)

install(DIRECTORY "${CMAKE_SOURCE_DIR}/xdiag/extern" # source directory
  DESTINATION include/xdiag
  FILES_MATCHING # install only matched files
  PATTERN "*" # select header files
)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

if (NOT DEFINED ${XDIAG_INSTALL_CMAKEDIR})
  set(XDIAG_INSTALL_CMAKEDIR "${CMAKE_INSTALL_LIBDIR}/cmake/${XDIAG_LIBRARY}"
    CACHE STRING "Path to ${XDIAG_LIBRARY} CMake files")
endif ()

install(TARGETS ${XDIAG_LIBRARY} EXPORT ${XDIAG_LIBRARY}_Targets
        RUNTIME COMPONENT ${XDIAG_LIBRARY}_Runtime
        LIBRARY COMPONENT ${XDIAG_LIBRARY}_Runtime
        NAMELINK_COMPONENT ${XDIAG_LIBRARY}_Development
        ARCHIVE COMPONENT ${XDIAG_LIBRARY}_Development
        INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")

if (BUILD_SHARED_LIBS)
    set(type shared)
else ()
    set(type static)
endif ()

install(EXPORT ${XDIAG_LIBRARY}_Targets
        DESTINATION "${XDIAG_INSTALL_CMAKEDIR}"
        NAMESPACE xdiag::
        FILE ${XDIAG_LIBRARY}-${type}-targets.cmake
        COMPONENT ${XDIAG_LIBRARY}_Development)

write_basic_package_version_file(
    cmake/${XDIAG_LIBRARY}ConfigVersion.cmake
    COMPATIBILITY SameMajorVersion)

install(FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/${XDIAG_LIBRARY}Config.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake/${XDIAG_LIBRARY}ConfigVersion.cmake"
    DESTINATION "${XDIAG_INSTALL_CMAKEDIR}"
    COMPONENT ${XDIAG_LIBRARY}_Development)

if(XDIAG_JULIA_WRAPPER)

  install(TARGETS xdiagjl
    EXPORT xdiagjl_Targets
    LIBRARY COMPONENT ${XDIAG_LIBRARY}_Runtime
    NAMELINK_COMPONENT ${XDIAG_LIBRARY}_Development
    ARCHIVE COMPONENT ${XDIAG_LIBRARY}_Development
    INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
  )
 
endif()
