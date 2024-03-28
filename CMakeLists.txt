cmake_minimum_required(VERSION 3.19)

project(
  hydra
  VERSION 1.0.0
  LANGUAGES CXX
)

option(HYDRA_DISTRIBUTED "Build the library for distibuted computing" OFF)
option(HYDRA_JULIA_WRAPPER "Build the Julia wrapper" OFF)
option(HYDRA_BUILD_TESTING "Build the tests" OFF)
option(HYDRA_DISABLE_THREADING "Determine whether library should be compiled without threading" OFF)
option(HYDRA_OPTIMIZE_FOR_NATIVE "Optimize for native architecture" OFF)

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
if (DEFINED HYDRA_SHARED_LIBS)
    set(BUILD_SHARED_LIBS "${HYDRA_SHARED_LIBS}")
endif ()

# determine which library is being built
if (HYDRA_JULIA_WRAPPER AND HYDRA_DISTRIBUTED)
  message(FATAL_ERROR "The julia wrapper cannot be build with distributed components")
elseif (HYDRA_JULIA_WRAPPER)
  set(HYDRA_LIBRARY "hydra")

  #library must be built as a shared library for the wrapper
  if (NOT HYDRA_SHARED_LIBS)
    message(MESSAGE "Forcing build of a shared library for the julia wrapper")
  endif()
  set(BUILD_SHARED_LIBS On)
elseif (HYDRA_DISTRIBUTED)
  set(HYDRA_LIBRARY "hydra_distributed")
else()
  set(HYDRA_LIBRARY "hydra")
endif()

# include main library
add_subdirectory(hydra)

## Generate the export header for SomeLib and attach it to the target
include(GenerateExportHeader)
generate_export_header(${HYDRA_LIBRARY} EXPORT_FILE_NAME include/hydra/export.h)
target_compile_definitions(
    ${HYDRA_LIBRARY} PUBLIC "$<$<NOT:$<BOOL:${BUILD_SHARED_LIBS}>>:HYDRA_STATIC_DEFINE>")
target_include_directories(
  ${HYDRA_LIBRARY} PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>")



#########################################################
# Testing
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    include(CTest)
endif()
enable_testing()
if((CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME OR HYDRA_BUILD_TESTING) AND BUILD_TESTING)
  message(STATUS "Building test programs")
  add_subdirectory(tests)
endif()

if(BUILD_JULIA)
  message(STATUS "Building Julia wrapper")
  set_property(TARGET hydra PROPERTY POSITION_INDEPENDENT_CODE ON)
  add_subdirectory(julia)
endif()

if(HYDRA_OPTIMIZE_FOR_NATIVE)
message(STATUS "Using native architecture optimizations")
target_compile_options(${HYDRA_LIBRARY} PRIVATE -march=native)
target_compile_options(${HYDRA_LIBRARY} PRIVATE -mtune=native)
endif()

# Install the headers
install(DIRECTORY "${CMAKE_SOURCE_DIR}/hydra" # source directory
  DESTINATION include
  FILES_MATCHING # install only matched files
  PATTERN "*.h" # select header files
)

install(DIRECTORY "${CMAKE_SOURCE_DIR}/hydra/extern" # source directory
  DESTINATION include/hydra
  FILES_MATCHING # install only matched files
  PATTERN "*" # select header files
)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

if (NOT DEFINED ${HYDRA_INSTALL_CMAKEDIR})
  set(HYDRA_INSTALL_CMAKEDIR "${CMAKE_INSTALL_LIBDIR}/cmake/${HYDRA_LIBRARY}"
    CACHE STRING "Path to ${HYDRA_LIBRARY} CMake files")
endif ()

install(TARGETS ${HYDRA_LIBRARY} EXPORT ${HYDRA_LIBRARY}_Targets
        RUNTIME COMPONENT ${HYDRA_LIBRARY}_Runtime
        LIBRARY COMPONENT ${HYDRA_LIBRARY}_Runtime
        NAMELINK_COMPONENT ${HYDRA_LIBRARY}_Development
        ARCHIVE COMPONENT ${HYDRA_LIBRARY}_Development
        INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")

if (BUILD_SHARED_LIBS)
    set(type shared)
else ()
    set(type static)
endif ()

install(EXPORT ${HYDRA_LIBRARY}_Targets
        DESTINATION "${HYDRA_INSTALL_CMAKEDIR}"
        NAMESPACE hydra::
        FILE ${HYDRA_LIBRARY}-${type}-targets.cmake
        COMPONENT ${HYDRA_LIBRARY}_Development)

write_basic_package_version_file(
    cmake/${HYDRA_LIBRARY}ConfigVersion.cmake
    COMPATIBILITY SameMajorVersion)

  install(FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/${HYDRA_LIBRARY}FindMKL.cmake"
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/${HYDRA_LIBRARY}Config.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake/${HYDRA_LIBRARY}ConfigVersion.cmake"
    DESTINATION "${HYDRA_INSTALL_CMAKEDIR}"
    COMPONENT ${HYDRA_LIBRARY}_Development)

# if(BUILD_JULIA)
#   install(TARGETS hydrajl
#     EXPORT hydrajlConfig
#     LIBRARY DESTINATION lib
#     ARCHIVE DESTINATION lib
#   )

#   install(EXPORT hydrajlConfig
#     FILE hydrajlConfig.cmake
#     DESTINATION lib/cmake/hydra
#   )
# endif()
