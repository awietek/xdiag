cmake_minimum_required(VERSION 3.15)

project(
  hydra
  VERSION 1.0
  LANGUAGES CXX
)

option(BUILD_TESTS "Build test programs" OFF)
option(BUILD_JULIA "Build the Julia wrapper" OFF)
option(BUILD_DISTRIBUTED "Build the library for distibuted computing" OFF)

option(OPTIMIZE_FOR_NATIVE "Optimize for native architecture" OFF)


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
else()
    message(
        STATUS
        "CMAKE_INSTALL_PREFIX was already set\n"
        "Current value: ${CMAKE_INSTALL_PREFIX}"
    )
endif()
  

if(BUILD_DISTRIBUTED)
  find_package(MPI REQUIRED)
  option(OMP_THREADING "Enable OpenMp threading" OFF)
  option(MKL_THREADING "Enable IntelMKL threading" OFF)
  if(BUILD_JULIA)
    message(FATAL_ERROR "Option clash, cannot build Julia library with Distributed, disable either BUILD_JULIA or BUILD_DISTRIBUTED")
  endif()
else()
  option(OMP_THREADING "Enable OpenMp threading" ON)
  option(MKL_THREADING "Enable IntelMKL threading" OFF)
endif()
option(SANITIZER "Build with -fsanitize=address" OFF)

message(STATUS "Compiler ID     : " ${CMAKE_CXX_COMPILER_ID})
message(STATUS "Compiler version: " ${CMAKE_CXX_COMPILER_VERSION})

# Set release build as default
if(NOT CMAKE_BUILD_TYPE)
 set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})

add_subdirectory(hydra)

target_compile_features(${HYDRA_LIBRARY} PUBLIC cxx_std_17)


if(BUILD_TESTING)
  message(STATUS "Building test programs")
  enable_testing()
  add_subdirectory(tests)
endif()

if(BUILD_JULIA)
  message(STATUS "Building Julia wrapper")
  set_property(TARGET hydra PROPERTY POSITION_INDEPENDENT_CODE ON)
  add_subdirectory(julia)
endif()

if(SANITIZER)
    message(STATUS "Compiling with sanitizers")
    if(BUILD_DISTRIBUTED)
	target_compile_options(hydra_distributed PUBLIC -fsanitize=address)
    	target_link_options(hydra_distributed PUBLIC -fsanitize=address)
    else()			      
    	target_compile_options(hydra PUBLIC -fsanitize=address)
    	target_link_options(hydra PUBLIC -fsanitize=address)
    endif()			      
endif()

if(OPTIMIZE_FOR_NATIVE)
message(STATUS "Using native architecture optimizations")
target_compile_options(${HYDRA_LIBRARY} PRIVATE -march=native)
target_compile_options(${HYDRA_LIBRARY} PRIVATE -mtune=native)
endif()

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

# install the target and create export-set
install(TARGETS ${HYDRA_LIBRARY}
    EXPORT "${HYDRA_LIBRARY}Targets"
    # these get default values from GNUInstallDirs, no need to set them
    #RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} # bin
    #LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} # lib
    #ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} # lib
    # except for public headers, as we want them to be inside a library folder
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} # include
)

install(EXPORT "${HYDRA_LIBRARY}Targets"
  FILE "${HYDRA_LIBRARY}Targets.cmake"
  # NAMESPACE ${namespace}::
  DESTINATION lib/cmake/${HYDRA_LIBRARY}
)

include(CMakePackageConfigHelpers)	

# generate the version file for the config file
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/${HYDRA_LIBRARY}ConfigVersion.cmake"
    VERSION "${version}"
    COMPATIBILITY AnyNewerVersion
)
# create config file
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/${HYDRA_LIBRARY}Config.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/${HYDRA_LIBRARY}Config.cmake"
    INSTALL_DESTINATION cmake
)

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/${HYDRA_LIBRARY}Config.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/${HYDRA_LIBRARY}ConfigVersion.cmake
  DESTINATION lib/cmake/${HYDRA_LIBRARY}
)

export(EXPORT ${HYDRA_LIBRARY}Targets
  FILE "${CMAKE_CURRENT_BINARY_DIR}/cmake/${HYDRA_LIBRARY}Targets.cmake"
  # NAMESPACE ${namespace}::
)

if(BUILD_JULIA)
  install(TARGETS hydrajl
    EXPORT hydrajlConfig
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
  )

  install(EXPORT hydrajlConfig
    FILE hydrajlConfig.cmake
    DESTINATION lib/cmake/hydra
  )
endif()
