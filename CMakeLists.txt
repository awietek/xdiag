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
  option(MKL_THREADING "Enable IntelMKL threading" ON)
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

if(BUILD_TESTS)
  message(STATUS "Building test programs")
  enable_testing()
  add_subdirectory(test)
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


install(DIRECTORY "${CMAKE_SOURCE_DIR}/hydra" # source directory
  DESTINATION include
  FILES_MATCHING # install only matched files
  PATTERN "*.h" # select header files
)

if(BUILD_DISTRIBUTED)
  install(TARGETS hydra_distributed
    EXPORT hydra_distributedConfig
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
  )

  include(CMakePackageConfigHelpers)	

  # generate the config file that includes the exports
  configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/hydraConfig.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/hydra_distributedConfig.cmake"
    INSTALL_DESTINATION "lib/cmake/hydra_distributed"
    NO_SET_AND_CHECK_MACRO
    NO_CHECK_REQUIRED_COMPONENTS_MACRO
  )
  
  write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/hydra_distributedConfigVersion.cmake"
    VERSION "${hydra_distributed_VERSION_MAJOR}.${hydra_distributed_VERSION_MINOR}"
    COMPATIBILITY AnyNewerVersion
  )

  install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/hydra_distributedConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/hydra_distributedConfigVersion.cmake
    DESTINATION lib/cmake/hydra_distributed
  )

  export(EXPORT hydra_distributedConfig
    FILE "${CMAKE_CURRENT_BINARY_DIR}/hydra_distributedConfig.cmake"
  )
else()
  install(TARGETS hydra
    EXPORT hydraConfig
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
  )

  include(CMakePackageConfigHelpers)
  
  # generate the config file that includes the exports
  configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/hydraConfig.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/hydraConfig.cmake"
    INSTALL_DESTINATION "lib/cmake/hydra"
    NO_SET_AND_CHECK_MACRO
    NO_CHECK_REQUIRED_COMPONENTS_MACRO
  )
  
  target_include_directories(hydra
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  )
  
  write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/hydraConfigVersion.cmake"
    VERSION "${hydra_VERSION_MAJOR}.${hydra_VERSION_MINOR}"
    COMPATIBILITY AnyNewerVersion
  )

  install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/hydraConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/hydraConfigVersion.cmake
    DESTINATION lib/cmake/hydra
  )

  export(EXPORT hydraConfig
    FILE "${CMAKE_CURRENT_BINARY_DIR}/hydraConfig.cmake"
  )

endif()

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
