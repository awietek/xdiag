cmake_minimum_required(VERSION 3.1...3.25)

project(
  hydra
  VERSION 1.0
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-O3")

include(cmake/sources.cmake)
add_library(hydra STATIC ${HYDRA_SOURCES})
target_include_directories(hydra PUBLIC .)

add_subdirectory(test)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

# Create configure file
find_package(Git)

# Get the Git Hash
execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
 WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
 OUTPUT_VARIABLE   HYDRA_GIT_HASH
 OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Compute HYDRA_HOSTNAME and HYDRA_COMPILEDBY
execute_process(COMMAND hostname 
    OUTPUT_VARIABLE HYDRA_HOSTNAME 
    RESULT_VARIABLE returncode 
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
execute_process(COMMAND id -un 
    OUTPUT_VARIABLE HYDRA_COMPILEDBY
    RESULT_VARIABLE returncode
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

configure_file(./cmake/config.h.in ${PROJECT_SOURCE_DIR}/hydra/config.h @ONLY)

# Linear Algebra 
find_package(MKL)
if(MKL_FOUND)
  message("MKL found")
  add_compile_definitions(HYDRA_USE_MKL)
  target_include_directories(hydra PUBLIC ${MKL_INCLUDE_DIR})
  set(lapackblas_libraries ${MKL_LIBRARIES})
else()
  message("MKL NOT found, trying conventional Lapack / Blas")

  find_package(BLAS)
  find_package(LAPACK)
  if(LAPACK_FOUND AND BLAS_FOUND)
    message("Lapack / Blas found")
    set(lapackblas_libraries "${BLAS_LIBRARIES} ${LAPACK_LIBRARIES}")
  else()
    message("Lapack / Blas NOT found")
  endif()

endif()

# OpenMP
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  message("OpenMP found")
  target_link_libraries(hydra PUBLIC OpenMP::OpenMP_CXX)
else()
  message("OpenMP NOT found")
endif()

# Hdf5
find_package(HDF5 COMPONENTS CXX)
if (HDF5_FOUND) 
  message("Found hdf5 library")
  target_include_directories(hydra PUBLIC ${HDF5_CXX_INCLUDE_DIRS})
  add_compile_definitions(HYDRA_USE_HDF5 ${HDF5_CXX_DEFINITIONS})
endif()


