include(../cmake/test_sources.cmake)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/../cmake")

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
  add_compile_definitions(HYDRA_USE_HDF5)
endif()

add_executable(tests ${HYDRA_TEST_SOURCES})
target_link_libraries(tests PUBLIC hydra ${lapackblas_libraries} ${HDF5_CXX_LIBRARIES} ${HDF5_LIBRARIES})
