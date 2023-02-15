list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}")

message( STATUS "-------- Determining Linear Algebra backend -------------")
# Linear Algebra 
find_package(MKL)
if(MKL_FOUND)
  message("MKL found")
  add_compile_definitions(HYDRA_USE_MKL)
  target_include_directories(hydra PUBLIC ${MKL_INCLUDE_DIR})
  set(lapack_libraries ${MKL_LIBRARIES})
else()
  message("MKL NOT found, trying conventional Lapack / Blas")

  find_package(BLAS)
  find_package(LAPACK)
  if(LAPACK_FOUND AND BLAS_FOUND)
    message("Lapack / Blas found")
    set(lapack_libraries "${BLAS_LIBRARIES} ${LAPACK_LIBRARIES}")
  else()
    message(FATAL_ERROR "No linear algebra backend (e.g. BLAS/Lapack) could be found")
  endif()
endif()
target_link_libraries(hydra PUBLIC ${lapack_libraries})


message( STATUS "-------- Determining if OpenMP is present -------------")
# OpenMP
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  message("OpenMP found")
  target_link_libraries(hydra PUBLIC OpenMP::OpenMP_CXX)
else()
  message("OpenMP NOT found")
endif()

message( STATUS "-------- Determining if HDF5 is present -------------")
# Hdf5
find_package(HDF5 COMPONENTS CXX)
if (HDF5_FOUND) 
  message("HDF5 found")
  target_include_directories(hydra PUBLIC ${HDF5_CXX_INCLUDE_DIRS})
  add_compile_definitions(HYDRA_USE_HDF5 ${HDF5_CXX_DEFINITIONS})
  set(hdf5_libraries ${HDF5_CXX_LIBRARIES} ${HDF5_LIBRARIES})
else()
  message("HDF5 NOT found")
endif()
