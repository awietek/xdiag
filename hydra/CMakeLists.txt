list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

set(HYDRA_DEFINITIONS "")
set(HYDRA_INCLUDE_DIRS "${CMAKE_SOURCE_DIR}")
set(HYDRA_LIBRARIES "${PROJECT_BINARY_DIR}/hydra/libhydra.a")

add_library(hydra STATIC
  utils/iochecks.cpp
  utils/print.cpp
  utils/say_hello.cpp
  bitops/bitops.cpp
  
  parallel/omp/omp_utils.cpp
  
  algebra/algebra.cpp
  algebra/matrix.cpp

  io/args.cpp
  io/args_handler.cpp
  io/file_toml.cpp
  io/toml/file_toml_handler.cpp
  io/toml/toml_conversion.cpp
  io/file_h5.cpp
  io/hdf5/file_h5_handler.cpp
  io/hdf5/utils.cpp
  io/hdf5/write.cpp
  io/hdf5/types.cpp
  
  combinatorics/binomial.cpp
  combinatorics/subsets.cpp
  combinatorics/subsets_index.cpp
  combinatorics/bit_patterns.cpp
  combinatorics/combinations.cpp
  combinatorics/combinations_index.cpp
  
  indexing/indexing_variants.cpp
  indexing/lin_table.cpp
  indexing/fermi_table.cpp
  indexing/combinations_indexing.cpp
  indexing/subsets_indexing.cpp
  
  indexing/spinhalf/indexing_sz.cpp
  indexing/spinhalf/indexing_no_sz.cpp
  indexing/spinhalf/indexing_symmetric_sz.cpp
  indexing/spinhalf/indexing_symmetric_no_sz.cpp
  indexing/spinhalf/indexing_sublattice.cpp
  indexing/spinhalf/symmetric_iterator.cpp

  indexing/tj/indexing_np.cpp
  indexing/tj/indexing_symmetric_np.cpp

  indexing/electron/indexing_np.cpp
  indexing/electron/indexing_no_np.cpp
  indexing/electron/indexing_symmetric_np.cpp
  indexing/electron/indexing_symmetric_no_np.cpp

  blocks/blocks.cpp
  blocks/utils/block_utils.cpp
  blocks/spinhalf/spinhalf.cpp
  blocks/spinhalf/spinhalf_matrix.cpp
  blocks/spinhalf/spinhalf_apply.cpp
  blocks/spinhalf/terms/compile.cpp
  blocks/spinhalf/terms/qns.cpp

  blocks/electron/electron.cpp
  blocks/electron/electron_matrix.cpp
  blocks/electron/electron_apply.cpp
  blocks/electron/terms/compile.cpp

  blocks/tj/tj.cpp
  blocks/tj/tj_matrix.cpp
  blocks/tj/tj_apply.cpp
  blocks/tj/terms/compile.cpp

  symmetries/qn.cpp
  symmetries/operations/symmetry_operations.cpp	
  symmetries/permutation.cpp
  symmetries/permutation_group.cpp
  symmetries/generated_group.cpp
  symmetries/representation.cpp
  symmetries/operations/fermi_sign.cpp

  symmetries/group_action/group_action.cpp
  symmetries/group_action/group_action_lookup.cpp
  symmetries/group_action/group_action_sublattice.cpp
  symmetries/group_action/sublattice_stability.cpp

  operators/bond.cpp
  operators/bondlist.cpp
  operators/bondlist_handler.cpp
  operators/compiler.cpp
  operators/symmetrized_operator.cpp
  operators/non_branching_bonds.cpp

  states/gpwf_spinhalf.cpp
  states/product_state.cpp
  states/random_state.cpp
  states/state.cpp

  random/random_utils.cpp
  random/hashes.cpp

  algorithms/lanczos/lanczos_convergence.cpp
  algorithms/lanczos/tmatrix.cpp
  algorithms/lanczos/lanczos_eigenvector.cpp
  algorithms/lanczos/lanczos_eigenvalues.cpp
  algorithms/sparse_diag.cpp
  algorithms/norm_estimate.cpp
  algorithms/time_evolution/zahexpv.cpp
  algorithms/time_evolution/time_evolution.cpp
  algorithms/time_evolution/pade_matrix_exponential.cpp
)


target_include_directories(hydra PUBLIC ${CMAKE_SOURCE_DIR})
target_compile_features(hydra PUBLIC cxx_std_17)

if(OPTIMIZE_FOR_NATIVE)
message(STATUS "Using -march=native optimizations")
target_compile_options(hydra PRIVATE -march=native)
target_compile_options(hydra PRIVATE -mtune=native)
endif()

target_compile_options(hydra PRIVATE -Wall -Wextra)
  
message(STATUS "-------- Determining Linear Algebra backend -------------")
find_package(MKL)
find_package(BLAS)
find_package(LAPACK)
if(MKL_FOUND)
  message(STATUS "Employing IntelMKL as linear algebra backend")
  target_compile_definitions(hydra PUBLIC HYDRA_USE_MKL MKL_ILP64)
  target_include_directories(hydra PUBLIC ${MKL_INCLUDE_DIR})
  target_link_libraries(hydra PUBLIC ${MKL_LIBRARIES})

  list(APPEND HYDRA_DEFINITIONS HYDRA_USE_MKL MKL_ILP64)
  list(APPEND HYDRA_INCLUDE_DIRS ${MKL_INCLUDE_DIR})
  list(APPEND HYDRA_LIBRARIES ${MKL_LIBRARIES})
  
elseif(LAPACK_FOUND AND BLAS_FOUND)
  message(STATUS "Emploging Lapack/Blas as linear algebra backend")
  target_link_libraries(hydra PUBLIC ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

  list(APPEND HYDRA_LIBRARIES ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
else()
  message(FATAL_ERROR "No linear algebra backend found")
endif()


message(STATUS "--------  Determining if OpenMP is present  -------------")
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  message(STATUS "Employing OpenMP")
  target_link_libraries(hydra PUBLIC OpenMP::OpenMP_CXX)
  list(APPEND HYDRA_LIBRARIES OpenMP::OpenMP_CXX)
endif()


message(STATUS "--------   Determining if HDF5 is present   -------------")
find_package(HDF5 COMPONENTS CXX)
if(HDF5_FOUND)
  message(STATUS "Employing HDF5")
  target_compile_definitions(hydra PUBLIC HYDRA_USE_HDF5 ${HDF5_CXX_DEFINITIONS})
  target_include_directories(hydra PUBLIC ${HDF5_CXX_INCLUDE_DIRS})
  target_link_libraries(hydra PUBLIC ${HDF5_CXX_LIBRARIES} ${HDF5_LIBRARIES})

  list(APPEND HYDRA_DEFINITIONS HYDRA_USE_HDF5 ${HDF5_CXX_DEFINITIONS})
  list(APPEND HYDRA_INCLUDE_DIRS ${HDF5_CXX_INCLUDE_DIRS})
  list(APPEND HYDRA_LIBRARIES ${HDF5_CXX_LIBRARIES} ${HDF5_LIBRARIES})
endif()

message(STATUS "--------         Creating config.h          -------------")
find_package(Git REQUIRED)
execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
 WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
 OUTPUT_VARIABLE   HYDRA_GIT_HASH
 OUTPUT_STRIP_TRAILING_WHITESPACE
)
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
configure_file(${CMAKE_SOURCE_DIR}/cmake/config.h.in ${CMAKE_SOURCE_DIR}/hydra/config.h @ONLY)
configure_file(${CMAKE_SOURCE_DIR}/cmake/hydraConfig.cmake.in ${CMAKE_SOURCE_DIR}/cmake/hydraConfig.cmake @ONLY)
