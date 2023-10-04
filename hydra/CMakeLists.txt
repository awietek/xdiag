list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

set(HYDRA_DEFINITIONS "")
set(HYDRA_INCLUDE_DIRS "${CMAKE_SOURCE_DIR}")
set(HYDRA_LIBRARIES "${PROJECT_BINARY_DIR}/hydra/libhydra.a")

set(HYDRA_SOURCES
  utils/iochecks.cpp
  utils/print.cpp
  utils/error.cpp
  utils/say_hello.cpp
  utils/read_vectors.cpp
  bits/bitops.cpp
  
  parallel/omp/omp_utils.cpp  
  algebra/algebra.cpp
  algebra/matrix.cpp
  algebra/apply.cpp

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
  combinatorics/combinations_indexing.cpp
  combinatorics/subsets_indexing.cpp
  combinatorics/lin_table.cpp
  combinatorics/fermi_table.cpp
    
  basis/basis.cpp

  basis/spinhalf/basis_sz.cpp
  basis/spinhalf/basis_no_sz.cpp
  basis/spinhalf/basis_symmetric_sz.cpp
  basis/spinhalf/basis_symmetric_no_sz.cpp
  basis/spinhalf/basis_sublattice.cpp

  basis/tj/basis_np.cpp
  basis/tj/basis_symmetric_np.cpp

  basis/electron/basis_np.cpp
  basis/electron/basis_no_np.cpp
  basis/electron/basis_symmetric_np.cpp
  basis/electron/basis_symmetric_no_np.cpp

  blocks/blocks.cpp

  blocks/spinhalf/spinhalf.cpp
  blocks/spinhalf/spinhalf_matrix.cpp
  blocks/spinhalf/spinhalf_apply.cpp
  blocks/spinhalf/compile.cpp
  blocks/spinhalf/qns.cpp

  blocks/electron/electron.cpp
  blocks/electron/electron_matrix.cpp
  blocks/electron/electron_apply.cpp
  blocks/electron/compile.cpp

  blocks/tj/tj.cpp
  blocks/tj/tj_matrix.cpp
  blocks/tj/tj_apply.cpp
  blocks/tj/compile.cpp

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

  states/gpwf.cpp
  states/product_state.cpp
  states/random_state.cpp
  states/state.cpp
  
  random/random_utils.cpp
  random/hash.cpp

  algorithms/lanczos/lanczos_convergence.cpp
  algorithms/lanczos/tmatrix.cpp
  algorithms/lanczos/eigvals_lanczos.cpp
  algorithms/lanczos/eigs_lanczos.cpp
  algorithms/sparse_diag.cpp
  algorithms/arnoldi/arnoldi_to_disk.cpp
  algorithms/gram_schmidt/gram_schmidt.cpp
  algorithms/gram_schmidt/orthogonalize.cpp

  algorithms/norm_estimate.cpp
  algorithms/time_evolution/zahexpv.cpp
  algorithms/time_evolution/exp_sym_v.cpp
  algorithms/time_evolution/time_evolution.cpp
  algorithms/time_evolution/pade_matrix_exponential.cpp
)

set(HYDRA_DISTRIBUTED_SOURCES
  parallel/mpi/allreduce.cpp
  parallel/mpi/alltoall.cpp
  parallel/mpi/buffer.cpp
  parallel/mpi/communicator.cpp
  parallel/mpi/datatype.cpp
  parallel/mpi/cdot_distributed.cpp
  parallel/mpi/timing_mpi.cpp
  parallel/mpi/logger_mpi.cpp
)

if(BUILD_DISTRIBUTED)
  message(STATUS "-------- BUILDING DISTRIBUTED HYDRA LIBRARY -------------")
  add_library(hydra_distributed STATIC ${HYDRA_SOURCES} ${HYDRA_DISTRIBUTED_SOURCES})
  target_compile_definitions(hydra_distributed PUBLIC HYDRA_USE_MPI)
  set(HYDRA_LIBRARY hydra_distributed)
else()
  message(STATUS "-------- BUILDING NORMAL HYDRA LIBRARY -------------")
  add_library(hydra STATIC ${HYDRA_SOURCES})
  set(HYDRA_LIBRARY hydra)
endif()

target_include_directories(${HYDRA_LIBRARY} PUBLIC ${CMAKE_SOURCE_DIR})
target_compile_features(${HYDRA_LIBRARY} PUBLIC cxx_std_17)

if(OPTIMIZE_FOR_NATIVE)
message(STATUS "Using native architecture optimizations")
target_compile_options(${HYDRA_LIBRARY} PRIVATE -march=native)
target_compile_options(${HYDRA_LIBRARY} PRIVATE -mtune=native)
endif()

target_compile_options(${HYDRA_LIBRARY} PRIVATE -Wall -Wextra)
  
message(STATUS "-------- Determining Linear Algebra backend -------------")
find_package(MKL)
if(MKL_FOUND)
  message(STATUS "Employing IntelMKL as linear algebra backend")

  target_compile_definitions(${HYDRA_LIBRARY} PUBLIC HYDRA_USE_MKL MKL_ILP64)
  list(APPEND HYDRA_DEFINITIONS HYDRA_USE_MKL MKL_ILP64)

  target_include_directories(${HYDRA_LIBRARY} PUBLIC ${MKL_INCLUDE_DIR})
  list(APPEND HYDRA_INCLUDE_DIRS ${MKL_INCLUDE_DIR})


  # Choose correct MKL libraries
  if(OMP_THREADING)

    target_compile_options(${HYDRA_LIBRARY} PUBLIC -fopenmp)
    
    if(MKL_THREADING)
      if(CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
	set(MKL_LIBRARIES ${MKL_INTERFACE_LIBRARY} ${MKL_THREAD_LAYER_LIBRARY} ${MKL_CORE_LIBRARY})
      elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
	set(MKL_LIBRARIES ${MKL_INTERFACE_LIBRARY} ${MKL_THREAD_LAYER_LIBRARY} ${MKL_CORE_LIBRARY})
      else()
	message(STATUS "OMP and MKL threading only works for GNU and Intel compilers")
	set(MKL_LIBRARIES ${MKL_INTERFACE_LIBRARY} ${MKL_SEQUENTIAL_LAYER_LIBRARY} ${MKL_CORE_LIBRARY})
      endif()
    else()
      set(MKL_LIBRARIES ${MKL_INTERFACE_LIBRARY} ${MKL_SEQUENTIAL_LAYER_LIBRARY} ${MKL_CORE_LIBRARY})
    endif()

  else()
    if(MKL_THREADING)
      set(MKL_LIBRARIES ${MKL_INTERFACE_LIBRARY} ${MKL_THREAD_LAYER_LIBRARY} ${MKL_CORE_LIBRARY})
    else()
      set(MKL_LIBRARIES ${MKL_INTERFACE_LIBRARY} ${MKL_SEQUENTIAL_LAYER_LIBRARY} ${MKL_CORE_LIBRARY})
    endif()
  
  endif()

  if(NOT WIN32 AND NOT APPLE)
    # Added -Wl block to avoid circular dependencies.
    # https://stackoverflow.com/questions/5651869/what-are-the-start-group-and-end-group-command-line-options
    # https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor
    set(MKL_LIBRARIES -Wl,--start-group ${MKL_LIBRARIES} -Wl,--end-group -liomp5 -lpthread -lm -ldl)
  elseif(APPLE)
    # MacOS does not support --start-group and --end-group
    set(MKL_LIBRARIES -Wl,${MKL_LIBRARIES} -Wl, -liomp5 -lpthread -lm -ldl)
  endif()
  
  message(STATUS "MKL libraries: " ${MKL_LIBRARIES})
  target_link_libraries(${HYDRA_LIBRARY} PUBLIC ${MKL_LIBRARIES})
  list(APPEND HYDRA_LIBRARIES ${MKL_LIBRARIES})

# Regular Blas/Lapack backend
else()
  find_package(BLAS)
  find_package(LAPACK)
  
  if(LAPACK_FOUND AND BLAS_FOUND)
    message(STATUS "Emploging Lapack/Blas as linear algebra backend")
    target_link_libraries(${HYDRA_LIBRARY} PUBLIC ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
    list(APPEND HYDRA_LIBRARIES ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
  else()
    message(FATAL_ERROR "No linear algebra backend found")
  endif()

  if(OMP_THREADING)
    message(STATUS "--------  Determining if OpenMP is present  -------------")
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
      message(STATUS "Employing OpenMP")
      target_link_libraries(${HYDRA_LIBRARY} PUBLIC OpenMP::OpenMP_CXX)
      list(APPEND HYDRA_LIBRARIES OpenMP::OpenMP_CXX)
    endif()
  else()
    message(STATUS "OpenMP disabled")
  endif()
  
endif()

message(STATUS "--------   Determining if HDF5 is present   -------------")
find_package(HDF5 COMPONENTS CXX)
if(HDF5_FOUND)
  message(STATUS "Employing HDF5")
  target_compile_definitions(${HYDRA_LIBRARY} PUBLIC HYDRA_USE_HDF5 ${HDF5_CXX_DEFINITIONS})
  target_include_directories(${HYDRA_LIBRARY} PUBLIC ${HDF5_CXX_INCLUDE_DIRS})
  target_link_libraries(${HYDRA_LIBRARY} PUBLIC ${HDF5_CXX_LIBRARIES} ${HDF5_LIBRARIES})

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
