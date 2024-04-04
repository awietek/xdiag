set(XDIAG_SOURCES
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
  algorithms/time_evolution/exp_sym_v.cpp
  algorithms/time_evolution/time_evolution.cpp
  algorithms/time_evolution/pade_matrix_exponential.cpp
)

set(XDIAG_DISTRIBUTED_SOURCES
  parallel/mpi/allreduce.cpp
  parallel/mpi/alltoall.cpp
  parallel/mpi/communicator.cpp
  parallel/mpi/datatype.cpp
  parallel/mpi/cdot_distributed.cpp
  parallel/mpi/timing_mpi.cpp
  parallel/mpi/buffer.cpp

  basis/tj_distributed/basis_np.cpp
  blocks/tj_distributed/tj_distributed.cpp
  blocks/tj_distributed/tj_distributed_apply.cpp
)

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
  target_compile_options(${XDIAG_LIBRARY} PUBLIC -fp-model precise)
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

  
message(STATUS "-------- Determining Linear Algebra backend -------------")

###########################################################################
# try to get the Intel MKL if possible and link statically
find_package(MKL)
if(MKL_FOUND)
  message(STATUS "Employing IntelMKL as linear algebra backend")
  target_compile_definitions(${XDIAG_LIBRARY} PUBLIC XDIAG_USE_MKL)
  target_include_directories(${XDIAG_LIBRARY} PUBLIC ${MKL_INCLUDE_DIR})

  # the following lines are needed to link MKL statically correctly
  add_library(mkl_intel_ilp64 STATIC IMPORTED)
  add_library(mkl_sequential STATIC IMPORTED)
  add_library(mkl_core STATIC IMPORTED)
  
  set_target_properties(mkl_intel_ilp64 PROPERTIES
    IMPORTED_LOCATION ${MKL_INTERFACE_LIBRARY}
    IMPORTED_LINK_INTERFACE_MULTIPLICITY 3)
  
  set_target_properties(mkl_sequential PROPERTIES
    IMPORTED_LOCATION ${MKL_SEQUENTIAL_LAYER_LIBRARY}
    IMPORTED_LINK_INTERFACE_MULTIPLICITY 3)
  
  set_target_properties(mkl_core PROPERTIES
    IMPORTED_LOCATION ${MKL_CORE_LIBRARY}
    IMPORTED_LINK_INTERFACE_MULTIPLICITY 3)

  target_compile_definitions(mkl_intel_ilp64 INTERFACE MKL_ILP64)
  set_target_properties(mkl_sequential PROPERTIES INTERFACE_LINK_LIBRARIES mkl_core)
  set_target_properties(mkl_core       PROPERTIES INTERFACE_LINK_LIBRARIES mkl_sequential)

  foreach(mkl_lib IN LISTS mkl_intel_ilp6 mkl_sequential mkl_core)
    target_include_directories(mkl_lib PUBLIC ${MKL_INCLUDE_DIR})
  endforeach()  
  
  target_link_libraries(${XDIAG_LIBRARY} PUBLIC mkl_intel_ilp64 mkl_sequential mkl_core dl pthread)

  get_filename_component(MKL_LIB_DIR ${MKL_CORE_LIBRARY} DIRECTORY)
  target_link_directories(${XDIAG_LIBRARY} PUBLIC ${MKL_LIB_DIR})
  
###########################################################################
# Regular Blas/Lapack backend
else()
  find_package(LAPACK)
  
  if(LAPACK_FOUND)
    message(STATUS "Employing Lapack/Blas as linear algebra backend")
    target_link_libraries(${XDIAG_LIBRARY} PUBLIC ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
  else()
    message(FATAL_ERROR "No linear algebra backend found")
  endif()
endif()


###########################################################################
# OpenMP
if(XDIAG_DISABLE_OPENMP)
  message(STATUS "-------   OpenMP support has been disabled    -----------")
elseif(XDIAG_DISTRIBUTED)
  message(STATUS "------- OpenMP disabled for distributed library ---------")
else()
  message(STATUS "--------  Determining if OpenMP is present  -------------")
  find_package(OpenMP)
  
  if(OpenMP_CXX_FOUND)
    message(STATUS "Employing OpenMP")
    target_link_libraries(${XDIAG_LIBRARY} PUBLIC OpenMP::OpenMP_CXX)
  else()
    message(STATUS "OpenMP not found")
  endif()
endif()

###########################################################################
# HDF5
if(XDIAG_DISABLE_HDF5)
  message(STATUS "--------   HDF5 support has been disabled   -------------")
else()
  message(STATUS "--------   Determining if HDF5 is present   -------------")
  find_package(HDF5 COMPONENTS CXX)
  if(HDF5_FOUND)
    message(STATUS "Employing HDF5")
    target_compile_definitions(${XDIAG_LIBRARY} PUBLIC XDIAG_USE_HDF5 ${HDF5_CXX_DEFINITIONS})
    target_include_directories(${XDIAG_LIBRARY} PUBLIC ${HDF5_CXX_INCLUDE_DIRS})
    target_link_libraries(${XDIAG_LIBRARY} PUBLIC ${HDF5_CXX_LIBRARIES} ${HDF5_LIBRARIES})
  else()
    message(STATUS "HDF5 not found")
  endif()
endif()

###########################################################################
# creating a config header with source code information
message(STATUS "--------         Creating config.hpp        -------------")
find_package(Git REQUIRED)
execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
 WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
 OUTPUT_VARIABLE   XDIAG_GIT_HASH
 OUTPUT_STRIP_TRAILING_WHITESPACE
)

execute_process(COMMAND hostname 
    OUTPUT_VARIABLE XDIAG_HOSTNAME 
    RESULT_VARIABLE returncode 
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
execute_process(COMMAND id -un 
    OUTPUT_VARIABLE XDIAG_COMPILEDBY
    RESULT_VARIABLE returncode
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
configure_file(${CMAKE_SOURCE_DIR}/cmake/config.hpp.in ${CMAKE_SOURCE_DIR}/xdiag/config.hpp @ONLY)
