sources+= hydra/utils/iochecks.cpp

sources+= hydra/combinatorics/binomial.cpp
sources+= hydra/combinatorics/subsets.cpp
sources+= hydra/combinatorics/bit_patterns.cpp
sources+= hydra/combinatorics/combinations.cpp
sources+= hydra/combinatorics/up_down_hole.cpp

sources+= hydra/parameters/parameter_value.cpp
sources+= hydra/parameters/parameters.cpp
sources+= hydra/parameters/parser.cpp

sources+= hydra/indexing/lintable.cpp
sources+= hydra/indexing/indexing_symmetric.cpp
sources+= hydra/indexing/indexing_symmetric_fermionic.cpp

sources+= hydra/models/utils/model_utils.cpp

sources+= hydra/models/spinhalf/spinhalf.cpp
sources+= hydra/models/spinhalf/spinhalf_matrix.cpp
sources+= hydra/models/spinhalf/spinhalf_apply.cpp

sources+= hydra/models/spinhalf_symmetric/spinhalf_symmetric.cpp
sources+= hydra/models/spinhalf_symmetric/spinhalf_symmetric_matrix.cpp
sources+= hydra/models/spinhalf_symmetric/spinhalf_symmetric_apply.cpp

sources+= hydra/models/tj/tj.cpp
sources+= hydra/models/tj/tj_utils.cpp
sources+= hydra/models/tj/tj_matrix.cpp
sources+= hydra/models/tj/tj_apply.cpp

sources+= hydra/models/tj_symmetric/tj_symmetric.cpp
sources+= hydra/models/tj_symmetric/tj_symmetric_matrix.cpp
sources+= hydra/models/tj_symmetric/tj_symmetric_apply.cpp

sources+= hydra/models/electron/electron.cpp
sources+= hydra/models/electron/electron_matrix.cpp
sources+= hydra/models/electron/electron_apply.cpp

sources+= hydra/models/electron_symmetric_simple/electron_symmetric_simple.cpp
sources+= hydra/models/electron_symmetric_simple/electron_symmetric_simple_matrix.cpp
sources+= hydra/models/electron_symmetric_simple/electron_symmetric_simple_apply.cpp

sources+= hydra/models/electron_symmetric/electron_symmetric.cpp
sources+= hydra/models/electron_symmetric/electron_symmetric_matrix.cpp
sources+= hydra/models/electron_symmetric/electron_symmetric_apply.cpp


sources+= hydra/symmetries/symmetry_utils.cpp
sources+= hydra/symmetries/permutation_group.cpp
sources+= hydra/symmetries/permutation_group_action.cpp
sources+= hydra/symmetries/permutation_group_lookup.cpp
sources+= hydra/symmetries/representation.cpp
sources+= hydra/symmetries/fermi_sign.cpp
sources+= hydra/symmetries/symmetric_operator.cpp

sources+= hydra/operators/bond.cpp
sources+= hydra/operators/bondlist.cpp
sources+= hydra/operators/couplings.cpp

sources+= hydra/wavefunctions/gpwf_spinhalf.cpp

sources+=hydra/linalg/lanczos/lanczos_convergence.cpp
sources+=hydra/linalg/lanczos/tmatrix.cpp

mpisources+= hydra/mpi/datatype.cpp
mpisources+= hydra/mpi/allreduce.cpp
mpisources+= hydra/mpi/alltoall.cpp
mpisources+= hydra/mpi/timing_mpi.cpp
mpisources+= hydra/mpi/stable_dot.cpp
mpisources+= hydra/mpi/communicator.cpp

mpisources+= hydra/models/model_utils_mpi.cpp

mpisources+= hydra/models/spinhalf_mpi/terms/get_prefix_postfix_mixed_bonds.cpp
mpisources+= hydra/models/spinhalf_mpi/spinhalf_mpi.cpp
mpisources+= hydra/models/spinhalf_mpi/spinhalf_mpi_apply.cpp

mpisources+= hydra/models/electron_mpi/electron_mpi.cpp
mpisources+= hydra/models/electron_mpi/electron_mpi_apply.cpp


testsources+= test/tests.cpp
testsources+= test/combinatorics/test_binomial.cpp
testsources+= test/combinatorics/test_subsets.cpp
testsources+= test/combinatorics/test_bit_patterns.cpp
testsources+= test/combinatorics/test_combinations.cpp

testsources+= test/indexing/test_lintable.cpp

testsources+= test/symmetries/test_permutation_group_action.cpp
testsources+= test/symmetries/test_permutation_group_lookup.cpp
testsources+= test/symmetries/test_fermi_sign.cpp
testsources+= test/symmetries/test_symmetric_operator.cpp

testsources+= test/models/spinhalf/testcases_spinhalf.cpp
testsources+= test/models/spinhalf/test_spinhalf_matrix.cpp
testsources+= test/models/spinhalf/test_spinhalf_apply.cpp

testsources+= test/models/spinhalf_symmetric/test_spinhalf_symmetric.cpp
testsources+= test/models/spinhalf_symmetric/test_spinhalf_symmetric_matrix.cpp
testsources+= test/models/spinhalf_symmetric/test_spinhalf_symmetric_apply.cpp

testsources+= test/models/tj/testcases_tj.cpp
testsources+= test/models/tj/test_tj_matrix.cpp
testsources+= test/models/tj/test_tj_apply.cpp

testsources+= test/models/tj_symmetric/test_tj_symmetric.cpp
testsources+= test/models/tj_symmetric/test_tj_symmetric_matrix.cpp
testsources+= test/models/tj_symmetric/test_tj_symmetric_apply.cpp

testsources+= test/models/electron/testcases_electron.cpp
testsources+= test/models/electron/test_electron_matrix.cpp
testsources+= test/models/electron/test_electron_apply.cpp

testsources+= test/models/electron_symmetric_simple/test_electron_symmetric_simple.cpp
testsources+= test/models/electron_symmetric_simple/test_electron_symmetric_simple_matrix.cpp
testsources+= test/models/electron_symmetric_simple/test_electron_symmetric_simple_apply.cpp

testsources+= test/models/electron_symmetric/test_electron_symmetric.cpp
testsources+= test/models/electron_symmetric/test_electron_symmetric_matrix.cpp
testsources+= test/models/electron_symmetric/test_electron_symmetric_apply.cpp

testsources+= test/linalg/lanczos/test_lanczos_generic.cpp
testsources+= test/linalg/lanczos/test_lanczos_eigenvalues.cpp
testsources+= test/linalg/lanczos/test_lanczos_eigenvector.cpp

testmpisources+= testmpi/tests.cpp
testmpisources+= testmpi/mpi/test_stable_dot.cpp
testmpisources+= test/models/spinhalf/testcases_spinhalf.cpp
testmpisources+= testmpi/models/spinhalf_mpi/test_spinhalf_mpi.cpp
testmpisources+= testmpi/models/spinhalf_mpi/test_spinhalf_mpi_apply.cpp


appsources+= hydra/applications/hubbarded/hubbarded.cpp
appsources+= hydra/applications/hubbardfulled/hubbardfulled.cpp
appsources+= hydra/applications/tjed/tjed.cpp
appsources+= hydra/applications/tjfulled/tjfulled.cpp
appsources+= hydra/applications/tjentanglement/tjentanglement.cpp
