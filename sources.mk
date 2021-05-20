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

sources+= hydra/models/spinhalf/spinhalf.cpp
sources+= hydra/models/electron/electron.cpp
sources+= hydra/models/electron/electron_matrix.cpp
sources+= hydra/models/electron/electron_apply.cpp
sources+= hydra/models/electron/electron_utils.cpp
sources+= hydra/models/electron/electron_symmetric.cpp
sources+= hydra/models/electron/electron_symmetric_matrix.cpp
sources+= hydra/models/electron/electron_symmetric_apply.cpp
sources+= hydra/models/tj/tj.cpp
sources+= hydra/models/tj/tj_utils.cpp
sources+= hydra/models/tj/tj_matrix.cpp

sources+= hydra/symmetries/spacegroup.cpp
sources+= hydra/symmetries/spinflip.cpp
sources+= hydra/symmetries/representation.cpp
sources+= hydra/symmetries/spacegroup_operations.cpp
sources+= hydra/symmetries/spacegroup_operator.cpp
sources+= hydra/symmetries/fermi_sign.cpp

sources+= hydra/operators/bond.cpp
sources+= hydra/operators/bondlist.cpp
sources+= hydra/operators/couplings.cpp

# sources+= hydra/combinatorics/up_down_hole.cpp
# sources+= hydra/symmetries/charactertable.cpp

# sources+= hydra/utils/iochecks.cpp

# sources+= hydra/states/state_spinhalf.cpp
# sources+= hydra/states/state_electron.cpp
# sources+= hydra/states/state_tj.cpp

# sources+= hydra/bases/basis_spinhalf.cpp
# sources+= hydra/bases/basis_tj.cpp
# sources+= hydra/bases/basis_electron.cpp

# sources+= hydra/indexing/index_electron.cpp
# sources+= hydra/indexing/index_spinhalf.cpp

# # sources+= hydra/indexing/index_symmetrized.cpp
# # sources+= hydra/indexing/index_table.cpp


# sources+= hydra/models/heisenbergmodel.cpp
# sources+= hydra/models/hubbardmodel.cpp
# sources+= hydra/models/hubbardmodelmpi.cpp
# sources+= hydra/models/hubbardmodeldetail.cpp
# sources+= hydra/models/tjmodel.cpp
# sources+= hydra/models/tjmodelmpi.cpp

# # sources+= hydra/models/spinlessfermions.cpp

# sources+= hydra/mpi/communicator.cpp



# sources+= hydra/parameters/palm_exception.cpp


# sources+= hydra/entanglement/reduced_density_matrix.cpp
# sources+= hydra/entanglement/entanglement_entropy.cpp
# sources+= hydra/symmetries/symmetrydetail.cpp
# sources+= hydra/thermodynamics/thermodynamics_detail.cpp
# sources+= hydra/thermodynamics/thermodynamics_exact.cpp
# sources+= hydra/thermodynamics/thermodynamics_tpq.cpp


testsources+= test/tests.cpp

testsources+= test/combinatorics/test_binomial.cpp
testsources+= test/combinatorics/test_subsets.cpp
testsources+= test/combinatorics/test_bit_patterns.cpp
testsources+= test/combinatorics/test_combinations.cpp

testsources+= test/indexing/test_lintable.cpp

testsources+= test/symmetries/test_spacegroup_operator.cpp
testsources+= test/symmetries/test_fermi_sign.cpp

testsources+= test/models/test_spinhalf.cpp
# testsources+= test/models/test_electron.cpp
testsources+= test/models/testcases_electron.cpp
testsources+= test/models/test_electron_matrix.cpp
testsources+= test/models/test_electron_apply.cpp
testsources+= test/models/test_electron_symmetric.cpp
testsources+= test/models/test_electron_symmetric_matrix.cpp
testsources+= test/models/test_electron_symmetric_apply.cpp

testsources+= test/models/testcases_tj.cpp
testsources+= test/models/test_tj_matrix.cpp


# testsources+= test/symmetries/test_spacegroup.cpp
# testsources+= test/symmetries/test_representation.cpp
# testsources+= test/symmetries/test_charactertable.cpp
# testsources+= test/symmetries/test_symmetry_operations.cpp


# testsources+= test/test_indexhubbard.cpp

# testsources+= test/combinatorics/test_up_down_hole.cpp

# testsources+= test/qns/test_qns.cpp
# testsources+= test/states/test_states.cpp
# testsources+= test/operators/test_bondlist.cpp



# testsources+= test/bases/test_basis_spinhalf.cpp
# testsources+= test/bases/test_basis_tj.cpp
# testsources+= test/bases/test_basis_electron.cpp

# testsources+= test/indexing/test_index_table.cpp
# testsources+= test/indexing/test_index_search.cpp


# # testsources+= test/test_charactertable.cpp
# # testsources+= test/test_indexhubbard.cpp
# # testsources+= test/test_spacegroup.cpp

# testsources+= test/models/test_tjmodel.cpp
# testsources+= test/models/test_tjmodelmpi.cpp
# testsources+= test/models/test_hubbardmodel.cpp
# testsources+= test/models/test_hubbardmodelmpi.cpp

# testsources+= test/entanglement/test_entanglement_entropy.cpp
# testsources+= test/entanglement/test_entanglement_entropy_mpi.cpp
# testsources+= test/entanglement/test_reduced_density_matrix_mpi.cpp


# appsources+= hydra/applications/hubbardopticalftlm/hubbardopticalftlm.cpp
# appsources+= hydra/applications/hubbarddynamics/hubbarddynamics.cpp
# appsources+= hydra/applications/hubbarddynamicsmpi/hubbarddynamicsmpi.cpp


appsources+= hydra/applications/hubbarded/hubbarded.cpp
appsources+= hydra/applications/hubbardfulled/hubbardfulled.cpp
appsources+= hydra/applications/tjed/tjed.cpp
appsources+= hydra/applications/tjfulled/tjfulled.cpp
appsources+= hydra/applications/tjentanglement/tjentanglement.cpp
