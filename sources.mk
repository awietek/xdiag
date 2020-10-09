sources+= hydra/utils/combinatorics.cpp
sources+= hydra/utils/iochecks.cpp

sources+= hydra/states/state_spinhalf.cpp
sources+= hydra/states/state_electron.cpp
sources+= hydra/states/state_tj.cpp

sources+= hydra/bases/basis_electron.cpp
sources+= hydra/bases/basis_spinhalf.cpp

sources+= hydra/indexing/index_electron.cpp
sources+= hydra/indexing/index_search.cpp
sources+= hydra/indexing/index_spinhalf.cpp
# sources+= hydra/indexing/index_symmetrized.cpp
sources+= hydra/indexing/index_table.cpp
sources+= hydra/indexing/lintable.cpp


sources+= hydra/models/heisenbergmodel.cpp
sources+= hydra/models/hubbardmodel.cpp
sources+= hydra/models/hubbardmodelmpi.cpp
sources+= hydra/models/hubbardmodeldetail.cpp
sources+= hydra/models/tjmodel.cpp
sources+= hydra/models/tjmodelmpi.cpp

# sources+= hydra/models/spinlessfermions.cpp


sources+= hydra/operators/bond.cpp
sources+= hydra/operators/bondlist.cpp
sources+= hydra/operators/couplings.cpp

sources+= hydra/parameters/palm_exception.cpp
sources+= hydra/parameters/parameter_value.cpp
sources+= hydra/parameters/parameters.cpp
sources+= hydra/parameters/parser.cpp


# sources+= hydra/symmetries/charactertable.cpp
# sources+= hydra/symmetries/spacegroup.cpp
# sources+= hydra/symmetries/symmetrydetail.cpp
# sources+= hydra/thermodynamics/thermodynamics_detail.cpp
# sources+= hydra/thermodynamics/thermodynamics_exact.cpp
# sources+= hydra/thermodynamics/thermodynamics_tpq.cpp


testsources+= test/tests.cpp
testsources+= test/test_qns.cpp
testsources+= test/test_states.cpp

testsources+= test/test_bondlist.cpp
# testsources+= test/test_charactertable.cpp
testsources+= test/test_combinatorics.cpp
testsources+= test/test_hubbard.cpp
# testsources+= test/test_indexhubbard.cpp
# testsources+= test/test_indextable.cpp
# testsources+= test/test_lintable.cpp

testsources+= test/test_spacegroup.cpp
testsources+= test/test_spinhalf.cpp
testsources+= test/test_tjmodel.cpp
testsources+= test/test_tjmodelmpi.cpp

testsources+= test/test_hubbardmodel.cpp
testsources+= test/test_hubbardmodelmpi.cpp
# testsources+= test/test_up_down_hole.cpp

# appsources+= hydra/applications/hubbardopticalftlm/hubbardopticalftlm.cpp
# appsources+= hydra/applications/hubbarddynamics/hubbarddynamics.cpp
# appsources+= hydra/applications/hubbarddynamicsmpi/hubbarddynamicsmpi.cpp


appsources+= hydra/applications/hubbarded/hubbarded.cpp
appsources+= hydra/applications/hubbardfulled/hubbardfulled.cpp
appsources+= hydra/applications/tjed/tjed.cpp
appsources+= hydra/applications/tjfulled/tjfulled.cpp

