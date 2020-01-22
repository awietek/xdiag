local_dir  := test
local_test_src  := $(addprefix $(local_dir)/,test_combinatorics.cpp test_spinhalf.cpp test_indextable.cpp test_hubbard.cpp test_indexhubbard.cpp test_bondlist.cpp test_spacegroup.cpp test_charactertable.cpp)
test_sources    += $(local_test_src)
tests   += $(addprefix $(local_dir)/,tests)
