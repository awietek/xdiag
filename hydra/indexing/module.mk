local_dir  := hydra/indexing
local_src  := $(addprefix $(local_dir)/,indextable.cpp indexsearch.cpp indexspinhalf.cpp indexsymmetrized.cpp)
sources    += $(local_src)
