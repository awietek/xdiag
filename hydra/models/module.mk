local_dir  := hydra/models
local_src  := $(addprefix $(local_dir)/, hubbardmodeldetail.cpp hubbardmodel.cpp heisenbergmodel.cpp spinlessfermions.cpp hubbardmodelmpi.cpp tjmodelmpi.cpp)
sources    += $(local_src)
