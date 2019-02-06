local_dir  := hydra/models
local_src  := $(addprefix $(local_dir)/,hubbardmodel.cpp heisenbergmodel.cpp spinlessfermions.cpp)
sources    += $(local_src)
