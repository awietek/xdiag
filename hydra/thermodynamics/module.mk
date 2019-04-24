local_dir  := hydra/thermodynamics
local_src  := $(addprefix $(local_dir)/,thermodynamics_exact.cpp thermodynamics_tpq.cpp thermodynamics_detail.cpp)
sources    += $(local_src)
