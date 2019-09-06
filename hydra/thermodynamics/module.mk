local_dir  := hydra/thermodynamics
local_src  := $(addprefix $(local_dir)/,thermodynamics_exact.cpp thermodynamics_detail.cpp thermodynamics_tpq.cpp)
sources    += $(local_src)

# thermodynamics_tpq.cpp 
