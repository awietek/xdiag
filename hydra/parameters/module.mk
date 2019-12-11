local_dir  := hydra/parameters
local_src  := $(addprefix $(local_dir)/, parameters.cpp parameter_value.cpp parser.cpp palm_exception.cpp)
sources    += $(local_src)
