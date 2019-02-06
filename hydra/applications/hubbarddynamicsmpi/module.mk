local_dir  := hydra/applications/hubbarddynamicsmpi
local_src  := $(addprefix $(local_dir)/,)
sources    += $(local_src)
mpiprograms   += $(addprefix $(local_dir)/,hubbarddynamicsmpi)
