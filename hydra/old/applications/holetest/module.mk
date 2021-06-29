local_dir  := hydra/applications/holetest
local_src  := $(addprefix $(local_dir)/,)
sources    += $(local_src)
programs   += $(addprefix $(local_dir)/,tjed)
