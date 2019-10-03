local_dir  := hydra/applications/hubbardopticalftlm
local_src  := $(addprefix $(local_dir)/,)
sources    += $(local_src)
programs   += $(addprefix $(local_dir)/,hubbardopticalftlm)
