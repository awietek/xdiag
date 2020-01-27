srcdir = hydra
appdir = hydra/applications
tstdir = test
arch = flatiron_linux

ifeq ($(arch), flatiron_linux)
	lapack        = -llapack -lblas
	lapack        = -lmkl_rt -DLILA_USE_MKL
	lila_dir =/mnt/home/awietek/Research/Software/lila
	lime_dir =/mnt/home/awietek/Research/Software/lime
	clara_dir=/mnt/home/awietek/Research/Software/Clara/include
endif
ifeq ($(arch), flatiron_gordon)
	lapack        = -llapack -lblas
	lapack        = -lmkl_rt -DLILA_USE_MKL
	lila_dir =/home/awietek/Research/Software/lila
	lime_dir =/home/awietek/Research/Software/lime
	clara_dir=/home/awietek/Research/Software/Clara/include
endif
ifeq ($(arch), osx)
	lapack        = -framework accelerate
	lila_dir =/Users/awietek/Research/Software/lila
	lime_dir =/Users/awietek/Research/Software/lime
	clara_dir=/Users/awietek/Research/Software/Clara/include
endif
ifeq ($(arch), hshackle_linux)
	lapack        = -llapack -lblas
	lila_dir =/home/hshackle/ed/lila
	lime_dir =/home/hshackle/ed/lime
	clara_dir=/home/hshackle/ed/clara/include
endif
modules = hilbertspaces utils indexing models operators symmetries dynamics thermodynamics parameters

apps= holetest hubbarded hubbarddynamicsmpi #hubbarddynamics #hubbardopticalftlm  # hubbardthermotpq hubbardopticalmpi # hubbarddynamicsmpi   heisenberged spinlessfermioned   hubbardthermo  heisenbergthermo  hubbardopticaltsl hubbarded 

# -g -fsanitize=address

CC         = mpicxx
CCOPT         = -O3
CCARCH        = -std=c++11 -Wall -pedantic -m64
programs     :=
mpiprograms  :=
sources      :=
libraries    := $(lapack) -lhydra -Llib -llime -L$(lime_dir)/lib -lhdf5
extra_clean  :=
CPPFLAGS     += $(addprefix -I ,$(include_dirs))
RM     := rm -f
MKDIR  := mkdir -p

# Set the variable sources
app_dirs = $(addprefix $(appdir)/,$(apps))
app_makefiles = $(addsuffix /module.mk,$(app_dirs))
module_dirs = $(addprefix $(srcdir)/,$(modules))
module_makefiles = $(addsuffix /module.mk,$(module_dirs))

include $(app_makefiles)
include $(module_makefiles)

objects = $(subst .cpp,.o,$(sources))
depends = $(subst .cpp,.d,$(sources))
depflags = -MT $@ -MMD -MP -MF $*.d


includes = $(addprefix -I,$(module_dirs)) -I$(lila_dir) -I. -I$(clara_dir) -I$(lime_dir)

.PHONY: all
all:  $(objects) lib

apps: $(apps)

$(depends):
include $(depends)

lib: $(objects)
	ar rcs lib/libhydra.a $(objects)

.PHONY: clean
clean:
	$(RM) -r $(objects) $(depends)

.PHONY: rebuild
rebuild: clean all lib

%.o: %.cpp
%.o: %.cpp %.d
	$(CC) $(CCOPT) $(CCARCH) $(depflags) -c $< -o $@ $(includes)

$(depdir): ; @mkdir -o $@

$(apps): $(appdir)/$@ lib
	$(CC) $(CCOPT) $(CCARCH) $(appdir)/$@/$@.cpp -o bin/$@ $(includes) $(libraries)

$(tests): $@ lib
	$(CC) $(CCOPT) $(CCARCH) $@.cpp -o $@ $(includes) $(libraries) $(test_objects)
