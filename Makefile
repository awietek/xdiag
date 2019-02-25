hydra_module_dirs = hydra/hilbertspaces hydra/utils hydra/indexing test  hydra/models hydra/operators hydra/thermodynamics hydra/symmetries hydra/applications/hubbarddynamics hydra/dynamics hydra/applications/hubbarddynamicsmpi


apps=hydra/applications/heisenberged hydra/applications/hubbarded hydra/applications/spinlessfermioned hydra/applications/hubbardthermo hydra/applications/heisenbergthermo

lila_dir = /mnt/home/awietek/Research/Software/lila
clara_dir=/mnt/home/awietek/Research/Software/Clara/include

# CC            = mpicxx 
# MPICC         = mpicxx 
CC            = mpicxx -cxx=icpc
MPICC         = mpicxx -cxx=icpc
CCOPT         = -O3
CCARCH        = -std=c++11 -Wall -pedantic

lapack        = -llapack -lblas
lapack        = -mkl -DLILA_USE_MKL
programs     :=
mpiprograms  :=
sources      :=
libraries    := $(lapack)
extra_clean  :=

objects      = $(subst .cpp,.o,$(sources))
test_objects = $(subst .cpp,.o,$(test_sources))
# dependencies = $(subst .cpp,.d,$(sources))

include_dirs := lib include
CPPFLAGS     += $(addprefix -I ,$(include_dirs))
vpath %.h $(include_dirs)

RM     := rm -f
MKDIR  := mkdir -p

hydra_module_makefiles = $(addsuffix /module.mk,$(hydra_module_dirs))
hydra_build_dirs = $(subst hydra,build,$(hydra_module_dirs))
hydra_includes = $(addprefix -I,$(hydra_module_dirs)) -I$(lila_dir) -I. -I$(clara_dir)

# all:

include $(hydra_module_makefiles)

.PHONY: all
all:  $(programs)

.PHONY: mpi
mpi:  $(mpiprograms)

# .PHONY: libraries
# libraries: $(libraries)

.PHONY: lib
lib: $(objects)
	ar rcs lib/libhydra.a $(objects)
.PHONY: clean
clean:
	$(RM) $(objects) $(programs) $(extra_clean) $(test_objects)

.PHONY: rebuild
rebuild: clean all lib

%.o : %.cpp
	$(CC) $(CCOPT) $(CCARCH) -c $< -o $@ $(hydra_includes)

$(programs): $(objects)	
	$(CC) $(CCOPT) $(CCARCH) $@.cpp -o bin/$(notdir $@) $(objects) $(hydra_includes) $(libraries)  

$(mpiprograms): $(objects)	
	$(MPICC) $(CCOPT) $(CCARCH) $@.cpp -o bin/$(notdir $@) $(objects) $(hydra_includes) $(libraries)  


tests: $(test_objects) $(objects)
	$(CC) $(CCOPT) $(CCARCH) test/tests.cpp -o bin/tests $(test_objects) $(objects) $(hydra_includes) $(libraries) 
