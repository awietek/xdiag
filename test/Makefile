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

CC         = mpicxx
CCOPT         = -O3
CCARCH        = -std=c++11 -Wall -pedantic -m64
programs     :=
mpiprograms  :=
sources      :=
libraries    := $(lapack) -lhydra -Llib -llime -L$(lime_dir)/lib -lhdf5
extra_clean  :=
CPPFLAGS     += $(addprefix -I ,$(include_dirs))

exe = tests
objects = tests.o test_hubbardmodelmpi.o test_bondlist.o test_charactertable.o test_combinatorics.o test_hubbard.o test_indexhubbard.o test_indextable.o test_spacegroup.o test_spinhalf.o test_lintable.o test_up_down_hole.o



includes = $(addprefix -I,$(module_dirs)) -I$(lila_dir) -I. -I$(clara_dir) -I$(lime_dir)


$(exe):  $(objects) 
	$(CC) $(CCOPT) $(CCARCH) -o $(exe) $(objects) $(lapack) -L../lib -lhydra

.PHONY: clean
clean:
	$(RM) -r $(objects) $(depends)

.PHONY: rebuild
rebuild: clean $(exe)

%.o: %.cpp
	$(CC) $(CCOPT) $(CCARCH) -c $< -o $@ $(includes) -I../
