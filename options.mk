ifeq ($(arch), flatiron_linux)
cc         = mpicxx
ccopt      = -O3 -mavx -DLILA_USE_MKL
ccarch     = -std=c++17 -Wall -pedantic -m64 -Wno-return-type-c-linkage
liladir    = /mnt/home/awietek/Research/Software/lila
limedir    = /mnt/home/awietek/Research/Software/lime
claradir    = /mnt/home/awietek/Research/Software/Clara
includes   = -I. -I$(liladir) -I$(limedir) -I$(claradir)/include
libraries  = -lhdf5 -lmkl_rt -L$(limedir)/lib -llime -lhdf5
endif

ifeq ($(arch), osx)
cc         = mpicxx
ccopt      = -O3 -mavx -DLILA_USE_ACCELERATE
ccarch     = -std=c++17 -Wall -pedantic -m64 -Wno-return-type-c-linkage
liladir    = /Users/awietek/Research/Software/lila
limedir    = /Users/awietek/Research/Software/lime
claradir    = /Users/awietek/Research/Software/Clara
includes   = -I. -I$(liladir) -I$(limedir) -I$(claradir)/include
libraries  = -framework Accelerate -L$(limedir)/lib -llime -lhdf5
endif

ifeq ($(arch), flatiron_gordon)
cc         = mpicxx
ccopt      = -O3 -mavx -DLILA_USE_MKL
ccarch     = -std=c++17 -Wall -pedantic -m64 -Wno-return-type-c-linkage
libraries  = -L/opt/hdf5/gnu/mvapich2_ib/lib -lhdf5 -lmkl_rt -DLILA_USE_MKL
liladir    = /home/awietek/Research/Software/lila
limedir    = /home/awietek/Research/Software/lime
claradir    = /home/awietek/Research/Software/Clara
includes   = -I. -I$(liladir) -I$(limedir) -I$(claradir)/include
endif

rm     := rm -f
mkdir  := mkdir -p
