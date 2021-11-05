include options.mk
include sources.mk

rm     := rm -f
mkdir  := mkdir -p

includes = -I$(liladir)

# Additional parameters for compilation
ccstd   = -std=c++17 
ccwarn  = -Wall -pedantic -Wno-return-type-c-linkage

objects = $(subst .cpp,.o,$(sources))
depends = $(subst .cpp,.d,$(sources))

mpiobjects = $(subst .cpp,.o,$(mpisources))
mpidepends = $(subst .cpp,.d,$(mpisources))

appobjects = $(subst .cpp,.o,$(appsources))
appdepends = $(subst .cpp,.d,$(appsources))
appbinaries = $(subst .cpp,,$(appsources))

testobjects = $(subst .cpp,.o,$(testsources))
testdepends = $(subst .cpp,.d,$(testsources))

testmpiobjects = $(subst .cpp,.o,$(testmpisources))
testmpidepends = $(subst .cpp,.d,$(testmpisources))

depflags = -MT $@ -MMD -MP -MF $*.d


.PHONY: all 
all:  $(objects) lib

$(depends):
include $(depends)

$(mpidepends):
include $(mpidepends)

$(appdepends):
include $(appdepends)

$(testdepends):
include $(testdepends)

$(testmpidepends):
include $(testmpidepends)

# different building suites
.PHONY: test 
test:  $(objects) $(testobjects) lib 
	$(cc) $(ccopt) $(ccstd) $(ccwarn) $(depflags) $(libraries) -Llib -lhydra $(objects) $(testobjects) -o test/tests 

.PHONY: mpi
mpi: $(objects) $(mpiobjects) libmpi


.PHONY: testmpi 
testmpi:  $(objects) $(mpiobjects) $(testmpiobjects) libmpi 
	$(mpicc) $(ccopt) $(ccstd) $(ccwarn) $(depflags) $(libraries) -Llib -lhydra_mpi $(objects) $(mpiobjects) $(testmpiobjects) -o testmpi/tests 

.PHONY: apps
apps: $(objects) $(appobjects) $(appbinaries) lib

# Libraries
.PHONY: lib
lib: $(objects)
	ar rcs lib/libhydra.a $(objects)

.PHONY: libmpi
libmpi: $(objects) $(mpiobjects)
	ar rcs lib/libhydra_mpi.a $(objects) $(mpiobjects)

$(appbinaries):
	$(cc) $(ccopt) $(ccstd) $(ccwarn) $(depflags) $@.o -Llib -lhydra $(libraries) -o bin/$(notdir $@)



.PHONY: clean
clean:
	$(RM) -r $(objects) $(mpiobjects) $(appobjects) $(testobjects) $(testmpiobjects) $(depends) $(appdepends) $(testdepends) $(testmpidepends) $(mpidepends)

.PHONY: rebuild
rebuild: clean all lib

%.o: %.cpp %.d
%.o: %.cpp 
	$(cc) $(ccopt) $(ccstd) $(ccwarn) $(depflags) $(lilabackend) -c $< -o $@ $(includes) -I.
