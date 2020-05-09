arch = osx

include options.mk
include sources.mk

objects = $(subst .cpp,.o,$(sources))
depends = $(subst .cpp,.d,$(sources))

appobjects = $(subst .cpp,.o,$(appsources))
appdepends = $(subst .cpp,.d,$(appsources))
appbinaries = $(subst .cpp,,$(appsources))

testobjects = $(subst .cpp,.o,$(testsources))
testdepends = $(subst .cpp,.d,$(testsources))
depflags = -MT $@ -MMD -MP -MF $*.d

.PHONY: all 
all:  $(objects) lib

.PHONY: test 
test:  $(objects) $(testobjects) lib 
	$(cc) $(ccopt) $(ccarch) $(depflags) $(libraries) $(testobjects) $(objects) -o test/tests 

.PHONY: apps
apps: $(objects) $(appobjects) $(appbinaries) lib

.PHONY: lib
lib: $(objects)
	ar rcs lib/libhydra.a $(objects)

$(appbinaries):
	$(cc) $(ccopt) $(ccarch) $(depflags) -Llib -lhydra $(libraries) $@.o -o bin/$(notdir $@)

$(depends):
include $(depends)

$(appdepends):
include $(appdepends)

$(testdepends):
include $(testdepends)


.PHONY: clean
clean:
	$(RM) -r $(objects) $(appobjects) $(testobjects) $(depends) $(appdepends) $(testdepends)

.PHONY: rebuild
rebuild: clean all lib

%.o: %.cpp %.d
%.o: %.cpp 
	$(cc) $(ccopt) $(ccarch) $(depflags) -c $< -o $@ $(includes)
