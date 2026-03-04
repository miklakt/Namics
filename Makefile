#Platform detection
UNAME_S     := $(shell uname -s)

#Compiler and Linker
# Optional parallel backends:
#   Windows: MSVC uses its own builtin thread pool (no TBB needed)
#   Linux:   install TBB (apt install libtbb-dev)
#   MacOS:   install gcc and tbb via Homebrew (brew install gcc tbb)
#            Apple clang does not support the <execution> header.
ifeq ($(UNAME_S),Darwin)
	# MacOS specific settings
    BREW_PREFIX := $(shell brew --prefix)
    GCC_VERSION := $(shell ls $(BREW_PREFIX)/bin/g++-* | sort -t- -k2 -n | tail -1)
    CC          := $(if $(GCC_VERSION),$(GCC_VERSION),g++)
else
    CC          := g++
endif

# c++20 baseline for modern language/library features across the codebase
CXX_STD     := c++20

#The Target Binary Program
TARGET      := namics

#The Directories, Source, Includes, Objects, Binary and Resources
SRCDIR      := src
INCDIR      := inc
BUILDDIR    := obj
TARGETDIR   := bin
RESDIR      := res
SRCEXT      := cpp
DEPEXT      := d
OBJEXT      := o

#flat DOUBLE

#Flags, Libraries and Includes
CFLAGS      := -Wall -O3 -ffast-math -std=$(CXX_STD) -march=native
LIB         := -lm -lpthread
INC         := -I$(SRCDIR) -Iexternal/LBFGSpp/include -I/usr/include/eigen3 -I/usr/local/include/eigen3 -I/usr/local/include -I/usr/include

# MacOS: add Homebrew paths for headers and libraries
ifeq ($(UNAME_S),Darwin)
ifdef BREW_PREFIX
    INC     += -I$(BREW_PREFIX)/include
    LIB     += -L$(BREW_PREFIX)/lib
endif
endif

ifdef PAR_STL
	CFLAGS += -DPAR_STL
	LIB += -ltbb
endif

# PAR_STL sanity checks
ifdef PAR_STL
HAS_TBB := $(shell pkg-config --exists tbb && echo yes || echo no)
ifeq ($(filter c++17 c++20 c++23 c++2a c++2b,$(CXX_STD)),)
    $(error PAR_STL requires C++17 or higher, but CXX_STD is $(CXX_STD). Set CXX_STD := c++17 in the Makefile.)
endif
ifeq ($(UNAME_S),Darwin)
ifneq ($(findstring g++,$(CC)),g++)
    $(error PAR_STL on macOS requires GCC. Apple clang does not support <execution>. Install GCC via: brew install gcc)
endif
endif
ifneq ($(UNAME_S),Windows_NT)
ifneq ($(HAS_TBB),yes)
    $(error PAR_STL requires Intel TBB but it was not found. Install with: apt install libtbb-dev (Linux) or brew install tbb (macOS))
endif
endif
endif

#Build configuration summary
ifeq ($(filter clean cleaner,$(MAKECMDGOALS)),)
$(info )
$(info Platform:    $(UNAME_S))
$(info C++ standard: $(CXX_STD))
$(info Compiler:    $(CC))
$(info Eigen:       required)
$(info LBFGSpp:     required)
ifdef PAR_STL
ifeq ($(HAS_TBB),yes)
$(info TBB:         found)
$(info Parallel:    C++17 STL (Intel TBB))
else ifeq ($(UNAME_S),Windows_NT)
$(info Parallel:    C++17 STL (MSVC thread pool))
else
$(info Parallel:    PAR_STL requested but TBB NOT found -- install with: apt install libtbb-dev (Linux) or brew install tbb (MacOS))
endif
else
$(info Parallel:    serial)
endif
$(info )
endif

#---------------------------------------------------------------------------------
#DO NOT EDIT BELOW THIS LINE
#---------------------------------------------------------------------------------

SOURCES     := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS     := $(filter $(BUILDDIR)/%, $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.$(OBJEXT))))

#Defauilt Make
all: resources $(TARGET)

#Remake
remake: cleaner all

#Copy Resources from Resources Directory to Target Directory
resources: directories
#	@cp $(RESDIR)/* $(TARGETDIR)/

#Make the Directories
directories:
	@mkdir -p $(TARGETDIR)
	@mkdir -p $(BUILDDIR)

#Clean only Objects
clean:
	@$(RM) -rf $(BUILDDIR)

#Full Clean, Objects and Binaries
cleaner: clean
	@$(RM) -rf $(TARGETDIR)

#Pull in dependency info for *existing* .o files
-include $(OBJECTS:.$(OBJEXT)=.$(DEPEXT))

#Link
$(TARGET): $(OBJECTS)
	$(CC) -o $(TARGETDIR)/$(TARGET) $^ $(LIB)

#Compile
$(BUILDDIR)/%.$(OBJEXT): $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<
	@$(CC) $(CFLAGS) $(INC) -MM $(SRCDIR)/$*.$(SRCEXT) > $(BUILDDIR)/$*.$(DEPEXT)
	@cp -f $(BUILDDIR)/$*.$(DEPEXT) $(BUILDDIR)/$*.$(DEPEXT).tmp
	@sed -e 's|.*:|$(BUILDDIR)/$*.$(OBJEXT):|' < $(BUILDDIR)/$*.$(DEPEXT).tmp > $(BUILDDIR)/$*.$(DEPEXT)
	@sed -e 's/.*://' -e 's/\\$$//' < $(BUILDDIR)/$*.$(DEPEXT).tmp | fmt -1 | sed -e 's/^ *//' -e 's/$$/:/' >> $(BUILDDIR)/$*.$(DEPEXT)
	@rm -f $(BUILDDIR)/$*.$(DEPEXT).tmp


#Non-File Targets
test-homopolymer-adsorption:
	@./tests/homopolymer_adsorption_test.sh

test-frozen-range-input-file:
	@./tests/test_frozen_range_input_file.sh

test-micelle-self-assembly:
	@./tests/test_micelle_self_assembly.sh

test-particle-in-cyl-coordinates:
	@./tests/test_particle_in_cyl_coordinates.sh

test-all:
	@./tests/test_all.sh

.PHONY: all remake clean cleaner resources test-homopolymer-adsorption test-frozen-range-input-file test-micelle-self-assembly test-particle-in-cyl-coordinates test-all
