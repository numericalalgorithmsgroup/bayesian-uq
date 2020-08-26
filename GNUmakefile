# Makefile to build project within the current directory.
# It is expected that external dependencies (e.g., NAG) are provided
# as a library/modules to link to.
#
# It is assumed that the build directory is not necessarily the checkout
# directory, typical setting might be:
#   [DIR]/project/        ... repository checkout with this GNUmakefile
#   [DIR]/project_build   ... build area with various builds
# The build is triggered from the build area, i.e., [DIR]/project_build
# is the current directory and make is invoked as
#   make -f ../project/GNUmakefile [option1=xx option2=yy ...] [target1 target2 ...]
# This behaviour can be modified if $(PROJECTROOT) (default ../$(PROJECT))
# and $(BUILDROOTDIR) (default .) are given.
#
# The main options are (defaults in brackets):
#   BITS (64), debugmode (0), CCOMPILER (icl), STORERES (1),
#   INC holds an include file to modify the defaults and provide location
#       of the libraries to link to (guessed by system)
# The typical targets are:
#   cleanall, lib, test1.exe, test1.r (where test1 sits in project/tests/)
#
# Call target help to see full details.

##################################################
# set defaults, they can be always overridden as a command line argument,
# e.g., BITS=32 debugmode=1...

# Should the makefile introduce itself?
PRINTMYNAME := 0

SHELL=/bin/sh

ifeq ($(PRINTMYNAME),1)
  $(info This is $(lastword $(MAKEFILE_LIST)))
  $(info Current location: $(shell pwd))
endif

# set project name to be name of current working directory
#PROJECT := $(notdir $(shell pwd))
PROJECT := bayesian-uq
MAKEFILEPATH := $(abspath $(lastword $(MAKEFILE_LIST)))
PROJECTPATH := $(dir $(MAKEFILEPATH))

# PROJECTROOT should point at the top-level "engine" directory
# of a checked-out copy of the the engine.
# Set PROJECTROOT to the root of the engine build directory,
# i.e. a checked-out copy of the repository.
#
# MUST be absolute or relative to the build directory!!!
ifndef PROJECTROOT
  PROJECTROOT := ../$(PROJECT)
endif
#
# Need absolute path? But careful in Cygwin
#PROJECTROOT  := $(shell (cd $(PROJECTROOT); pwd))

# PROJ is another name for the project and it serves only for BUILD_DIR_NAME,
# typically it would be set via the 'shortcut' ./makefile
# This gets useful when building from various checkouts or projects
# in the same build directory.
# use the same as PROJECT unless set otherwise
ifndef PROJ
  PROJ := $(PROJECT)
endif

# Unless specified by the including GNUmakefile, or overridden on
# the command line, build in scratch directory. All
# objects, executables, results etc. end up in a tree under here.
# All objects get built under ../scratch/BUILDROOTDIR
#BUILDROOTDIR := ../scratch
BUILDROOTDIR := .
# Need absolute path?
#BUILDROOTDIR := $(shell (cd $(BUILDROOTDIR); pwd))

# Specify VERBOSE on the command line to echo every hidden command
ifndef VERBOSE
  QUIET := @
endif

# Default to assuming a 64-bit build.
ifndef BITS
  BITS := 64
endif

# No OMP alternative code by default
# (Same routines expected as in the normal build, i.e., no new names for now.)
ifndef OMP
  OMP := 0
endif

# optimized build by default
ifndef
  debugmode := 0
endif

# icl by default
ifndef
	CCOMPILER := gcc
endif

# don't store results in result file, rather print them
ifndef STORERES
  STORERES := 0
endif


##################################################
# find out system/platform/user for convenience
# of INC files and similar decision

# find out operating system
uname_os := $(shell uname -s)

ifeq ($(findstring CYGWIN,$(uname_os)),CYGWIN)
  platos := cygwin
else ifeq ($(uname_os),Linux)
  platos := linux
else
  platos := $(uname_os)
endif

# find out machine name (e.g., olney.nag.co.uk, jf-horley, ...)
machine := $(shell uname -n)

# find out user
user := $(shell whoami)
ifeq ($(user),Philip)
  user := philipm
else ifeq ($(user),sy834292)
  user := philipm
else ifeq ($(user),philip.maybank)
  user := philipm
endif

$(info INFO machine $(machine) with OS $(platos) as user $(user))
$(info INFO compiler settings $(CCOMPILER))

##################################################
# Load settings of the compiler

COMPILERINC := $(PROJECTROOT)/makefiles/compilers/GNUmakefile.$(CCOMPILER).inc

include $(COMPILERINC)

$(COMPILERINC) :
	$(error Is settings of CCOMPILER as expected? $(COMPILERINC) does not exist)

##################################################
# settings for scalar (AD) types and seed for Stan RNG

# default is to use dco types
scalar_type := dco
dco_version := dcl6i34ngl

stan_seed := 1736
stan_num_samples := 1000
finite_difference_type := cfd

##################################################
# Load personal settings (location of the libraries to link to etc.)

ifndef USERINC
  ifndef INC
    # guess based on the username
    INC := $(user)
  endif
  USERINC := $(PROJECTROOT)/makefiles/GNUmakefile.$(INC).inc
endif

include $(USERINC)
STAN_USER_HEADER := stan/spec-fun.hpp
USER_CINC := $(USER_CINC) -I$(EIGENDIR) -I$(DCOINCDIR) -I$(FFTWINCDIR)
USER_LIBS := $(USER_LIBS) $(FFTWLIB) $(DCOLIB)

$(USERINC) :
	$(error Is settings of INC as expected? $(USERINC) does not exist)

##################################################
# Set user defines e.g., USE_DCO_TYPES

ifeq ($(scalar_type),basic)
  USER_DEF := $(USER_DEF) -DUSE_BASIC_TYPES
endif

ifeq ($(finite_difference_type),cfd)
  USER_DEF := $(USER_DEF) -DUSE_CFD
endif

ifeq ($(finite_difference_type),ffd)
  USER_DEF := $(USER_DEF) -DUSE_FFD
endif

ifeq ($(scalar_type),dco)
#  USER_DEF := $(USER_DEF) -DUSE_DCO_TYPES -DDCO_DISABLE_AUTO_WARNING -DDCO_DISABLE_AVX2_WARNING -march=native
  USER_DEF := $(USER_DEF) -DUSE_DCO_TYPES -DDCO_DISABLE_AUTO_WARNING -DDCO_DISABLE_AVX2_WARNING -DDCO_MAX_ALLOCATION=1e6 -DDCO_MEM_RATIO=0.5
  USER_LIBS := $(USER_LIBS) $(DCO_LIB)
endif


# ======================================================================
# Definition of macros
# ======================================================================

####################
# Macro to link to object
# $(1) is the objects needed for linking (if any)
# $(2) additional stuff (e.g., libraries)
# $(3) is the output exe-name
# Usage: $(call clink,objs-in,libs,exe-out)
define clink
    $(QUIET) echo "Linking: $(CLINK) $(CLFLAGS) $(1) $(NAMEEXE)$(3) $(CLFLAGS_POST) $(2) $(USER_LIBS) $(OTHER_SYS_LIBS)"
	$(QUIET) $(CLINK) $(CLFLAGS) $(1) $(NAMEEXE)$(3) $(CLFLAGS_POST) $(2) $(USER_LIBS) $(OTHER_SYS_LIBS)
endef


####################
# Macro to list objects in the variable either with full
# path or just names
# $(1) things to print
# Usage: $(call list,things-to-print)
ifeq ($(LISTINFULL),1)
  list = $(1)
else
  list = $(notdir $(1))
endif

####################
# Macro to run an executable; it tests if the datafile exists
# and if so it redirects it to stdin, the results may be redirected
# to the results file
# $(1) is the executable
# $(2) is an data file (tested for existance)
# $(3) is a redirection to the result file (optional)
# Usage: $(call runex,program.exe,program.d,> program.r)
define runex
     $(QUIET) ( if [ -f $(2) ]; then        \
       echo "Running $(1) < $(2) $(3)";     \
       time $(1) < $(2) $(3);                    \
     else                                   \
       echo "Running $(1) $(3)";            \
       time $(1) $(3);                           \
     fi )
endef

# ======================================================================
# Definition of Directories
# ======================================================================

# Where to really build?
ifneq ($(debugmode),0)
  EXTRADIRNAME := $(EXTRADIRNAME)_DEBUG$(debugmode)
endif
ifeq ($(OMP),1)
  EXTRADIRNAME := $(EXTRADIRNAME)_OMP
endif

EXTRADIRNAME := $(EXTRADIRNAME)

BUILD_DIR_NAME := BUILD_$(PROJ)_$(platos)_$(notdir $(basename $(CCOMPILER)))_$(BITS)$(EXTRADIRNAME)
BUILDDIR := $(BUILDROOTDIR)/$(BUILD_DIR_NAME)

$(info INFO building in $(BUILDDIR))

##### SOURCES DIRECTORIES #####

# Source dir
BASESRCDIR := $(PROJECTROOT)/source
TESTSRCDIR := $(PROJECTROOT)/examples
STANSRCDIR := $(PROJECTROOT)/stan
INCDIR := $(PROJECTROOT)/include

BASEDEPSDIR := $(PROJECTROOT)/deps/source
TESTDEPSDIR := $(PROJECTROOT)/deps/examples


# tests extra dirs
TESTBASERESULTS := $(TESTSRCDIR)/baseresults_$(filter basic dco, $(scalar_type))
TESTDATA := $(TESTSRCDIR)/data

##### BUILD DIRECTORIES #####

ifeq ($(scalar_type),basic)
  # A directory to hold compiled common objects
  OBJDIR       := $(BUILDDIR)/OBJECTS_$(scalar_type)
  # A directory to hold the libraries
  LIB_DIR := $(BUILDDIR)/LIBS_$(scalar_type)
  # compiled tests dirs
  TEST_BUILD_DIR := $(BUILDDIR)/examples_$(scalar_type)
  # compiled Stan executable
  STANEXEDIR := ./stan/build_$(scalar_type)
endif

ifeq ($(scalar_type),dco)
  # A directory to hold compiled common objects
  OBJDIR       := $(BUILDDIR)/OBJECTS_$(scalar_type)_$(dco_version)
  # A directory to hold the libraries
  LIB_DIR := $(BUILDDIR)/LIBS_$(scalar_type)_$(dco_version)
  # compiled tests dirs
  TEST_BUILD_DIR := $(BUILDDIR)/examples_$(scalar_type)_$(dco_version)
  # compiled Stan executable
  STANEXEDIR := ./stan/build_$(scalar_type)_$(dco_version)
endif

TESTOBJDIR := $(TEST_BUILD_DIR)/OBJECTS
TESTEXEDIR := $(TEST_BUILD_DIR)/EXECUTABLES
TESTRESDIR := $(TEST_BUILD_DIR)/results
TESTDIFFDIR := $(TEST_BUILD_DIR)/diffs


# ======================================================================
# Special files named and other derivatives
# ======================================================================

# main project library name
LIBMAIN    := lib$(PROJECT).$(LIBEXT)
FULLLIBMAIN:= $(LIB_DIR)/$(LIBMAIN)


# ======================================================================
# Gathering files
# ======================================================================

.SUFFIXES:
.SUFFIXES: .exe $(OBJEXT) .cpp .depends

# The base source files
BASE_CSOURCES := $(wildcard $(BASESRCDIR)/*.cpp)
BASE_CSOURCES_IGN := $(patsubst %.ignore,%,$(wildcard $(BASESRCDIR)/*.cpp.ignore))
BASE_CSOURCES := $(filter-out $(BASE_CSOURCES_IGN),$(BASE_CSOURCES))
BASE_CSOURCES := $(sort $(BASE_CSOURCES))

# dependency files for source directory
BASE_CDEPS := $(patsubst %.cpp,%.depends,$(notdir $(BASE_CSOURCES)))
BASE_CDEPS := $(patsubst %,$(BASEDEPSDIR)/%,$(BASE_CDEPS))
BASE_CDEPS := $(sort $(BASE_CDEPS))

# Test sources
TEST_CSOURCES := $(wildcard $(TESTSRCDIR)/*.cpp)

# Test sources
STAN_SOURCES := $(wildcard $(STANSRCDIR)/*.stan)

# dependency files for tests directory
TEST_CDEPS := $(patsubst %.cpp,%.depends,$(notdir $(TEST_CSOURCES)))
TEST_CDEPS := $(patsubst %,$(TESTDEPSDIR)/%,$(TEST_CDEPS))
TEST_CDEPS := $(sort $(TEST_CDEPS))

# The base files objects
BASE_COBJECTS := $(patsubst %.cpp,%$(OBJEXT),$(notdir $(BASE_CSOURCES)))
BASE_COBJECTS := $(patsubst %,$(OBJDIR)/%,$(BASE_COBJECTS))
BASE_COBJECTS := $(sort $(BASE_COBJECTS))

# Test objects
TEST_COBJECTS := $(sort $(addprefix $(TESTOBJDIR)/,$(notdir $(TEST_CSOURCES:.cpp=$(OBJEXT)))))

# Test Executables
TEST_EXES := $(sort $(addprefix $(TESTEXEDIR)/,$(notdir $(TEST_CSOURCES:.cpp=.exe))))

# Test results
TEST_RES := $(sort $(addprefix $(TESTRESDIR)/,$(notdir $(TEST_EXES:.exe=.r))))

# Test diffs
TEST_DIFFS := $(sort $(addprefix $(TESTDIFFDIR)/,$(notdir $(TEST_EXES:.exe=.x))))

# Stan Executables
STAN_EXES := $(sort $(addprefix $(STANEXEDIR)/,$(notdir $(STAN_SOURCES:.stan=.exe))))

# All objects to put into the library
ALLLIB_OBJECTS := $(BASE_COBJECTS)

# ======================================================================
# Compilation rules
# ======================================================================

# Phony targets
.PHONY: stan first module modules clean cleanall cleanexe cleanres cleandocs cleandepends cleancoverage all lib libs libengine help listvars listmakevars printvars docs tags ctags results coverage $(notdir $(TEST_RES)) $(notdir $(TEST_DIFFS)) $(notdir $(TEST_EXES)) $(notdir $(STAN_EXES)) $(notdir $(ALLLIB_OBJECTS)) $(notdir $(TEST_OBJECTS) $(TEST_COBJECTS))


ifeq ($(OMP),1)
  VPATH = $(PROJECTROOT)/source-omp : $(MODSRCDIR) : $(BASESRCDIR) : $(TESTSRCDIR) : $(STANSRCDIR)
else
  VPATH = $(MODSRCDIR) : $(BASESRCDIR) : $(TESTSRCDIR) : $(STANSRCDIR)
endif

# if no arguments are supplied to make just build shared library from code in source directory
libs : $(FULLLIBMAIN)


# Dependency files
# list all the source files as prerequisites to the file containing depndency information (.depends)
# if any file changes that would cause the object file to appear out of date (e.g. adding an additional #include in source),
# the .depends file should also be rebuilt

# regex comments
# \([^:]*\) # matches string up to first colon and writes result to \1
# :         # matches colon
# [^\/]*    # matches string up to next forward slash
# \/        # matches forward slash
# \(.*\)    # matches rest of line and writes result to \2
#%.depends : %.cpp

$(BASE_CDEPS) : $(BASEDEPSDIR)/%.depends : $(BASESRCDIR)/%.cpp
	$(QUIET) $(MKDIR) $(BASEDEPSDIR)
	$(CC) $(CDEPOPT) $(CFLAGS) $(USER_DEF) $(USER_CINC) -I$(INCDIR) $< | grep -i '$(INCDIR)\|\.o:' | sed 's#\([^:]*\):[^\/]*\/[^\/]*\/\(.*\)#$$(OBJDIR)/\1: \2#' > $@

$(TEST_CDEPS) : $(TESTDEPSDIR)/%.depends : $(TESTSRCDIR)/%.cpp
	$(QUIET) $(MKDIR) $(TESTDEPSDIR)
	$(CC) $(CDEPOPT) $(CFLAGS) $(USER_DEF) $(USER_CINC) -I$(INCDIR) $< | grep -i '$(INCDIR)\|\.o:' | sed 's#\([^:]*\):[^\/]*\/[^\/]*\/\(.*\)#$$(TESTOBJDIR)/\1: \2#' > $@

# include makefiles containing rules that cause object files to be rebuilt when header is out of date
include $(TEST_CDEPS)
include $(BASE_CDEPS)
include custom.depends

# Objects

$(BASE_COBJECTS): $(OBJDIR)/%$(OBJEXT): $(BASESRCDIR)/%.cpp $(BASEDEPSDIR)/%.depends
	$(QUIET) $(MKDIR) $(OBJDIR)
	$(CC) $(CFLAGS) $(USER_DEF) $(USER_CINC) -I$(INCDIR) $< $(NAMECOBJ)$@

$(FULLLIBMAIN) : $(ALLLIB_OBJECTS)
	$(QUIET) $(MKDIR) $(LIB_DIR)
	$(AR) $(ARFLAGS)$@ $(AROFLAG)$^

$(TEST_COBJECTS): $(TESTOBJDIR)/%$(OBJEXT): $(TESTSRCDIR)/%.cpp $(TESTDEPSDIR)/%.depends
	$(QUIET) $(MKDIR) $(TESTOBJDIR)
	$(CC) $(CFLAGS) $(USER_DEF) $(USER_CINC) -I$(INCDIR) $< $(NAMECOBJ)$@

# Exes

# All our examples need the library created from files in source directory
$(TEST_EXES): $(FULLLIBMAIN)

$(TEST_EXES): $(TESTEXEDIR)/%.exe: $(TESTOBJDIR)/%$(OBJEXT)
	$(QUIET) $(MKDIR) $(TESTEXEDIR)
	$(QUIET) $(call clink,$<,$(FULLLIBMAIN),$@)

# Results, depending on $(STORERES) either store them to a result
# file and create a diff or just run them with the output to stdout
ifeq ($(STORERES),1)
# store results as usual and make diff

$(TEST_RES): $(TESTRESDIR)/%.r : $(TESTEXEDIR)/%.exe
	$(QUIET) $(MKDIR) $(TESTRESDIR)
	$(QUIET) $(MKDIR) $(TESTDIFFDIR)
	$(QUIET) $(call runex,$<,$(TESTDATA)/$*.d,> $@)
	-$(DIFF) $@ $(TESTBASERESULTS)/$*.r > $(TESTDIFFDIR)/$*.x  2>&1

# diffs - done with results
$(TEST_DIFFS) : $(TESTDIFFDIR)/%.x : $(TESTRESDIR)/%.r

else
# don't store results
.PHONY: $(TEST_RES)

$(TEST_RES): $(TESTRESDIR)/%.r : $(TESTEXEDIR)/%.exe
	$(QUIET) $(call runex,$<,$(TESTDATA)/$*.d, )

endif

$(STAN_EXES): $(FULLLIBMAIN)

$(STAN_EXES): $(STANEXEDIR)/%.exe : $(STANSRCDIR)/%.stan
	$(RM) $(PROJECTROOT)/$(STANSRCDIR)/$(basename $@)
	cd ${STANDIR} && \
		EIGEN=$(EIGENDIR) CXXFLAGS="-I$(DCOINCDIR) -I$(PROJECTPATH)/include $(USER_DEF)" \
		LDLIBS="$(PROJECTPATH)/$(FULLLIBMAIN) $(DCOLIB)" \
		$(MAKE) $(PROJECTPATH)/$(STANSRCDIR)/$(notdir $(basename $@)) STANCFLAGS="--allow_undefined" \
			USER_HEADER=$(PROJECTPATH)/$(STAN_USER_HEADER)
	$(QUIET) $(MKDIR) $(STANEXEDIR)
	$(MV) $(PROJECTPATH)/$(STANSRCDIR)/$(notdir $(basename $@)) $(PROJECTPATH)/$@

# cleaning

clean:
	$(RMDIR) $(BUILDDIR)
	$(RMDIR) $(STANEXEDIR)

cleanexe:
	$(RMDIR) $(TESTEXEDIR)


cleanres:
	$(RMDIR) $(TESTRESDIR)
	$(RMDIR) $(TESTDIFFDIR)

cleandepends:
	$(RMDIR) $(BASEDEPSDIR)
	$(RMDIR) $(TESTDEPSDIR)

cleanall: clean

# shortcuts
all: $(FULLLIBMAIN) $(TEST_EXES)

STANRESDIR := stan/results_$(scalar_type)

clean-stan:
	$(RMDIR) $(STANRESDIR)
	$(RMDIR) $(STANEXEDIR)

results: $(TEST_RES)

# hard-coded recipes for Stan patch, sampling and results analysis
stan_patch:
	patch $(STANDIR)/stan/lib/stan_math/stan/math/rev/fun/Eigen_NumTraits.hpp stan_math_eigen.patch


$(STANRESDIR)/spectral-inference_samples.csv: $(STAN_EXES)

$(STANRESDIR)/spectral-inference_samples.csv: $(STANEXEDIR)/spectral-inference.exe
	$(QUIET) $(MKDIR) $(STANRESDIR)
	$(STANEXEDIR)/spectral-inference.exe sample num_samples=$(stan_num_samples) \
		data file=stan/spectral-inference.data.R \
		random seed=$(stan_seed) \
		output file=$@

stan_tools:
	cd ${STANDIR} && EIGEN=$(EIGENDIR) \
		$(MAKE) build

stansummary_nuts : stan_tools $(STANRESDIR)/spectral-inference_samples.csv

stansummary_nuts :
	$(RM) $(STANRESDIR)/spectral-inference_summary.csv
	${STANDIR}/bin/stansummary --csv_file=$(STANRESDIR)/spectral-inference_summary.csv $(STANRESDIR)/spectral-inference_samples.csv
	python tabulate_stansummary.py $(STANRESDIR)/spectral-inference_summary.csv

examples/samples/harmonic-oscillator/output_$(scalar_type).csv: $(TEST_EXES)
	$(QUIET) $(MKDIR) examples/samples/harmonic-oscillator
	time $(TEST_BUILD_DIR)/EXECUTABLES/smMALA-harmonic-oscillator.exe

stansummary_smMALA : stan_tools examples/samples/harmonic-oscillator/output_$(scalar_type).csv
	$(RM) examples/samples/harmonic-oscillator/summary_$(scalar_type).csv
	${STANDIR}/bin/stansummary --csv_file=examples/samples/harmonic-oscillator/summary_$(scalar_type).csv \
                                        examples/samples/harmonic-oscillator/output_$(scalar_type).csv
	python tabulate_stansummary.py examples/samples/harmonic-oscillator/summary_$(scalar_type).csv

#	Rule to allow executable files to be specified without directory
$(notdir $(TEST_EXES)) : %.exe : $(TESTEXEDIR)/%.exe

# Shortcut for results & diffs
$(notdir $(TEST_RES)) : %.r : $(TESTRESDIR)/%.r
$(notdir $(TEST_DIFFS)) : %.x : $(TESTDIFFDIR)/%.x

# Shortcut for objects
$(notdir $(TEST_OBJECTS) $(TEST_COBJECTS)) : %$(OBJEXT) : $(TESTOBJDIR)/%$(OBJEXT)

# Shortcut for Stan executable and results
$(notdir $(STAN_EXES)) : %.exe : $(STANEXEDIR)/%.exe

# shortcut for dependencies
# $(notdir $(TEST_CDEPS)) : %.depends :  $(TESTDEPSDIR)/%.depends

listmakevars:
	@echo "============================================================="
	@echo "List of 'collected' make vars"
	@echo "(to display full path in the file listings, set LISTINFULL=1)"
	@echo "-------------------------------------------------------------"
	@echo "Settings..."
	@echo "  platos            = $(platos)"
	@echo "  CCOMPILER         = $(CCOMPILER)"
	@echo "  CLINK             = $(CLINK)"
	@echo "  CC                = $(CC)"
	@echo "  COPT              = $(COPT)"
	@echo "  CFLAGS            = $(CFLAGS)"
	@echo "  CLFLAGS           = $(CLFLAGS)"
	@echo "  USER_DEF          = $(USER_DEF)"
	@echo "Dependencies ..."
	@echo "  BASE_CDEPS        = $(call list,$(BASE_CDEPS))"
	@echo "  TEST_CDEPS        = $(call list,$(TEST_CDEPS))"
	@echo "Dirs (deps) ..."
	@echo "  INSTALLDIR        = $(INSTALLDIR)"
	@echo "  STANDIR           = $(STANDIR)"
	@echo "Dirs (sources) ..."
	@echo "  DEPENDSINC        = $(DEPENDSINC)"
	@echo "  TESTSRCDIR        = $(TESTSRCDIR)"
	@echo "  TESTBASERESULTS   = $(TESTBASERESULTS)"
	@echo "  TESTDATA          = $(TESTDATA)"
	@echo "Dirs (building) ..."
	@echo "  TESTOBJDIR        = $(TESTOBJDIR)"
	@echo "  TESTEXEDIR        = $(TESTEXEDIR)"
	@echo "  TESTRESDIR        = $(TESTRESDIR)"
	@echo "  TESTDIFFDIR       = $(TESTDIFFDIR)"
	@echo "Sources..."
	@echo "  VPATH             = $(VPATH)"
	@echo "  BASE_CSOURCES     = $(call list,$(BASE_CSOURCES))"
	@echo "  TEST_CSOURCES     = $(call list,$(TEST_CSOURCES))"
	@echo "  STAN_SOURCES      = $(call list,$(STAN_SOURCES))"
	@echo "Objects..."
	@echo "  ALLLIB_OBJECTS    = $(call list,$(ALLLIB_OBJECTS))"
	@echo "  SUBLIB_OBJECTS    = $(call list,$(SUBLIB_OBJECTS))"
	@echo "Libs..."
	@echo "  LIBMAIN           = $(LIBMAIN)"
	@echo "  FULLLIBMAIN       = $(FULLLIBMAIN)"
	@echo "  FULLLIBMAIN       = $(FULLLIBMAIN)"
	@echo "Exes..."
	@echo "  TEST_EXES         = $(call list,$(TEST_EXES))"
	@echo "  STAN_EXES         = $(call list,$(STAN_EXES))"
	@echo "Results..."
	@echo "  TEST_RES          = $(call list,$(TEST_RES))"
	@echo "============================================================="
