ifeq ($(machine),LAPTOP-HOH0NDO7)
  ifeq ($(platos),linux)
    INSTALLDIR := /mnt/c/Users/Philip
    STANDIR  := $(INSTALLDIR)/stan-dev/cmdstan
  else ifeq ($(platos),cygwin)
    INSTALLDIR := ..
    STANDIR  := $(INSTALLDIR)/stan-dev/cmdstan
  endif
else ifeq ($(machine),pm-toab)
  ifeq ($(platos),linux)
    INSTALLDIR := /mnt/c/cygwin64/home/philip.maybank
    STANDIR  := $(INSTALLDIR)/stan-dev/cmdstan
  else ifeq ($(platos),cygwin)
    INSTALLDIR := ..
    STANDIR  := $(INSTALLDIR)/stan-dev/cmdstan
  endif
else ifeq ($(machine),olney.nag.co.uk)
  INSTALLDIR := /home/nag/olney/$(user)
  STANDIR  := $(INSTALLDIR)/stan-dev/cmdstan
else ifeq ($(machine),oban)
  INSTALLDIR := /home/nag/oban/$(user)
  STANDIR  := $(INSTALLDIR)/stan-dev/cmdstan
endif


EIGENDIR  := $(INSTALLDIR)/eigen-ad
# EIGENDIR  := $(STANDIR)/stan/lib/stan_math/lib/eigen_3.3.3
ifeq ($(platos),linux)
  FFTWINCDIR:= $(INSTALLDIR)/fftw-3.3.8/api
  FFTWLIB   := $(INSTALLDIR)/fftw-3.3.8/.libs/libfftw3.a
else ifeq ($(platos),cygwin)
  FFTWINCDIR:= $(INSTALLDIR)/fftw-3.3.5-dll64
  FFTWLIB   := $(INSTALLDIR)/fftw-3.3.5-dll64/libfftw3-3.lib
endif

STAN_USER_HEADER := stan/spec-fun.hpp

ifeq ($(dco_version), dcl6i34ngl)
  DCODIR   := /home/philipm/NAG/dcl6i34ngl
  DCOINCDIR := $(DCODIR)/include
  DCOLIB := $(DCODIR)/lib/libdcoc.a
endif

ifeq ($(dco_version), develop)
  DCODIR   := $(INSTALLDIR)/dco_cpp_dev/src
  DCOINCDIR := $(DCODIR)/
  DCOLIB :=
endif

ifeq ($(dco_version),dcw6i34nel)
  # DCODIR   := /cygdrive/c/Progra~1/NAG/DC34/dcw6i34nel
  DCODIR := "C:\Progra~1\NAG\DC34\dcw6i34nel"
  DCOINCDIR := $(DCODIR)/include
  DCOLIB := $(DCODIR)/lib/libdcoc_MT.lib
endif
