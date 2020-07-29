ifeq ($(machine),LAPTOP-HOH0NDO7)
  INSTALLDIR := /mnt/c/Users/Philip
  STANDIR  := $(INSTALLDIR)/stan-dev/cmdstan
else ifeq ($(machine),olney.nag.co.uk)
  INSTALLDIR := /home/nag/olney/$(user)
else ifeq ($(machine),oban)
  INSTALLDIR := /home/nag/oban/$(user)
  STANDIR  := $(INSTALLDIR)/cmdstan-2.23.0
endif

EIGENDIR  := $(INSTALLDIR)/eigen-ad
FFTWINCDIR:= $(INSTALLDIR)/fftw-3.3.8/api
FFTWLIB   := $(INSTALLDIR)/fftw-3.3.8/.libs/libfftw3.a
STAN_USER_HEADER := stan/spec-fun.hpp

ifeq ($(dco_version), dcl6i34ngl)
  DCODIR   := /home/philipm/NAG/dcl6i34ngl
  DCOINCDIR := $(DCODIR)/include
  DCOLIB := $(DCODIR)/lib/libdcoc.a
endif

ifeq ($(dco_version), develop)
  DCODIR   := $(INSTALLDIR)/dco_cpp_dev/src
  DCOINCDIR := $(DCODIR)
  DCOLIB :=
endif
