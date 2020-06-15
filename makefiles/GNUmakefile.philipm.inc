
# Must be set by user

ifeq ($(machine),LAPTOP-HOH0NDO7)
  INSTALLDIR := /mnt/c/Users/Philip
  STANDIR  := $(INSTALLDIR)/stan-dev/cmdstan
else ifeq ($(machine),olney.nag.co.uk)
  INSTALLDIR := /home/nag/olney/$(user)
else ifeq ($(machine),oban)
  INSTALLDIR := /home/nag/oban/$(user)
  STANDIR  := $(INSTALLDIR)/cmdstan-2.23.0
endif

EIGENDIR := $(INSTALLDIR)/eigen-ad
DCODIR   := $(INSTALLDIR)/dco_cpp_dev/src
FFTWDIR  := $(INSTALLDIR)/fftw-3.3.8
STAN_USER_HEADER := stan/spec-fun.hpp

PROJECTPATH := $(INSTALLDIR)/$(PROJECT)
USER_CINC := $(USER_CINC) -I$(EIGENDIR) -I$(DCODIR) -I$(FFTWDIR)/api
# USER_CINC := $(USER_CINC) -I$(EIGENDIR) -I$(FFTWDIR)/api
# USER_LIBS := $(USER_LIBS) $(FFTWDIR)/.libs/libfftw3.a $(DCODIR)/libdco.a
USER_LIBS := $(USER_LIBS) $(FFTWDIR)/.libs/libfftw3.a
