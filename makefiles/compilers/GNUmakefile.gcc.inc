# gcc compiler flags

CC := gcc
CLINK := g++

RM := rm -f
RMDIR := rm -rf
MV := mv
MKDIR := mkdir -p
CP := cp
PERL := perl
PYTHON := python
DIFF := diff

CDEPOPT := -MM

ifneq ($(debugmode),0)
  # debug flags
  OPT := -g -Wall -std=c++14 -O0 -fvar-tracking -fvar-tracking-assignments
  COPT := -g -Wall -std=c++14 -O0 -fvar-tracking -fvar-tracking-assignments
else
  # optimized flags
  OPT := -O3 -std=c++14 
  COPT := -O3 -std=c++14 
endif

# Object-naming switches
NAMEOBJ     := -o # Leave space at end of line
NAMECOBJ     := -o # Leave space at end of line
OBJEXT      := .o

# Linking switches
NOLNK       := -c
LINK        :=
LINKFLAGS   := -m$(BITS)
NAMEEXE := -o # Leave space at end of line

# archiving
AR          := ar
LIBEXT      := a
ARFLAGS     := r #extra space required here
AREXTRACTFLAGS := x
# flag for supplying objects list
AROFLAG     :=

# making objects
CFLAGS := $(COPT) $(NOLNK)

# making executables
CLFLAGS   := $(LINKFLAGS) $(COPT)
