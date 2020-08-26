# cl compiler flags

CC := cl
CLINK := cl

RM := rm -f
RMDIR := rm -rf
MV := mv
MKDIR := mkdir -p
CP := cp
PERL := perl
PYTHON := python
DIFF := diff

ifneq ($(debugmode),0)
  # debug flags
  OPT := /Od /Z7 /EHsc
  COPT := /Od /Z7 /EHsc /bigobj
else
  # optimized flags
  OPT := /Ox /EHsc /MT /bigobj
  COPT := /Ox /EHsc /MT /bigobj /wd4661
endif

# Object-naming switches
NAMEOBJ     := /Fo
NAMECOBJ     := /Fo
OBJEXT      := .obj

# Linking switches
NOLNK       := -c
LINK        :=
LINKFLAGS   := -m$(BITS)
NAMEEXE := /Fe

# archiving
AR          := lib
LIBEXT      := a
# no space after colon
ARFLAGS     := /OUT:
AREXTRACTFLAGS := x
# flag for supplying objects list
AROFLAG     :=

# making objects
CFLAGS := $(COPT) $(NOLNK)

# making executables
CLFLAGS   := $(LINKFLAGS) $(COPT)
OTHER_SYS_LIBS := user32.lib advapi32.lib winspool.lib /link /FORCE
