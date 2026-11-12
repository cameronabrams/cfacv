comma := ,

# CFACV -- a simple collective-variables module for NAMD
#
# Cameron F Abrams 2016-2018, Drexel University
#
# Builds lib/cfacv.so, lib/genericdataspace.so and the three binaries
# bin/ForcesFromLog, bin/Reconstruct, bin/FEPFromForces.
#
# ---------------------------------------------------------------------------
# This makefile discovers Tcl, GSL and swig instead of hardcoding one host's
# paths. Override anything on the command line:
#
#     make all                     # autodetect; -O3, portable binary
#     make all OPT='-O3 -march=native'   # tuned to THIS cpu, not portable
#     make all TCL_INC=/path/to/dir GSL_PREFIX=/path
#     make PROFILE=1 all           # -pg, for gprof
#     make DEBUG=1 all             # -O0 -g, no optimisation
#
# FLOATING POINT: -ffast-math is deliberately NOT used. This code computes
# forces fed back into NAMD; reassociation and the loss of strict IEEE
# semantics can change a trajectory. Do not add it without re-validating.
# ---------------------------------------------------------------------------

# make predefines CC=cc, so ?= never fires; only a real override wins
ifeq ($(origin CC),default)
  CC := gcc
endif

# ---- optimisation ---------------------------------------------------------
OPT ?= -O3
WARN ?= -Wall
ifdef DEBUG
  OPT := -O0 -g
endif
ifdef PROFILE
  OPT += -pg
  LDEXTRA += -pg
endif

# ---- Tcl ------------------------------------------------------------------
# first hit wins; TCL_INC on the command line overrides discovery entirely
TCL_CANDIDATES := /usr/include /usr/include/tcl /usr/include/tcl8.6 \
                  $(CONDA_PREFIX)/include $(HOME)/miniconda3/include
TCL_INC ?= $(firstword $(foreach d,$(TCL_CANDIDATES),$(if $(wildcard $(d)/tcl.h),$(d))))
ifeq ($(strip $(TCL_INC)),)
  $(error tcl.h not found. Pass TCL_INC=/dir/containing/tcl.h)
endif

# ---- GSL (Reconstruct and FEPFromForces only) -----------------------------
GSL_CONFIG := $(shell command -v gsl-config 2>/dev/null)
ifneq ($(GSL_CONFIG),)
  GSL_CFLAGS := $(shell $(GSL_CONFIG) --cflags)
  GSL_RAW    := $(shell $(GSL_CONFIG) --libs)
  # gsl-config gives -L but no -rpath; without it the binaries fail at run time
  GSL        := $(GSL_RAW) $(patsubst -L%,-Wl$(comma)-rpath$(comma)%,$(filter -L%,$(GSL_RAW)))
else ifneq ($(strip $(GSL_PREFIX)),)
  GSL_CFLAGS := -I$(GSL_PREFIX)/include
  GSL        := -L$(GSL_PREFIX)/lib -Wl,-rpath,$(GSL_PREFIX)/lib -lgsl -lgslcblas
else ifneq ($(wildcard $(CONDA_PREFIX)/include/gsl),)
  GSL_CFLAGS := -I$(CONDA_PREFIX)/include
  GSL        := -L$(CONDA_PREFIX)/lib -Wl,-rpath,$(CONDA_PREFIX)/lib -lgsl -lgslcblas
else
  GSL_CFLAGS :=
  GSL        := -lgsl -lgslcblas
endif

# ---- swig -----------------------------------------------------------------
SWIG ?= $(shell command -v swig 2>/dev/null)
ifeq ($(strip $(SWIG)),)
  $(error swig not found. Pass SWIG=/path/to/swig, or conda install swig)
endif

# ---- assembled flags ------------------------------------------------------
CFLAGS  += $(OPT) $(WARN) -I$(TCL_INC) $(GSL_CFLAGS)
LIBS    += -lm
LDFLAGS += $(LDEXTRA)
USER_FLAGS ?=

OBJ   += measurements.o cfacv.o centers.o
SSOBJ += basisrep.o

.PHONY: all clean distclean config
all: cfacv.so genericdataspace.so ForcesFromLog Reconstruct FEPFromForces

config:
	@echo "CC       = $(CC)"
	@echo "OPT      = $(OPT)"
	@echo "TCL_INC  = $(TCL_INC)"
	@echo "SWIG     = $(SWIG)"
	@echo "GSL      = $(GSL)"
	@echo "CFLAGS   = $(CFLAGS)"

../lib ../bin:
	mkdir -p $@

cfacv.so: $(OBJ) cfacv_wrap.o | ../lib
	$(CC) $(CFLAGS) $(LDFLAGS) -shared -o $@ $^ $(LIBS)
	cp $@ ../lib/

genericdataspace.so: genericdataspace.o genericdataspace_wrap.o | ../lib
	$(CC) $(CFLAGS) $(LDFLAGS) -shared -o $@ $^ $(LIBS)
	cp $@ ../lib/

%.o: %.c
	$(CC) -c -fPIC $(CFLAGS) -o $@ $<

%_wrap.c: %.i %.c %.h
	$(SWIG) -tcl8 $<

ForcesFromLog: ForcesFromLog.c | ../bin
	$(CC) $(CFLAGS) $(USER_FLAGS) $(LDFLAGS) -o $@ $^ $(LIBS)
	cp $@ ../bin/

Reconstruct: $(SSOBJ) Reconstruct.o | ../bin
	$(CC) $(CFLAGS) $(USER_FLAGS) $(LDFLAGS) -o $@ $^ $(LIBS) $(GSL)
	cp $@ ../bin/

FEPFromForces: FEPFromForces.c | ../bin
	$(CC) $(CFLAGS) $(USER_FLAGS) $(LDFLAGS) -o $@ $^ $(LIBS) $(GSL)
	cp $@ ../bin/

# ---- header dependencies --------------------------------------------------
measurements.o : measurements.h
centers.o      : centers.h
cfacv.o        : cfacv.h measurements.h centers.h wrapcoords.h
Reconstruct.o  : basisrep.h
basisrep.o     : basisrep.h

clean:
	rm -f *.o *.so *_wrap.c ForcesFromLog Reconstruct FEPFromForces

distclean: clean
	rm -f ../lib/*.so ../bin/ForcesFromLog ../bin/Reconstruct ../bin/FEPFromForces
