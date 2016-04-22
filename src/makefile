# CFACV -- a simple collective-variables module for NAMD
#
# Cameron F Abrams 2016
# Drexel University, Philadelphia, Pennsylvania
#
# cameronabrams@gmail.com
#
# makefile for cfacv.so, libgenericdataspace.so
# forcesFromLog and Reconstruct (for Single-Sweep)

CC=gcc

GSL = -lgsl -lgslcblas
LIBS += -lm

CFLAGS += -O3

OBJ += measurements.o
OBJ += cfacv.o
OBJ += centers.o

SSOBJ += basisrep.o

all: cfacv.so genericdataspace.so ForcesFromLog Reconstruct

cfacv.so:  $(OBJ) cfacv_wrap.o
	$(CC) $(CCFLAGS) -shared -o $@ $^ $(LIBS)
	cp $@ ../lib/

genericdataspace.so:  genericdataspace.o genericdataspace_wrap.o
	$(CC) $(CCFLAGS) -shared -o $@ $^ $(LIBS)
	cp $@ ../lib/

%.o: %.c
	$(CC) -c -fpic $(CCFLAGS) $< $(LIBS)

%_wrap.c: %.i %.c %.h
	swig -tcl8 $<

ForcesFromLog: ForcesFromLog.c
	$(CC) $(CCFLAGS) -o $@ $^ $(LIBS)
	cp $@ ../bin/

Reconstruct:  $(SSOBJ) Reconstruct.o
	$(CC) $(CCFLAGS) -o $@ $^ $(LIBS) $(GSL)
	cp $@ ../bin

measurements.o : measurements.h
centers.o : centers.h
cfacv.h : measurements.h centers.h wrapcoords.h 
cfacv.o : cfacv.h
Reconstruct.c: basisrep.h


.PHONY: clean force

clean:
	rm -f *.o *.so *_wrap.c

