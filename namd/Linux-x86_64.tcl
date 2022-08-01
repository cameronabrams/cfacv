
#TCLDIR=/Projects/namd2/tcl/tcl8.5.9-linux-x86_64
#TCLDIR=/Projects/namd2/tcl/tcl8.5.9-linux-x86_64-threaded
TCLDIR=/usr
TCLINCL=-I$(TCLDIR)/include
#TCLLIB=-L$(TCLDIR)/lib -ltcl8.5 -ldl
#TCLLIB=-L$(TCLDIR)/lib -ltcl8.5 -ldl -lpthread
TCLLIB=-L$(TCLDIR)/lib64 -ltcl8.6 -ldl
TCLFLAGS=-DNAMD_TCL
TCL=$(TCLINCL) $(TCLFLAGS)

