#
#  Let's compile the TJFC class!
#
#ifndef SKOFL_ROOT
#  SKOFL_ROOT = ../..
#endif

#include $(SKOFL_ROOT)/config.gmk
include ./config.gmk

# Different operating systems handle shared libraries diffently
# Also, encode the path in the OSX libraries so you dont have to set an enviorment
# variable.
ifeq ($(strip $(SYSTEM)),Darwin)
 SHAREDFLAG = -dynamiclib -undefined dynamic_lookup -install_name $(LIBDIR)/lib$*.so
else
 SHAREDFLAG = -shared
endif

ROOTCFLAGS   := $(shell $(ROOTSYS)/bin/root-config --cflags)  
ROOTLIBS     := $(shell $(ROOTSYS)/bin/root-config --glibs) -lTreePlayer 
ROOTLIBDIR   := $(shell $(ROOTSYS)/bin/root-config --libdir) 


# override to get the rootcint we want!!!
ROOTCINT = $(ROOTSYS)/bin/rootcint

#
#  Objects
#

SKROOT=  TJFC.o

SOBJS =  TJFC.so

INCFILES = TJFC.h

OBJS = $(SKROOT)

LIBNAME = TJFC

LIBDIR = ./lib/

#
#  Rules for building library 
#

.PHONY:  lib$(LIBNAME).a $(LIBDIR)lib$(LIBNAME).a

lib$(LIBNAME).a : $(OBJS)
	$(RM) $@
	$(AR) $@ $(OBJS) 
	$(RANLIB) $@

$(LIBDIR)lib$(LIBNAME).a : lib$(LIBNAME).a
	$(RM) $@
	$(INSTALL_LIB) $< $@

.cc.so: 
	$(ROOTCINT) -f $*Dict.cc -c -I. $*.h $*LinkDef.h
	$(CXX) $(CXXFLAGS) -c $*Dict.cc
	$(CXX) $(CXXFLAGS) -c $*.cc
#	$(CXX) $(CXXFLAGS) $(SHAREDFLAG) $*Dict.o $*.o -o lib$*.so
	$(CXX) $(CXXFLAGS) $(SHAREDFLAG) $*Dict.o $*.o -o $*.so
#	mv lib$*.so $(LIBDIR)
	rm -f $*Dict.cc $*Dict.h $*Dict.o $*.o
#
#  Targets 
#

.PHONY:  clean setup includes install.includes depend lib install.lib exec install.exec


emptyrule:: lib

clean::
	$(RM) *.o *~ *.a *Dict.cc *Dict.h *.exe *.so core* AutoDict*

setup::

includes:: $(INCFILES) 

install.includes:: $(INCFILES) 
	$(INSTALL_INC) -t $(INCDIR) $(INCFILES) 

depend::

lib:: lib$(LIBNAME).a $(SOBJS)

install.lib:: $(LIBDIR)lib$(LIBNAME).a 

exec::

install.exec:: 
