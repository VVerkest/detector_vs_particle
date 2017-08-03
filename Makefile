os = $(shell uname -s)

INCFLAGS      =  -I$(shell root-config --incdir) $(shell fastjet-config --cxxflags)  -I$(shell pythia8-config --includedir) -I$(BOOSTDIR)/include -I$(STARPICODIR) # -I$(BOOSTDIR)/include

ifeq ($(os),Linux)
CXXFLAGS      = -std=c++11
else
CXXFLAGS      = -O -std=c++11 -fPIC -pipe -Wall -Wno-deprecated-writable-strings -Wno-unused-variable -Wno-unused-private-field -Wno-gnu-static-float-init
## for debugging:
# CXXFLAGS      = -g -O0 -fPIC -pipe -Wall -Wno-deprecated-writable-strings -Wno-unused-variable -Wno-unused-private-field -Wno-gnu-static-float-init
endif

ifeq ($(os),Linux)
LDFLAGS       = -g
LDFLAGSS      = -g --shared
else
LDFLAGS       = -O -Xlinker -bind_at_load -flat_namespace
LDFLAGSS      = -flat_namespace -undefined suppress
LDFLAGSSS     = -bundle
endif

ifeq ($(os),Linux)
CXX          = g++
else
CXX          = clang
endif


ROOTLIBS      = $(shell root-config --libs)
FJLIBS        = $(shell fastjet-config --plugins=yes --libs)
PYTHIALIBS    = $(shell pythia8-config --ldflags)
LIBPATH       = -L$(FASTJETDIR)/lib -L$(STARPICODIR) $(shell root-config --libs)
LIBS          =  $(ROOTLIBS) $(FJLIBS) $(PYTHIALIBS) -lfastjet -lfastjettools -lTStarJetPico

# for cleanup
SDIR          = src
ODIR          = src/obj
BDIR          = bin


###############################################################################
################### Remake when these headers are touched #####################
###############################################################################


###############################################################################
# standard rules
$(ODIR)/%.o : $(SDIR)/%.cxx $(INCS)
	@echo
	@echo COMPILING
	$(CXX) $(CXXFLAGS) $(INCFLAGS) -c $< -o $@

$(BDIR)/%  : $(ODIR)/%.o
	@echo
	@echo LINKING
	$(CXX) $(LDFLAGS) $^ -o $@ $(LIBPATH) $(LIBS)

###############################################################################

###############################################################################
############################# Main Targets ####################################
###############################################################################
all : $(BDIR)/detector_vs_particle $(BDIR)/particle_vs_detector

#$(SDIR)/dict.cxx :	$(SDIR)/ktTrackEff.hh
#	cd $(SDIR); rootcint6 -f dict.cxx -c -I. ./ktTrackEff.hh

$(ODIR)/functions.o : $(SDIR)/functions.cxx $(SDIR)/functions.hh
#$(ODIR)/dict.o : $(SDIR)/dict.cxx
#$(ODIR)/ktTrackEff.o : $(SDIR)/ktTrackEff.cxx $(SDIR)/ktTrackEff.hh
$(ODIR)/detector_vs_particle.o : $(SDIR)/detector_vs_particle.cxx
$(ODIR)/particle_vs_detector.o : $(SDIR)/particle_vs_detector.cxx

#data analysis
$(BDIR)/detector_vs_particle :	$(ODIR)/detector_vs_particle.o	$(ODIR)/functions.o
$(BDIR)/particle_vs_detector :	$(ODIR)/particle_vs_detector.o	$(ODIR)/functions.o

###############################################################################
##################################### MISC ####################################
###############################################################################

clean :
	@echo
	@echo CLEANING
	rm -vf $(ODIR)/*.o
	rm -vf $(BDIR)/*
	rm -vf lib/*

