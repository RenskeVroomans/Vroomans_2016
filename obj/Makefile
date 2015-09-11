#
# General C++ Makefile by Anton Crombach, A.B.M.Crombach@bio.uu.nl
#

#CXX = g++4
CXX = g++

#Flag to see whether we compile normal code or analysis code
ANA= 
THO= 
VID=  

CPPFLAGS =  -Wall -O2 -march=x86-64 -mmmx -msse -msse2 -fprefetch-loop-arrays -DDSFMT_MEXP=521 -DHAVE_SSE2
#CPPFLAGS = -g -Wall -march=x86-64  -DDSFMT_MEXP=521 #-DHAVE_SSE2
#CPPFLAGS = -pg -Wall -march=x86-64

#LNFLAGS = -pg -Wl,-rpath,/usr/lib/gcc-4.0.2/lib:/home/khwjtuss/local/lib
LNFLAGS = -Wl,-rpath,/usr/lib/gcc-4.0.2/lib:/home/renske/local/lib


MYPATH = /home/renske/local
INCDIR = -I../include -I$(MYPATH)/include -I/usr/include 
LIBDIR = -L$(MYPATH)/lib 

LIBS = -lpng -lboost_program_options -lfftw3 #-lefence 
# Source/object paths
vpath %.cc ../src ../test
vpath %.d .
vpath %.hh ../include
vpath %.o .

BINPATH = ../bin
LIBPATH = ../lib

# Project vars
ifndef ANA
PROJECT = ancestrace #babybugs	#
LIBRARY = mor
ALL = Misc.o Cell.o ChromBB.o Gene.o TFBS.o Genome.o Vertex.o Edge.o Network.o  Agent.o AgentAnalysis.o Population.o World.o dSFMT.o Vars.o IO.o
else ifndef THO
PROJECT = shortanalyse
LIBRARY = ana
ALL= Misc.o Cell.o ChromBB.o Gene.o TFBS.o Genome.o Vertex.o Edge.o Network.o  Agent.o AgentAnalysis.o Population.o LessOutputAnalyseOneAgent.o dSFMT.o Vars.o IO.o
else ifndef VID
PROJECT = fourier
LIBRARY = fou
ALL= Misc.o Cell.o ChromBB.o Gene.o TFBS.o Genome.o Vertex.o Edge.o Network.o  Agent.o AgentAnalysis.o Population.o FourierAnalysis.o dSFMT.o Vars.o IO.o
else
PROJECT = video
LIBRARY = loo
ALL= Misc.o Cell.o ChromBB.o Gene.o TFBS.o Genome.o Vertex.o Edge.o Network.o  Agent.o AgentAnalysis.o Population.o MakeVideo.o dSFMT.o Vars.o IO.o
endif

OBJECTS = $(ALL)

OBJECTSTWO=$(ALLTWO)

# Targets
all: $(PROJECT)

$(PROJECT): $(OBJECTS) 
	$(CXX) $(LNFLAGS) $(LIBDIR) $^ -o $(BINPATH)/$@ $(LIBS)

$(LIBRARY): $(OBJECTS)
	$(CXX) -shared -Wl,-soname,$@.so \
	$(LIBDIR) $(LIBS) $^ -o $(LIBPATH)/$@.so

$(OBJECTS): %.o: %.cc
	$(CXX) -c $(CPPFLAGS) $(INCDIR) $< -o $@

%.d: %.cc
	@set -e; rm -f $@; \
	$(CXX) $(INCDIR) -MM $(CPPFLAGS) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$


ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),realclean)
-include $(OBJECTS:.o=.d)
endif
endif

.PHONY: clean realclean 
clean:
	@rm -f *.d.* *.o *.d 

realclean: clean 
	@cd ../src; rm -f *~ .*.swp *.bak; \
	cd ../include; rm -f *~ .*.swp *.bak
