CXX := g++
CXXFLAGS := -O3 -march=native -std=c++17
LDFLAGS :=

UNAME = $(shell uname)
USER_NAME = $(shell whoami)
HOST_NAME = $(shell hostname)

ifeq ($(UNAME), Darwin)
	BREW_INCDIR = /opt/homebrew/include
	BREW_LIBDIR = /opt/homebrew/lib
	CXXFLAGS = -I$(BREW_INCDIR) $(CXXFLAGS)
	LDFLAGS = -L$(BREW_LIBDIR) -ligraph $(LDFLAGS)

	GVE_OMP_CXXFLAGS = -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include
	GVE_OMP_LDFLAGS  = -L/opt/homebrew/opt/libomp/lib -lomp
	RABBIT_EXTRA_LDFLAGS =
else ifeq ($(UNAME), Linux)
	ifeq ($(HOST_NAME), epyc)
		IGP_INCDIR = /opt/igraph-1.0.1/include
		IGP_LIBDIR = /opt/igraph-1.0.1/lib
		CXXFLAGS = $(CXXFLAGS) -I$(IGP_INCDIR) 
		LDFLAGS = $(LDFLAGS) -L$(IGP_LIBDIR) -ligraph -lopenblas
	endif

	GVE_OMP_CXXFLAGS = -fopenmp
	GVE_OMP_LDFLAGS  = -fopenmp
	RABBIT_EXTRA_LDFLAGS = -ltcmalloc_minimal -lnuma
endif

TARGET = gmc abmc louvain leiden leiden_cpm gve-leiden gve-leiden-cpm eblock ite_leiden ite_leiden_cpm
HDRS := common/Timer.hpp

ASRCS := abmc.cpp
AOBJS := $(ASRCS:.cpp=.o)

LSRCS := louvain.cpp
LOBJS := $(LSRCS:.cpp=.o)

L2SRCS := leiden.cpp
L2OBJS := $(L2SRCS:.cpp=.o)
L2CPMOBJS := leiden_cpm.o

GVESRCS := gve-leiden.cpp
GVEOBJS := $(GVESRCS:.cpp=.o)
GVEHDRS := $(wildcard gve-leiden-inc/*.hxx) common/Timer.hpp

CPMSRCS := gve-leiden-cpm.cpp
CPMOBJS := $(CPMSRCS:.cpp=.o)
CPMHDRS := $(wildcard gve-leiden-inc/*.hxx) $(wildcard inc/*.hxx) common/Types.hpp common/Timer.hpp

GSRCS := gmc.cpp
GOBJS := $(GSRCS:.cpp=.o)

RSRCS := rabmc.cpp
ROBJS := $(RSRCS:.cpp=.o)
RHDRS := $(RSRCS:.cpp=.hpp)

BSRCS := rabbit.cpp
BOBJS := $(BSRCS:.cpp=.o)
BHDRS := $(BSRCS:.cpp=.hpp)

MSRCS := rcm.cpp
MOBJS := $(MSRCS:.cpp=.o)
MHDRS := $(MSRCS:.cpp=.hpp)

CSRCS := common/mm_io.cpp common/Coloring.cpp common/BlockIO.cpp
COBJS := $(CSRCS:.cpp=.o)
CHDRS := $(CSRCS:.cpp=.hpp)

ESRCS := eblock.cpp
EOBJS := $(ESRCS:.cpp=.o)

ITESRCS := ite_leiden.cpp
ITEOBJS := $(ITESRCS:.cpp=.o)
ITECPMOBJS := ite_leiden_cpm.o

all: $(TARGET)

abmc: $(AOBJS) $(COBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

louvain: $(LOBJS) $(COBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

leiden: $(L2OBJS) $(COBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

leiden_cpm: $(L2CPMOBJS) $(COBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

gve-leiden: $(GVEOBJS) $(COBJS)
	$(CXX) $(CXXFLAGS) $(GVE_OMP_CXXFLAGS) -o $@ $^ $(LDFLAGS) $(GVE_OMP_LDFLAGS)

gve-leiden-cpm: $(CPMOBJS)
	$(CXX) $(CXXFLAGS) $(GVE_OMP_CXXFLAGS) -o $@ $^ $(GVE_OMP_LDFLAGS)

gmc: $(GOBJS) $(COBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

rabmc: $(ROBJS) $(COBJS)
	$(CXX) $(CXXFLAGS) $(GVE_OMP_CXXFLAGS) -o $@ $^ $(LDFLAGS) $(GVE_OMP_LDFLAGS) $(RABBIT_EXTRA_LDFLAGS)

rabbit: $(BOBJS) $(COBJS)
	$(CXX) $(CXXFLAGS) $(GVE_OMP_CXXFLAGS) -o $@ $^ $(LDFLAGS) $(GVE_OMP_LDFLAGS) $(RABBIT_EXTRA_LDFLAGS)

rcm: $(MOBJS) $(COBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

rabmc.o: CXXFLAGS += $(GVE_OMP_CXXFLAGS)
rabbit.o: CXXFLAGS += $(GVE_OMP_CXXFLAGS)

gve-leiden.o: gve-leiden.cpp $(GVEHDRS)
	$(CXX) $(CXXFLAGS) $(GVE_OMP_CXXFLAGS) -c $< -o $@

leiden_cpm.o: leiden.cpp common/Timer.hpp
	$(CXX) $(CXXFLAGS) -DCPM -c $< -o $@

gve-leiden-cpm.o: gve-leiden-cpm.cpp $(CPMHDRS)
	$(CXX) $(CXXFLAGS) -I. $(GVE_OMP_CXXFLAGS) -c $< -o $@

eblock: $(EOBJS) $(COBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

common/%.o: common/%.cpp $(CHDRS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

ite_leiden: $(ITEOBJS) $(COBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

ite_leiden_cpm: $(ITECPMOBJS) $(COBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

ite_leiden_cpm.o: ite_leiden.cpp common/Timer.hpp
	$(CXX) $(CXXFLAGS) -DCPM -c $< -o $@
	
clean:
	rm -f $(TARGET) $(AOBJS) $(LOBJS) $(L2OBJS) $(L2CPMOBJS) $(GOBJS) $(ROBJS) $(BOBJS) $(MOBJS) $(GVEOBJS) $(CPMOBJS) $(COBJS) $(EOBJS) $(ITEOBJS) $(ITECPMOBJS)
