CXX = g++
# CXXFLAGS = -g -O0 -std=c++17
# CXXFLAGS = -O3 -ffast-math -march=native -std=c++17
CXXFLAGS = -O3 -march=native -std=c++17
# CXXFLAGS = -Wall -Wextra -Wcast-align -Wcast-qual -Wconversion -Wfloat-equal \
# 	    -Wformat=2 -Winit-self -Wmissing-declarations \
# 	    -Wmissing-include-dirs -Wpointer-arith -Wredundant-decls \
# 	    -Wswitch-default -Wuninitialized -Wwrite-strings \
# 	    -Wno-sign-conversion -Wno-unused-function \
#         -Wno-missing-declarations \
#         -std=c++14 -mcx16 -O3 -DNDEBUG
# LDFLAGS = -lboost_system -lboost_filesystem -lboost_iostreams -ligraph

UNAME = $(shell uname)
ifeq ($(UNAME), Darwin)
	BREW_INCDIR = /opt/homebrew/include
	BREW_LIBDIR = /opt/homebrew/lib
	CXXFLAGS := -I$(BREW_INCDIR) $(CXXFLAGS)
	LDFLAGS := -L$(BREW_LIBDIR) -lboost_filesystem -lboost_iostreams -ligraph
endif

UNAME = $(shell uname)
ifeq ($(UNAME), Linux)
	IGP_INCDIR = /opt/igraph-1.0.1/include
	IGP_LIBDIR = /opt/igraph-1.0.1/lib
	CXXFLAGS := -I$(IGP_INCDIR) $(CXXFLAGS)
# 	LDFLAGS := -L$(IGP_LIBDIR) -lboost_filesystem -lboost_iostreams -ligraph -lopenblas
	LDFLAGS := -L$(IGP_LIBDIR) -ligraph -lopenblas
endif

ifeq ($(UNAME), Darwin)
	GVE_OMP_CXXFLAGS = -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include
	GVE_OMP_LDFLAGS  = -L/opt/homebrew/opt/libomp/lib -lomp
	RABBIT_EXTRA_LDFLAGS =
else
	GVE_OMP_CXXFLAGS = -fopenmp
	GVE_OMP_LDFLAGS  = -fopenmp
	RABBIT_EXTRA_LDFLAGS = -ltcmalloc_minimal -lnuma
endif

TARGET = gmc abmc louvain leiden leiden_cpm gve-leiden gve-leiden-cpm eblock rcm
HDRS := common/Timer.hpp

ASRCS := abmc.cpp
AOBJS := $(ASRCS:.cpp=.o)
AHDRS := $(ASRCS:.cpp=.hpp)

LSRCS := louvain.cpp
LOBJS := $(LSRCS:.cpp=.o)
LHDRS := $(LSRCS:.cpp=.hpp)

L2SRCS := leiden.cpp
L2OBJS := $(L2SRCS:.cpp=.o)
L2HDRS := $(L2SRCS:.cpp=.hpp)
L2CPMOBJS := leiden_cpm.o

GVESRCS := gve-leiden.cpp
GVEOBJS := $(GVESRCS:.cpp=.o)
GVEHDRS := $(wildcard gve-leiden-inc/*.hxx) common/Timer.hpp

CPMSRCS := gve-leiden-cpm.cpp
CPMOBJS := $(CPMSRCS:.cpp=.o)
CPMHDRS := $(wildcard gve-leiden-inc/*.hxx) $(wildcard inc/*.hxx) common/Types.hpp common/Timer.hpp

GSRCS := gmc.cpp
GOBJS := $(GSRCS:.cpp=.o)
GHDRS := $(GSRCS:.cpp=.hpp)

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

install: $(TARGET)
	install -d $(HOME)/.local/bin
	install -m 755 $(TARGET) $(HOME)/.local/bin/

clean:
	rm -f $(TARGET) $(AOBJS) $(LOBJS) $(L2OBJS) $(L2CPMOBJS) $(GOBJS) $(ROBJS) $(BOBJS) $(MOBJS) $(GVEOBJS) $(CPMOBJS) $(COBJS) $(EOBJS)
