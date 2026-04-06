CXX = g++
# CXXFLAGS = -g -O0 -std=c++17
CXXFLAGS = -O3 -ffast-math -march=native -std=c++17
# CXXFLAGS = -Wall -Wextra -Wcast-align -Wcast-qual -Wconversion -Wfloat-equal \
# 	    -Wformat=2 -Winit-self -Wmissing-declarations \
# 	    -Wmissing-include-dirs -Wpointer-arith -Wredundant-decls \
# 	    -Wswitch-default -Wuninitialized -Wwrite-strings \
# 	    -Wno-sign-conversion -Wno-unused-function \
#         -Wno-missing-declarations \
#         -std=c++14 -mcx16 -O3 -DNDEBUG
LDFLAGS = -lboost_system -lboost_filesystem -lboost_iostreams -ligraph

UNAME = $(shell uname)
ifeq ($(UNAME), Darwin)
	BREW_INCDIR = /opt/homebrew/include
	BREW_LIBDIR = /opt/homebrew/lib
	CXXFLAGS := -I$(BREW_INCDIR) $(CXXFLAGS)
	LDFLAGS := -L$(BREW_LIBDIR) -lboost_filesystem -lboost_iostreams -ligraph
endif

TARGET = abmc louvain gmc rabmc eblock

ASRCS := abmc.cpp
AOBJS := $(ASRCS:.cpp=.o)
AHDRS := $(ASRCS:.cpp=.hpp)

LSRCS := louvain.cpp
LOBJS := $(LSRCS:.cpp=.o)
LHDRS := $(LSRCS:.cpp=.hpp)

GSRCS := gmc.cpp
GOBJS := $(GSRCS:.cpp=.o)
GHDRS := $(GSRCS:.cpp=.hpp)

RSRCS := rabmc.cpp
ROBJS := $(RSRCS:.cpp=.o)

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

gmc: $(GOBJS) $(COBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

rabmc: $(ROBJS) $(COBJS)
	$(CXX) $(CXXFLAGS) -fopenmp -o $@ $^ $(LDFLAGS) -ltcmalloc_minimal -lnuma

%.o: %.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

rabmc.o: CXXFLAGS += -fopenmp

eblock: $(EOBJS) $(COBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

common/%.o: common/%.cpp $(CHDRS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

install: $(TARGET)
	install -d $(HOME)/.local/bin
	install -m 755 $(TARGET) $(HOME)/.local/bin/

clean:
	rm -f $(TARGET) $(AOBJS) $(LOBJS) $(GOBJS) $(ROBJS) $(COBJS)