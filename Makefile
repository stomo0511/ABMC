CXX = g++
CXXFLAGS = -Wno-nan-infinity-disabled -O3 -march=native -ffast-math -std=c++17
LDFLAGS = -lboost_system -lboost_filesystem -lboost_iostreams -ligraph

UNAME = $(shell uname)
ifeq ($(UNAME), Darwin)
	BREW_INCDIR = /opt/homebrew/include
	BREW_LIBDIR = /opt/homebrew/lib
	CXXFLAGS := -I$(BREW_INCDIR) $(CXXFLAGS)
	LDFLAGS := -L$(BREW_LIBDIR) -lboost_filesystem -lboost_iostreams -ligraph
endif

TARGET = abmc louvain gmc rblock

ASRCS := abmc.cpp 
AOBJS := $(ASRCS:.cpp=.o)
AHDRS := $(ASRCS:.cpp=.hpp)

LSRCS := louvain.cpp
LOBJS := $(LSRCS:.cpp=.o)
LHDRS := $(LSRCS:.cpp=.hpp)

GSRCS := gmc.cpp
GOBJS := $(GSRCS:.cpp=.o)
GHDRS := $(GSRCS:.cpp=.hpp)

RSRCS := rblock.cpp 
ROBJS := $(RSRCS:.cpp=.o)
RHDRS := $(RSRCS:.cpp=.hpp)

CSRCS := common/mm_io.cpp common/Coloring.cpp common/BlockIO.cpp
COBJS := $(CSRCS:.cpp=.o)
CHDRS := $(CSRCS:.cpp=.hpp)


all: $(TARGET)

abmc: $(AOBJS) $(COBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

louvain: $(LOBJS) $(COBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

gmc: $(GOBJS) $(COBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

rblock: $(ROBJS) $(COBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

# %.o: %.cpp $(HDRS)
# 	$(CXX) $(CXXFLAGS) -c $<

common/%.o: common/%.cpp $(CHDRS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

install: $(TARGET)
	install -d $(HOME)/.local/bin
	install -m 755 $(TARGET) $(HOME)/.local/bin/

clean:
	rm -f $(TARGET) $(AOBJS) $(LOBJS) $(GOBJS)$(ROBJS) $(COBJS)