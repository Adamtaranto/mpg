CXXFLAGS    += -std=c++11 -Wall -O3 -g -Isrc
CXX         ?= g++

all: randseq mpg-burnin

%:src/%.cc src/libmpg.cc
	$(CXX) $(CXXFLAGS) -o $@ $^
