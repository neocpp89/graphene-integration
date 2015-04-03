# Project: Graphene Integration.
CXXFLAGS = -fno-omit-frame-pointer -g -O3 -march=native --std=c++11
.PHONY: clean

all: gengrid scattering

tests: testsympoints testinterpolations

clean:
	rm gengrid scattering

gengrid: GenerateCartesianGrid.cpp
	$(CXX) -o gengrid $(CXXFLAGS) GenerateCartesianGrid.cpp

scattering: Scattering.cpp
	$(CXX) -o scattering $(CXXFLAGS) Scattering.cpp
