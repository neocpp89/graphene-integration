# Project: Graphene Integration.
CXXFLAGS = -fno-omit-frame-pointer -g -O3 -march=native --std=c++11
.PHONY: clean
all: gengrid scattering unit


gengrid: GenerateCartesianGrid.o
	$(CXX) -o $@ $(CXXFLAGS) $^

scattering: Scattering.o
	$(CXX) -o $@ $(CXXFLAGS) $^

TESTS = TestSymmetryCoordinates.o

unit: catchmain.o $(TESTS)
	$(CXX) -o $@ $(CXXFLAGS) $^

clean:
	rm gengrid scattering \
	   	GenerateCartesianGrid.o \
	   	Scattering.o \
		catchmain.o $(TESTS)
