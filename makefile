# Project: Graphene Integration.
CXXFLAGS = -fno-omit-frame-pointer -g -O3 -march=native --std=c++11 -Wall -Wextra
# CXXFLAGS = -fno-omit-frame-pointer -g -O0 --std=c++11 -Wall -Wextra
.PHONY: clean
all: gengrid scattering integrate unit


gengrid: GenerateCartesianGrid.o
	$(CXX) -o $@ $(CXXFLAGS) $^

scattering: Scattering.o
	$(CXX) -o $@ $(CXXFLAGS) $^

integrate: CalculateThermalConductivity.o
	$(CXX) -o $@ $(CXXFLAGS) $^

TESTS = TestSymmetryCoordinates.o

unit: catchmain.o $(TESTS)
	$(CXX) -o $@ $(CXXFLAGS) $^

clean:
	rm gengrid scattering integrate unit \
	   	GenerateCartesianGrid.o \
		CalculateThermalConductivity.o \
	   	Scattering.o \
		catchmain.o $(TESTS)
